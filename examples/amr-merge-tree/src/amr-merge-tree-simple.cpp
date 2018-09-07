#include "reeber-real.h"

// to print nice backtrace on segfault signal
#include <signal.h>
#include <execinfo.h>
#include <cxxabi.h>

#include <diy/master.hpp>
#include <diy/io/block.hpp>
#include <diy/decomposition.hpp>

#include <dlog/stats.h>
#include <dlog/log.h>
#include <opts/opts.h>

#include <reeber/box.h>

#include "fab-block.h"
#include "fab-tmt-block.h"
#include "reader-interfaces.h"
#include "diy/vertices.hpp"
#include "reeber/grid.h"
#include "amr-merge-tree-helper.h"

#include "../../local-global/output-persistence.h"

#include "read-npy.h"


// block-independent types
using AMRLink = diy::AMRLink;

using Bounds = diy::DiscreteBounds;
using AmrVertexId = r::AmrVertexId;
using AmrEdge = reeber::AmrEdge;


#define DIM 3

using FabBlockR = FabBlock<Real, DIM>;

using Block = FabTmtBlock<Real, DIM>;
using Vertex = Block::Vertex;
using Component = Block::Component;
using VertexNeighborMap = Block::TripletMergeTree::VertexNeighborMap;
using AmrTripletMergeTree = Block::TripletMergeTree;
using MaskedBox = Block::MaskedBox;
using GidVector = Block::GidVector;

using Neighbor = AmrTripletMergeTree::Neighbor;

constexpr bool abort_on_segfault = true;

struct IsAmrVertexLocal
{
    bool operator()(const Block& b, const Neighbor& from) const
    {
        return from->vertex.gid == b.gid;
    }
};

using OutputPairsR = OutputPairs<Block, IsAmrVertexLocal>;


std::set<diy::BlockID> link_unique(AMRLink *amr_link, int gid)
{
    std::set<diy::BlockID> result;
    for (int i = 0; i < amr_link->size(); ++i) {
        if (amr_link->target(i).gid != gid) {
            result.insert(amr_link->target(i));
        }
    }
    return result;
}

/**
 *
 * @param link *AMRLink
 * Link in which the block is searched
 *
 * @param gid int
 * gid of the block
 *
 * @return true, if link contains a block with given gid
 */
template<class Link>
bool link_contains_gid(Link *link, int gid)
{
    for (int i = 0; i < link->size(); ++i)
        if (link->target(i).gid == gid)
            return true;
    return false;
}


template<unsigned D>
void send_to_neighbors(FabTmtBlock<Real, D> *b, const diy::Master::ProxyWithLink& cp)
{
//    bool debug = (b->gid == 3) || (b->gid == 11) || (b->gid == 0) || (b->gid == 1);
    bool debug = false;
    if (debug) fmt::print("Called send_to_neighbors for block = {}\n", b->gid);

    auto *l = static_cast<AMRLink *>(cp.link());


    auto receivers = link_unique(l, b->gid);

    if (debug)
        fmt::print("In send_to_neighbors for block = {}, link size = {}, unique = {}\n", b->gid, l->size(),
                   receivers.size());


    for (const diy::BlockID& receiver : receivers) {
        int receiver_gid = receiver.gid;


        // if we have sent our tree to this receiver before, only send n_trees = 0
        // else send the tree and all outgoing edges
        int n_trees = (b->processed_receivers_.count(receiver_gid) == 0 and
                       b->new_receivers_.count(receiver_gid) == 1);

        if (debug)
            fmt::print("In send_to_neighbors for block = {}, sending to {}, n_trees = {}\n", b->gid, receiver_gid,
                       n_trees);

        cp.enqueue(receiver, n_trees);

        if (n_trees) {
            // send local tree and all outgoing edges that end in receiver
            cp.enqueue(receiver, b->original_tree_);
            cp.enqueue(receiver, b->original_vertex_to_deepest_);
            cp.enqueue(receiver, b->get_original_deepest_vertices());
            cp.enqueue(receiver, b->get_all_outgoing_edges());
            cp.enqueue(receiver, b->get_original_link_gids());

            if (debug)
                fmt::print("In send_to_neighbors for block = {}, receiver = {}, enqueued original_link_gids = {}\n",
                           b->gid, receiver.gid, container_to_string(b->get_original_link_gids()));

            diy::MemoryBuffer& out = cp.outgoing(receiver);
            diy::LinkFactory::save(out, l);

            // mark receiver_gid as processed
            b->new_receivers_.erase(receiver_gid);
            b->processed_receivers_.insert(receiver_gid);
        }
    }

    int done = b->done_;
    b->round_++;
    if (debug)
        fmt::print("In send_to_neighbors for block = {}, b->done = {}, b->round = {}\n", b->gid, done, b->round_);
}

/**
 *
 * @param b FabTmtBlock, used only to print debug information
 * @param cp communication proxy
 * @param l link of b, to be updated
 * @param received_links links received
 * @param received_original_gids
 */
void
expand_link(Block *b, const diy::Master::ProxyWithLink& cp, AMRLink *l, std::vector<AMRLink>& received_links,
            std::vector<std::vector<int>>& received_original_gids)
{
//    bool debug = (b->gid == 3) || (b->gid == 11) || (b->gid == 0) || (b->gid == 1);
    bool debug = false;
    if (debug) fmt::print("in expand_link for block = {}, round = {}, started updating link\n", b->gid, b->round_);
    int n_added = 0;
    assert(received_links.size() == received_original_gids.size());
    std::set<int> added_gids;

    for (size_t i = 0; i < received_links.size(); ++i) {
        const AMRLink& received_link = received_links[i];
        assert(not received_original_gids[i].empty());
        if (debug)
            fmt::print("in expand_link for block = {}, i = {}, received_links[i].size = {}\n", b->gid, i,
                       received_links[i].size());

        for (int j = 0; j < received_link.size(); ++j) {
            // if we are already sending to this block, skip it
            int candidate_gid = received_link.target(j).gid;
            if (link_contains_gid(l, candidate_gid))
                continue;

            if (debug)
                fmt::print("in expand_link for block = {}, candidate_gid = {}, received_original_links[{}] = {}\n",
                           b->gid,
                           candidate_gid, i, container_to_string(received_original_gids[i]));

            // skip non-original gids (we only include the original link)
            if (std::find(received_original_gids[i].begin(), received_original_gids[i].end(), candidate_gid) ==
                received_original_gids[i].end()) {
                if (debug)
                    fmt::print("in expand_link for block = {}, gid = {} not in original gids, skipping\n", b->gid,
                               candidate_gid);
                continue;
            }

            n_added++;
            l->add_neighbor(received_link.target(j));
            added_gids.insert(candidate_gid);
            l->add_bounds(received_link.level(j), received_link.refinement(j), received_link.core(j),
                          received_link.bounds(j));
            if (debug) fmt::print("in expand_link for block = {}, added gid = {}\n", b->gid, candidate_gid);
        }
    }

    if (debug)
        fmt::print(
                "In expand_link for block = {}, round = {}, b->done_ = {}, n_added = {}, new link size = {}, new link size_unqie = {}, added_gids = {}\n",
                b->gid, b->round_, b->done_, n_added, l->size(), l->size_unique(), container_to_string(added_gids));
    cp.master()->add_expected(n_added);
}

/**
 *
 * @tparam D dimension, unsigned int template parameter
 * @param b FabTmtBlock
 * @param cp Communication proxy
 */

template<unsigned D>
void get_from_neighbors_and_merge(FabTmtBlock<Real, D> *b, const diy::Master::ProxyWithLink& cp)
{
//    bool debug = (b->gid == 3) || (b->gid == 11) || (b->gid == 0) || (b->gid == 1);
    bool debug = false;

    //    if (debug) fmt::print("Called get_from_neighbors_and_merge for block = {}\n", b->gid);

    using Block = FabTmtBlock<Real, D>;
    using AmrTripletMergeTree = typename Block::TripletMergeTree;
    using AmrEdgeVector = typename Block::AmrEdgeContainer;
    using VertexVertexMap = typename Block::VertexVertexMap;
    //    using VertexSizeMap = typename Block::VertexSizeMap;
    using LinkVector = std::vector<AMRLink>;


    auto *l = static_cast<AMRLink *>(cp.link());

    cp.collectives()->clear();

    auto senders = link_unique(l, b->gid);

    // TODO: delete this
    std::vector<int> sender_gids_debug;

    // TODO: delete this
    std::vector<int> debug_received_ntrees;

    std::vector<AmrTripletMergeTree> received_trees;
    std::vector<VertexVertexMap> received_vertex_to_deepest;
    std::vector<AmrEdgeVector> received_edges;
    std::vector<std::vector<AmrVertexId>> received_deepest_vertices;
    std::vector<std::vector<int>> received_original_gids;

    LinkVector received_links;

    if (debug) fmt::print("In get_from_neighbors_and_merge for block = {}, # senders = {}\n", b->gid, senders.size());

    for (const diy::BlockID& sender : senders) {
        int n_trees;
        cp.dequeue(sender, n_trees);

        debug_received_ntrees.push_back(n_trees);

        if (debug)
            fmt::print("In get_from_neighbors_and_merge for block = {}, dequeued from sender {} n_trees = {} \n",
                       b->gid, sender.gid, n_trees);

        if (n_trees > 0) {
            assert(n_trees == 1);
            sender_gids_debug.push_back(sender.gid);

            received_trees.emplace_back();
            received_vertex_to_deepest.emplace_back();
            received_deepest_vertices.emplace_back();
            received_edges.emplace_back();
            received_original_gids.emplace_back();

            cp.dequeue(sender, received_trees.back());
            cp.dequeue(sender, received_vertex_to_deepest.back());
            cp.dequeue(sender, received_deepest_vertices.back());
            cp.dequeue(sender, received_edges.back());
            cp.dequeue(sender, received_original_gids.back());

            if (debug)
                fmt::print(
                        "In get_from_neighbors_and_merge for block = {}, dequeued from sender {} original link gids = {}\n",
                        b->gid, sender.gid, container_to_string(received_original_gids.back()));

            diy::MemoryBuffer& in = cp.incoming(sender.gid);
            AMRLink *l = static_cast<AMRLink *>(diy::LinkFactory::load(in));
            received_links.push_back(*l);
            delete l;
        }
    }

    assert(received_trees.size() == received_vertex_to_deepest.size() and
           received_trees.size() == received_edges.size() and
           received_trees.size() == received_deepest_vertices.size() and
           received_trees.size() == received_original_gids.size());

    //#ifdef DEBUG
    //    // just for debug - check all received edges
    //    for (const AmrEdgeContainer& ec : received_edges)
    //    {
    //        for (const AmrEdge& e : ec)
    //        {
    //            AmrVertexId my_vertex = std::get<1>(e);

    //            assert(my_vertex.gid != b->gid or
    //                   b->local_.mask_by_index(my_vertex) == MaskedBox::ACTIVE or
    //                   b->local_.mask_by_index(my_vertex) == MaskedBox::LOW);
    //        }
    //    }
    //#endif

    //if (debug) fmt::print("In get_from_neighbors_and_merge for block = {}, dequeed all, edges checked OK\n", b->gid);

    // vertices in our processed neighbourhood, from which there is an edge going out to blocks we have not communicated with
    // if any of these vertices is in a connected component of original tree, we are not done
    std::vector<AmrVertexId> vertices_to_check;

    // merge all received trees
    for (size_t i = 0; i < received_trees.size(); ++i) {
        if (received_edges[i].empty())
            continue;
        AmrTripletMergeTree& rt = received_trees[i];
        r::merge(b->mt_, rt, received_edges[i], true);
        r::repair(b->mt_);
        if (debug)
            fmt::print(
                    "In get_from_neighbors_and_merge for block = {}, merge and repair OK for sender = {}, tree size = {}\n",
                    b->gid,
                    sender_gids_debug[i], b->mt_.size());

        // save information about vertex-component relation and component merging in block
        b->original_vertex_to_deepest_.insert(received_vertex_to_deepest[i].begin(), received_vertex_to_deepest[i].end());

        // make_set in disjoint-sets of components
        for (const AmrVertexId& new_deepest_vertex : received_deepest_vertices[i]) {
            b->add_component_to_disjoint_sets(new_deepest_vertex);
        }

        if (debug)
            fmt::print("In get_from_neighbors_and_merge for block = {}, add_component_to_disjoint_sets OK\n", b->gid);
        // update receivers - if a block from the 1-neighbourhood of sender has not received a tree from us
        // we must send it to this block in the next round

        for (const AMRLink& rl : received_links) {
            for (int k = 0; k < rl.size(); ++k) {
                int sender_neighbor_gid = rl.target(k).gid;

                auto& original_gids = received_original_gids[i];

                bool is_in_original_gids = std::find(original_gids.begin(), original_gids.end(), sender_neighbor_gid) !=
                                           original_gids.end();
                bool is_not_processed = b->processed_receivers_.count(sender_neighbor_gid) == 0;

                if (debug)
                    fmt::print(
                            "In get_from_neighbors_and_merge for block = {}, round = {}, sender_neighbor_gid = {}, is_in_original_gids = {}, is_not_processed = {}\n",
                            b->gid, b->round_, sender_neighbor_gid, is_in_original_gids, is_not_processed);

                if (is_not_processed and is_in_original_gids)
                    b->new_receivers_.insert(sender_neighbor_gid);
            }
        }
        if (debug) fmt::print("In get_from_neighbors_and_merge for block = {}, processed_receiveres_ OK\n", b->gid);
    }

    // update disjoint sets data structure (some components are now connected to each other)
    for (size_t i = 0; i < received_trees.size(); ++i) {
        for (const AmrEdge& e : received_edges[i]) {
            if (b->edge_exists(e)) {
                // edge e connects two vertices that we have, connect their components
                AmrVertexId deepest_a = b->deepest(std::get<0>(e));
                AmrVertexId deepest_b = b->deepest(std::get<1>(e));
                b->connect_components(deepest_a, deepest_b);
            } else {
                vertices_to_check.push_back(std::get<0>(e));
            }
        }
    }

    if (debug)
        fmt::print("In get_from_neighbors_and_merge for block = {}, disjoint sets updated OK, tree size = {}\n", b->gid,
                   b->mt_.size());

    b->done_ = b->is_done_simple(vertices_to_check);
    int n_undone = 1 - b->done_;

    cp.all_reduce(n_undone, std::plus<int>());

    if (debug)
        fmt::print("In get_from_neighbors_and_merge for block = {}, is_done_simple OK, vertices_to_check.size = {}\n",
                   b->gid, vertices_to_check.size());

    int old_size_unique = l->size_unique();
    int old_size = l->size();

    if (debug)
        fmt::print(
                "In get_from_neighbors_and_merge for block = {}, b->done_ = {}, old link size = {}, old link size_unqie = {}\n",
                b->gid, b->done_, old_size, old_size_unique);

    expand_link(b, cp, l, received_links, received_original_gids);

    if (debug) fmt::print("In get_from_neighbors_and_merge for block = {}, expand_link OK\n", b->gid);
}


inline bool ends_with(const std::string& s, const std::string& suffix)
{
    if (suffix.size() > s.size())
        return false;
    return std::equal(suffix.rbegin(), suffix.rend(), s.rbegin());
}

void read_from_file(std::string infn,
                    diy::mpi::communicator& world,
                    diy::Master& master_reader,
                    diy::ContiguousAssigner& assigner,
                    diy::MemoryBuffer& header,
                    diy::DiscreteBounds& domain,
                    int nblocks)
{
    if (ends_with(infn, ".npy")) {
//        fmt::print("read npy\n");
        read_from_npy_file<DIM>(infn, world, nblocks, master_reader, assigner, header, domain);
    } else {
//        fmt::print("read amr\n");
        diy::io::read_blocks(infn, world, assigner, master_reader, header, FabBlockR::load);
        diy::load(header, domain);
    }
}


void catch_sig(int signum)
{
    LOG_SEV_IF(true, fatal) << "caught signal " << signum;
    //    << ",  local group = " << pro {}, local rank = {}", signum, active_puppet, proc_map->group(), proc_map->local_rank());

    // print backtrace
    void *callstack[128];
    int frames = backtrace(callstack, 128);
    char **strs = backtrace_symbols(callstack, frames);

    size_t funcnamesize = 256;
    char *funcname = (char *) malloc(funcnamesize);

    // iterate over the returned symbol lines. skip the first, it is the
    // address of this function.
    for (int i = 1; i < frames; i++) {
        char *begin_name = 0, *begin_offset = 0, *end_offset = 0;

        // find parentheses and +address offset surrounding the mangled name:
        // ./module(function+0x15c) [0x8048a6d]
        for (char *p = strs[i]; *p; ++p) {
            if (*p == '(')
                begin_name = p;
            else if (*p == '+')
                begin_offset = p;
            else if (*p == ')' && begin_offset) {
                end_offset = p;
                break;
            }
        }

        if (begin_name && begin_offset && end_offset && begin_name < begin_offset) {
            *begin_name++ = '\0';
            *begin_offset++ = '\0';
            *end_offset = '\0';

            // mangled name is now in [begin_name, begin_offset) and caller
            // offset in [begin_offset, end_offset). now apply __cxa_demangle():

            int status;
            char *ret = abi::__cxa_demangle(begin_name, funcname, &funcnamesize, &status);
            if (status == 0) {
                funcname = ret; // use possibly realloc()-ed string
                LOG_SEV_IF(true, fatal) << "  " << strs[i] << " : " << funcname << "+" << begin_offset;
            } else {
                // demangling failed. Output function name as a C function with no arguments.
                //                logger->critical("  {} : {}()+{}", strs[i], begin_name, begin_offset);
                LOG_SEV_IF(true, fatal)  << "  " << strs[i] << " : " << funcname << "+" << begin_offset;
            }
        } else {
            // couldn't parse the line? print the whole line.
            LOG_SEV_IF(true, fatal)  << "  " << strs[i];
            //            logger->critical("  {}", strs[i]);
        }
    }

    free(funcname);
    free(strs);

    signal(signum, SIG_DFL);    // restore the default signal
    if (abort_on_segfault)
        MPI_Abort(MPI_COMM_WORLD, 1);
}


int main(int argc, char **argv)
{
    signal(SIGSEGV, catch_sig);
    signal(SIGABRT, catch_sig);

    diy::mpi::environment env(argc, argv);
    diy::mpi::communicator world;

    if (argc < 2) {
        fmt::print(std::cerr, "Usage: {} IN.amr rho\n", argv[0]);
        return 1;
    }

    int nblocks = world.size();
    std::string prefix = "./DIY.XXXXXX";
    int in_memory = -1;
    int threads = 1;
    std::string profile_path;
    std::string log_level = "info";

    // threshold
    Real rho = 81.66;
    Real integral_rho = 90.0;

    using namespace opts;

    opts::Options ops(argc, argv);
    ops
            >> Option('b', "blocks", nblocks, "number of blocks to use")
            >> Option('m', "memory", in_memory, "maximum blocks to store in memory")
            >> Option('j', "jobs", threads, "threads to use during the computation")
            >> Option('s', "storage", prefix, "storage prefix")
            >> Option('t', "threshold", rho, "threshold")
            >> Option("intthreshold", integral_rho, "integral threshold")
            >> Option('p', "profile", profile_path, "path to keep the execution profile")
            >> Option('l', "log", log_level, "log level");

    bool absolute =
            ops >> Present('a', "absolute", "use absolute values for thresholds (instead of multiples of mean)");
    bool negate = ops >> opts::Present('n', "negate", "sweep superlevel sets");
    bool split = ops >> Present("split", "use split IO");

    std::string input_filename, output_filename, output_diagrams_filename, output_integral_filename;

    if (ops >> Present('h', "help", "show help message") or
        not(ops >> PosOption(input_filename))
        or not(ops >> PosOption(output_filename))) {
        if (world.rank() == 0) {
            fmt::print("Usage: {} INPUT.AMR OUTPUT \n", argv[0]);
            fmt::print("Compute local-global tree from AMR data\n");
            fmt::print("{}", ops);
        }
        return 1;
    }

    bool write_diag = (ops >> PosOption(output_diagrams_filename));

    bool write_integral = (ops >> PosOption(output_integral_filename));

    if (write_integral) {
        if ((negate and integral_rho < rho) or (not negate and integral_rho > rho))
            throw std::runtime_error("Bad integral threshold");
    }

    diy::FileStorage storage(prefix);

    diy::Master master_reader(world, 1, in_memory, &FabBlockR::create, &FabBlockR::destroy);
    diy::Master master(world, threads, in_memory, &Block::create, &Block::destroy, &storage, &Block::save,
                       &Block::load);
    diy::ContiguousAssigner assigner(world.size(), nblocks);
    diy::MemoryBuffer header;
    diy::DiscreteBounds domain;

    dlog::add_stream(std::cerr, dlog::severity(log_level))
            << dlog::stamp() << dlog::aux_reporter(world.rank()) << dlog::color_pre() << dlog::level()
            << dlog::color_post() >> dlog::flush();

    world.barrier();
    dlog::Timer timer;
    LOG_SEV_IF(world.rank() == 0, info) << "Starting computation, input_filename = " << input_filename << ", nblocks = "
                                                                                     << nblocks
                                                                                     << ", rho = " << rho;
    world.barrier();

    read_from_file(input_filename, world, master_reader, assigner, header, domain, nblocks);

    world.barrier();

    LOG_SEV_IF(world.rank() == 0, info) << "Data read, local size = " << master.size();
    LOG_SEV_IF(world.rank() == 0, info) << "Time to read data:       " << dlog::clock_to_string(timer.elapsed());
    timer.restart();

    world.barrier();

    // copy FabBlocks to FabTmtBlocks
    // in FabTmtConstructor mask will be set and local trees will be computed
    // FabBlock can be safely discarded afterwards

    master_reader.foreach(
            [&master, &assigner, domain, rho, negate, absolute](FabBlockR *b, const diy::Master::ProxyWithLink& cp) {
                auto *l = static_cast<AMRLink *>(cp.link());
                AMRLink *new_link = new AMRLink(*l);

                // prepare neighbor box info to save in MaskedBox
                int local_ref = l->refinement();
                int local_lev = l->level();

                master.add(cp.gid(),
                           new Block(b->fab, local_ref, local_lev, domain, l->bounds(), l->core(), cp.gid(),
                                     new_link, rho, negate, absolute),
                           new_link);

            });


    if (absolute) {
        LOG_SEV_IF(world.rank() == 0, info) << "Time to compute local trees and components:  "
                << dlog::clock_to_string(timer.elapsed());
        timer.restart();
    } else {
        LOG_SEV_IF(world.rank() == 0, info) << "Time to construct FabTmtBlocks: "
                << dlog::clock_to_string(timer.elapsed());
        timer.restart();

        master.foreach([&rho](Block *b, const diy::Master::ProxyWithLink& cp) {
            cp.collectives()->clear();
            cp.all_reduce(b->sum_, std::plus<Real>());
            cp.all_reduce(static_cast<Real>(b->n_unmasked_) * b->scaling_factor(), std::plus<Real>());
        });

        master.exchange();

        const diy::Master::ProxyWithLink& proxy = master.proxy(master.loaded_block());

        Real mean = proxy.get<Real>() / proxy.get<Real>();
        rho *= mean;                                            // now rho contains absolute threshold
        integral_rho *= mean;

        world.barrier();
        LOG_SEV_IF(world.rank() == 0, info) << "Average = " << mean << ", rho = " << rho
                                                            << ", time to compute average: "
                                                            << dlog::clock_to_string(timer.elapsed());
        timer.restart();

        master.foreach([rho](Block *b, const diy::Master::ProxyWithLink& cp) {
            AMRLink *l = static_cast<AMRLink *>(cp.link());
            b->init(rho, l);
            cp.collectives()->clear();
        });


        world.barrier();
        LOG_SEV_IF(world.rank() == 0, info) <<
        "Time to initializee FabTmtBlocks (low vertices, local trees, components, outgoing edges): "
                << dlog::clock_to_string(timer.elapsed());
        timer.restart();
    }

    int global_n_undone = 1;

    master.foreach(&send_edges_to_neighbors<DIM>);
    master.exchange();
    master.foreach(&delete_low_edges<DIM>);

    LOG_SEV_IF(world.rank() == 0, info)  << "edges symmetrized, time elapsed "
            << dlog::clock_to_string(timer.elapsed());
    timer.restart();

    master.foreach([](Block *b, const diy::Master::ProxyWithLink& cp) {
        auto *l = static_cast<AMRLink *>(cp.link());

        std::set<diy::BlockID> receivers;
        for (int i = 0; i < l->size(); ++i) {
            if (l->target(i).gid != b->gid) {
                receivers.insert(l->target(i));
            }
        }

        for (const diy::BlockID& receiver : receivers) {
            int receiver_gid = receiver.gid;
            cp.enqueue(receiver, (int) (b->new_receivers_.count(receiver_gid)));
        }
    });

    master.exchange();

    master.foreach([](Block *b, const diy::Master::ProxyWithLink& cp) {
        auto *l = static_cast<AMRLink *>(cp.link());
        std::set<diy::BlockID> senders;
        for (int i = 0; i < l->size(); ++i) {
            if (l->target(i).gid != b->gid) {
                senders.insert(l->target(i));
            }
        }

        for (const diy::BlockID& sender : senders) {
            int t;
            cp.dequeue(sender, t);

            if (t != (int) (b->new_receivers_.count(sender.gid))) {
                throw std::runtime_error("Asymmetry");
            }

        }
    });

    LOG_SEV_IF(world.rank() == 0, info)  << "Symmetry checked in "
            << dlog::clock_to_string(timer.elapsed());
    timer.restart();


    int rounds = 0;
    while (global_n_undone) {
        rounds++;
        master.foreach(&send_to_neighbors<DIM>);
        master.exchange();
        master.foreach(&get_from_neighbors_and_merge<DIM>);
        master.exchange();
        global_n_undone = master.proxy(master.loaded_block()).read<int>();
        LOG_SEV_IF(world.rank() == 0, info) << "MASTER round " << rounds << ", global_n_undone = " << global_n_undone;
    }

    world.barrier();

    //    fmt::print("world.rank = {}, time for exchange = {}\n", world.rank(), dlog::clock_to_string(timer.elapsed()));

    LOG_SEV_IF(world.rank() == 0, info) << "Time for exchange:  " << dlog::clock_to_string(timer.elapsed());
    timer.restart();

    //    fmt::print("----------------------------------------\n");
    //    master.foreach([](Block* b, const diy::Master::ProxyWithLink& cp) {
    //        int my_nodes = 0;
    //        int not_my_nodes = 0;
    //        const auto& const_tree = b->mt;
    //        for(const auto& x : const_tree.nodes()) {
    //            if (x.second->vertex != x.first)
    //                fmt::print("Bad node vertex = {}, node vertex = {}\n", x.first, x.second->vertex);
    //            if (b->gid == x.second->vertex.gid) {
    //                my_nodes++;
    //            } else {
    //                not_my_nodes++;
    //            }
    //        }
    //        fmt::print("finished gid = {}, mt.size = {}, my_nodes = {}, not_my_nodes = {}\n", b->gid, b->mt.size(), my_nodes, not_my_nodes);
    //    });
    //    fmt::print("----------------------------------------\n");

    // save the result
    if (output_filename != "none") {
        if (!split)
            diy::io::write_blocks(output_filename, world, master);
        else
            diy::io::split::write_blocks(output_filename, world, master);
    }

    world.barrier();
    LOG_SEV_IF(world.rank() == 0, info) << "Time to write tree:  " << dlog::clock_to_string(timer.elapsed());
    timer.restart();

    bool verbose = false;

    if (write_diag) {
        bool ignore_zero_persistence = true;
        OutputPairsR::ExtraInfo extra(output_diagrams_filename, verbose);
        IsAmrVertexLocal test_local;
        master.foreach(
                [&extra, &test_local, ignore_zero_persistence, rho](Block *b, const diy::Master::ProxyWithLink& cp) {
                    output_persistence(b, cp, extra, test_local, rho, ignore_zero_persistence);
                });
    }

    world.barrier();
    LOG_SEV_IF(world.rank() == 0, info) << "Time to write diagrams:  " << dlog::clock_to_string(timer.elapsed());
    timer.restart();

    if (write_integral) {
        Real rho_min, rho_max;
        if (negate) {
            rho_min = rho;
            rho_max = integral_rho;
        } else {
            rho_min = integral_rho;
            rho_max = rho;
        }

        master.foreach([rho_min, rho_max, output_integral_filename](Block *b, const diy::Master::ProxyWithLink& cp) {

            b->compute_final_connected_components();

            b->compute_local_integral(rho_min, rho_max);

            AMRLink *l = static_cast<AMRLink *>(cp.link());

            for(const auto& vertex_value_pair : b->local_integral_)
            {
                int receiver_gid = vertex_value_pair.first.gid;
                if (receiver_gid == b->gid)
                    continue;
                auto receiver = l->target(l->find(receiver_gid));
                cp.enqueue(receiver, vertex_value_pair);
            }
        });

        master.exchange();

        master.foreach([rho_min, rho_max, output_integral_filename](Block *b, const diy::Master::ProxyWithLink& cp) {
            AMRLink *l = static_cast<AMRLink *>(cp.link());
            for(auto bid : l->neighbors())
            {
                while(cp.incoming(bid.gid))
                {
                    Block::LocalIntegral::value_type x;
                    cp.dequeue(bid, x);
                    assert(x.first.gid == b->gid);
                    b->local_integral_[x.first] += x.second;
                }
            }

            std::string integral_local_fname = fmt::format("{}-b{}.integral", output_integral_filename, b->gid);
            std::ofstream ofs(integral_local_fname);
            for (const auto& root_value_pair : b->local_integral_) {
                AmrVertexId root = root_value_pair.first;
                if (root.gid != b->gid)
                    continue;
                fmt::print(ofs, "{} {} {}\n", b->local_.global_position(root), root_value_pair.second, root);
        }
        });

        world.barrier();
        LOG_SEV_IF(world.rank() == 0, info) << "Time to compute and write integral:  " << dlog::clock_to_string(timer.elapsed());
        timer.restart();
    }

    master.foreach([](Block *b, const diy::Master::ProxyWithLink& cp) {
        auto sum_n_vertices_pair = b->get_local_stats();
        cp.collectives()->clear();
        cp.all_reduce(sum_n_vertices_pair.first, std::plus<Real>());
        cp.all_reduce(sum_n_vertices_pair.second, std::plus<size_t>());

    });

    master.exchange();

    world.barrier();

    const diy::Master::ProxyWithLink& proxy = master.proxy(master.loaded_block());

    Real total_value = proxy.get<Real>();
    size_t total_vertices = proxy.get<size_t>();

    LOG_SEV_IF(world.rank() == 0, info) << "Total value = " << total_value << ", total # vertices = " << total_vertices;

    return 0;
}
