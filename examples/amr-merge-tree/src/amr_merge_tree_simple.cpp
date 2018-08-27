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

#include "../../local-global/output_persistence.h"


// block-independent types
using Bounds = diy::DiscreteBounds;
using AmrVertexId = r::AmrVertexId;
using AmrEdge = reeber::AmrEdge;
using AmrEdgeContainer = reeber::AmrEdgeContainer;

#define DIM 3

using FabBlockR = FabBlock<Real, DIM>;

using Block = FabTmtBlock<Real, DIM>;
using Vertex = Block::Vertex;
using Component = Block::Component;
using VertexNeighborMap = Block::TripletMergeTree::VertexNeighborMap;
using AmrTripletMergeTree = Block::TripletMergeTree;
using MaskedBox = Block::MaskedBox;

using Neighbor = AmrTripletMergeTree::Neighbor;


struct IsAmrVertexLocal
{
    bool operator()(const Block& b, const Neighbor& from) const
    {
        return from->vertex.gid == b.gid;
    }
};

using OutputPairsR = OutputPairs<Block, IsAmrVertexLocal>;

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
bool link_contains_gid(Link* link, int gid)
{
    for (int i = 0; i < link->size(); ++i)
        if (link->target(i).gid == gid)
            return true;
    return false;
}

/**
 *
 * send outgoing edges computed in b to all its neighbors
 * used to symmetrize the edge set in the beginning of the algorithm
 *
 * @param b FabTmtBlock
 * block that will send its edges
 *
 * @param cp diy::Master::ProxyWithLink
 * communication proxy
 *
 */
template<unsigned D>
void send_edges_to_neighbors(FabTmtBlock<double, D>* b, const diy::Master::ProxyWithLink& cp)
{
    bool debug = false;

    if (debug) fmt::print("Called send_edges_to_neighbors for block = {}\n", b->gid);

    auto* l = static_cast<diy::AMRLink*>(cp.link());

    std::set<diy::BlockID> receivers;
    for (int i = 0; i < l->size(); ++i) {
        if (l->target(i).gid != b->gid) {
            receivers.insert(l->target(i));
        }
    }

    for (const diy::BlockID& receiver : receivers) {
        int receiver_gid = receiver.gid;
        if (b->new_receivers_.count(receiver_gid)) {
            if (debug) {
                fmt::print("In send_edges_to_neighbors for block = {}, sending to receiver= {}, cardinality = {}\n",
                           b->gid, receiver_gid, b->get_all_outgoing_edges().size());
            };
            cp.enqueue(receiver, b->get_all_outgoing_edges());
        } else {
            if (debug)
                fmt::print("In send_edges_to_neighbors for block = {}, sending to receiver= {} empty container\n",
                           b->gid, receiver_gid);
            cp.enqueue(receiver, r::AmrEdgeContainer());
        }
    }
}

/**
 *
 * @tparam D
 * @param b
 * @param cp
 */
template<unsigned D>
void delete_low_edges(FabTmtBlock<double, D>* b, const diy::Master::ProxyWithLink& cp)
{
    bool debug = false;

    if (debug) fmt::print("Called delete_low_edges for block = {}\n", b->gid);

    auto* l = static_cast<diy::AMRLink*>(cp.link());

    std::set<diy::BlockID> senders;
    for (int i = 0; i < l->size(); ++i) {
        if (l->target(i).gid != b->gid) {
            senders.insert(l->target(i));
        }
    }

    for (const diy::BlockID& sender : senders) {
        AmrEdgeContainer edges_from_neighbor;
        if (debug) fmt::print("In delete_low_edges for block = {}, dequeing from sender = {}\n", b->gid, sender.gid);

        cp.dequeue(sender, edges_from_neighbor);

        if (debug)
            fmt::print("In delete_low_edges for block = {}, dequed {} edges from sender = {}\n", b->gid,
                       edges_from_neighbor.size(), sender.gid);

        b->delete_low_edges(sender.gid, edges_from_neighbor);

        if (debug)
            fmt::print("In delete_low_edges for block = {}, from sender = {}, b->delete_low_edges OK\n", b->gid,
                       sender.gid);
    }

    b->adjust_outgoing_edges();
}


template<unsigned D>
void send_to_neighbors(FabTmtBlock<double, D>* b, const diy::Master::ProxyWithLink& cp)
{
    bool debug = false;
    //    if (debug) fmt::print("Called send_to_neighbors for block = {}\n", b->gid);
    auto* l = static_cast<diy::AMRLink*>(cp.link());

    cp.collectives()->clear();

    std::set<diy::BlockID> receivers;
    for (int i = 0; i < l->size(); ++i) {
        if (l->target(i).gid != b->gid) {
            receivers.insert(l->target(i));
        }
    }

    if (debug)
        fmt::print("In send_to_neighbors for block = {}, link size = {}, unique = {}\n", b->gid, l->size(),
                   receivers.size());


    for (const diy::BlockID& receiver : receivers) {
        int receiver_gid = receiver.gid;


        // if we have sent our tree to this receiver before, only send n_trees = 0
        // else send the tree and all outgoing edges
        int n_trees = (b->processed_receiveres_.count(receiver_gid) == 0);

        if (debug)
            fmt::print("In send_to_neighbors for block = {}, sending to {}, n_trees = {}\n", b->gid, receiver_gid,
                       n_trees);

        cp.enqueue(receiver, n_trees);

        if (n_trees) {
            // send local tree and all outgoing edges that end in receiver
            cp.enqueue(receiver, b->original_tree_);
            cp.enqueue(receiver, b->vertex_to_deepest_);
            cp.enqueue(receiver, b->get_deepest_vertices());
            cp.enqueue(receiver, b->get_all_outgoing_edges());
            cp.enqueue(receiver, b->get_original_link_gids());

            if (debug)
                fmt::print("In send_to_neighbors for block = {}, receiver = {}, enqueued original_link_gids = {}\n",
                           b->gid, receiver.gid, container_to_string(b->get_original_link_gids()));

            diy::MemoryBuffer& out = cp.outgoing(receiver);
            diy::LinkFactory::save(out, l);
        }
        // mark receiver_gid as processed
        b->new_receivers_.erase(receiver_gid);
        b->processed_receiveres_.insert(receiver_gid);
    }

    int done = b->done_;
    b->round_++;
    if (debug)
        fmt::print("In send_to_neighbors for block = {}, b->done = {}, b->round = {}\n", b->gid, done, b->round_);
    cp.all_reduce(done, std::logical_and<int>());
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
expand_link(Block* b, const diy::Master::ProxyWithLink& cp, diy::AMRLink* l, std::vector<diy::AMRLink>& received_links,
            std::vector<std::vector<int>>& received_original_gids)
{
    bool debug = false;
    if (debug) fmt::print("in expand_link for block = {}, started updating link\n", b->gid);
    int n_added = 0;
    assert(received_links.size() == received_original_gids.size());
    std::set<int> added_gids;

    for (size_t i = 0; i < received_links.size(); ++i) {
        const diy::AMRLink& received_link = received_links[i];
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
                "In expand_link for block = {}, b->done_ = {}, n_added = {}, new link size = {}, new link size_unqie = {}, added_gids = {}\n",
                b->gid, b->done_, n_added, l->size(), l->size_unique(), container_to_string(added_gids));
    cp.master()->add_expected(n_added);
}

/**
 *
 * @tparam D dimension, unsigned int template parameter
 * @param b FabTmtBlock
 * @param cp Communication proxy
 */

template<unsigned D>
void get_from_neighbors_and_merge(FabTmtBlock<double, D>* b, const diy::Master::ProxyWithLink& cp)
{
    bool debug = false;

    //    if (debug) fmt::print("Called get_from_neighbors_and_merge for block = {}\n", b->gid);

    using Block = FabTmtBlock<double, D>;
    using AmrTripletMergeTree = typename Block::TripletMergeTree;
    using AmrEdgeVector = typename Block::AmrEdgeContainer;
    using VertexVertexMap = typename Block::VertexVertexMap;
    //    using VertexSizeMap = typename Block::VertexSizeMap;
    using LinkVector = std::vector<diy::AMRLink>;


    auto* l = static_cast<diy::AMRLink*>(cp.link());

    std::set<diy::BlockID> senders;
    for (int i = 0; i < l->size(); ++i) {
        if (l->target(i).gid != b->gid) {
            senders.insert(l->target(i));
        }
    }

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
            diy::AMRLink* l = static_cast<diy::AMRLink*>(diy::LinkFactory::load(in));
            received_links.push_back(*l);
            delete l;
        }
    }

    assert(received_trees.size() == received_vertex_to_deepest.size() and
           received_trees.size() == received_edges.size() and
           received_trees.size() == received_deepest_vertices.size() and
           received_trees.size() == received_original_gids.size());

#ifdef DEBUG
    // just for debug - check all received edges
    for (const AmrEdgeContainer& ec : received_edges) {
        for (const AmrEdge& e : ec) {
            AmrVertexId my_vertex = std::get<1>(e);

            assert(my_vertex.gid != b->gid or
                   b->local_.mask_by_index(my_vertex) == MaskedBox::ACTIVE or
                   b->local_.mask_by_index(my_vertex) == MaskedBox::LOW);
        }
    }
#endif

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
            fmt::print("In get_from_neighbors_and_merge for block = {}, merge and repair OK for sender = {}\n", b->gid,
                       sender_gids_debug[i]);

        // save information about vertex-component relation and component merging in block
        b->vertex_to_deepest_.insert(received_vertex_to_deepest[i].begin(), received_vertex_to_deepest[i].end());

        // make_set in disjoint-sets of components
        for (const AmrVertexId& new_deepest_vertex : received_deepest_vertices[i]) {
            b->add_component_to_disjoint_sets(new_deepest_vertex);
        }

        if (debug)
            fmt::print("In get_from_neighbors_and_merge for block = {}, add_component_to_disjoint_sets OK\n", b->gid);
        // update receivers - if a block from the 1-neighbourhood of sender has not received a tree from us
        // we must send it to this block in the next round

        for (const diy::AMRLink& rl : received_links) {
            for (int k = 0; k < rl.size(); ++k) {
                int sender_neighbor_gid = rl.target(k).gid;
                if (b->processed_receiveres_.count(sender_neighbor_gid) == 0)
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
                assert(b->edge_goes_out(e));
                vertices_to_check.push_back(std::get<0>(e));
            }
        }
    }

    if (debug) fmt::print("In get_from_neighbors_and_merge for block = {}, disjoint sets updated OK\n", b->gid);

    b->done_ = b->is_done_simple(vertices_to_check);

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


/**
 *
 * @tparam D dimension, template parameter
 * @param b FabTmtBlock
 * @param cp communication proxy
 * @param first bool, true, if it is the first round
 *
 * if it is not the first round, then receive the information from the previous round
 * and update local merge trees.
 * in every round send to neighbors.
 */
template<unsigned D>
void amr_merge_tree_simple(FabTmtBlock<double, D>* b, const diy::Master::ProxyWithLink& cp, bool first)
{
    if (!first) {
        // receive from neighbors
        get_from_neighbors_and_merge<D>(b, cp);
    }
    // in all rounds send to neighbors
    send_to_neighbors<D>(b, cp);
}

void read_from_file(std::string infn,
                    diy::mpi::communicator& world,
                    diy::Master& master_reader,
                    diy::Master& master,
                    diy::ContiguousAssigner& assigner,
                    diy::MemoryBuffer& header,
                    diy::DiscreteBounds& domain)
{
    diy::io::read_blocks(infn, world, assigner, master_reader, header, FabBlockR::load);
    diy::load(header, domain);
    fmt::print("data read\n");
}

int main(int argc, char** argv)
{
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
    double rho = 1E54;

    using namespace opts;

    opts::Options ops(argc, argv);
    ops
            >> Option('b', "blocks", nblocks, "number of blocks to use")
            >> Option('m', "memory", in_memory, "maximum blocks to store in memory")
            >> Option('j', "jobs", threads, "threads to use during the computation")
            >> Option('s', "storage", prefix, "storage prefix")
            >> Option('t', "threshold", rho, "threshold")
            >> Option('p', "profile", profile_path, "path to keep the execution profile")
            >> Option('l', "log", log_level, "log level");

    bool negate = ops >> opts::Present('n', "negate", "sweep superlevel sets");
    bool split = ops >> Present("split", "use split IO");

    std::string infn, outfn, outdiagfn;

    if (ops >> Present('h', "help", "show help message") or
        not(ops >> PosOption(infn))
        or not(ops >> PosOption(outfn))) {
        if (world.rank() == 0) {
            fmt::print("Usage: {} INPUT.AMR OUTPUT \n", argv[0]);
            fmt::print("Compute local-global tree from AMR data\n");
            fmt::print("{}", ops);
        }
        return 1;
    }

    bool write_diag = (ops >> PosOption(outdiagfn));

    diy::FileStorage storage(prefix);

    diy::Master master_reader(world, 1, -1, FabBlockR::create);
    diy::Master master(world, threads, in_memory, &Block::create, &Block::destroy, &storage, &Block::save, &Block::load);
    diy::ContiguousAssigner assigner(world.size(), 0);
    diy::MemoryBuffer header;
    diy::DiscreteBounds domain;

    world.barrier();
    dlog::Timer timer;
    LOG_SEV_IF(world.rank() == 0, info) << "Starting computation";

    read_from_file(infn, world, master_reader, master, assigner, header, domain);

    world.barrier();

    LOG_SEV_IF(world.rank() == 0, info) << "Data read, local size = " << master.size();
    LOG_SEV_IF(world.rank() == 0, info) << "Time to read data:       " << dlog::clock_to_string(timer.elapsed());
    timer.restart();


    // copy FabBlocks to FabTmtBlocks
    // in FabTmtConstructor mask will be set and local trees will be computed
    // FabBlock can be safely discarded afterwards

    master_reader.foreach(
            [&master, &assigner, domain, rho, negate](FabBlockR* b, const diy::Master::ProxyWithLink& cp) {
                auto* l = static_cast<diy::AMRLink*>(cp.link());
                diy::AMRLink* new_link = new diy::AMRLink(*l);

                // prepare neighbor box info to save in MaskedBox
                int local_ref = l->refinement();
                int local_lev = l->level();


                fmt::print("copying block, bounds = {} - {}, core = {} - {}, gid = {}\n", l->bounds().min,
                           l->bounds().max, l->core().min, l->core().max, cp.gid());

                master.add(cp.gid(),
                           new Block(b->fab, local_ref, local_lev, domain, l->bounds(), l->core(), cp.gid(),
                                     new_link, rho, negate),
                           new_link);
            });

    fmt::print("FabBlocks copied\n");

    world.barrier();
    LOG_SEV_IF(world.rank() == 0, info) << "Time to compute local trees and components:  " << dlog::clock_to_string(timer.elapsed());
    timer.restart();

    int global_done = false;
    int rounds = 0;
    bool first = true;
    while (true) {
        rounds++;
        if (rounds == 1) {
            // remove non-existing edges that end in LOW vertex in neighbouring block
            master.foreach(&send_edges_to_neighbors<DIM>);
            master.exchange();
            master.foreach(&delete_low_edges<3>);
        } else {

            master.foreach(
                    [&](Block* b, const diy::Master::ProxyWithLink& cp) { amr_merge_tree_simple<DIM>(b, cp, first); });
            master.exchange();

            // get done flag from any available block, they are all the same

            global_done = master.proxy(master.loaded_block()).read<int>();

            if (master.communicator().rank() == 0) { fmt::print("MASTER round {}, global_done = {}\n", rounds, global_done); }

            if (global_done)
                break;
            first = false;
        }
    }

    world.barrier();
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
    if (outfn != "none") {
        if (!split)
            diy::io::write_blocks(outfn, world, master);
        else
            diy::io::split::write_blocks(outfn, world, master);
    }

    world.barrier();
    LOG_SEV_IF(world.rank() == 0, info) << "Time to write tree:  " << dlog::clock_to_string(timer.elapsed());
    timer.restart();

    bool verbose = false;

    if (write_diag) {
        OutputPairsR::ExtraInfo extra(outdiagfn, verbose);
        IsAmrVertexLocal test_local;
        master.foreach([&extra, &test_local](Block* b, const diy::Master::ProxyWithLink& cp) {
            output_persistence(b, cp, extra, test_local);
        });
    }

    world.barrier();
    LOG_SEV_IF(world.rank() == 0, info) << "Time to write diagrams:  " << dlog::clock_to_string(timer.elapsed());
    timer.restart();

    return 0;
}
