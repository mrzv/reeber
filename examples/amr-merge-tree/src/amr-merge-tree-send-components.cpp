#include "reeber-real.h"

#define SEND_COMPONENTS

#include <diy/master.hpp>
#include <diy/io/block.hpp>
#include <opts/opts.h>

#include <reeber/box.h>

#include "fab-block.h"
#include "fab-tmt-block.h"
#include "amr-merge-tree-helper.h"
#include "reader-interfaces.h"
#include "diy/vertices.hpp"
#include "reeber/grid.h"


// block-independent types
using Bounds = diy::DiscreteBounds;
using AmrVertexId = r::AmrVertexId;
using AmrEdge = reeber::AmrEdge;
using AmrEdgeContainer = reeber::AmrEdgeContainer;


using FabBlockR = FabBlock<Real, 3>;

using Block = FabTmtBlock<Real, 3>;
using Vertex = Block::Vertex;
using AmrVertexContainer = Block::AmrVertexContainer;
using GidContainer = Block::GidContainer;
using Component = Block::Component;
using VertexNeighborMap = Block::TripletMergeTree::VertexNeighborMap;
using AmrTripletMergeTree = Block::TripletMergeTree;
using MaskedBox = Block::MaskedBox;

/**
 *
 * @param link *AMRLink
 * Link in which the block is searched
 *
 * @param gid int
 * gid of the block
 *
 * @return BlockID of block with given gid from link, throw, if not found
 */
template<class Link>
diy::BlockID get_block_id_by_gid(Link* link, int gid)
{
    for (int i = 0; i < link->size(); ++i)
        if (link->target(i).gid == gid) {
            return link->target(i);
        }
    throw std::runtime_error("Cannot find block with gid = " + std::to_string(gid));
}

std::string get_link_gids(diy::AMRLink* l)
{
    std::set<int> gids;
    for (int i = 0; i < l->size(); ++i) {
        gids.insert(l->target(i).gid);
    }
    return container_to_string(gids);
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
bool link_contains_gid(Link* link, int gid)
{
    for (int i = 0; i < link->size(); ++i)
        if (link->target(i).gid == gid)
            return true;
    return false;
}

void
expand_link(Block* b, const diy::Master::ProxyWithLink& cp, diy::AMRLink* l, std::vector<diy::AMRLink>& received_links)
{
    bool debug = true;
    if (debug) fmt::print("In expand_link for block = {}, started updating link", b->gid);
    int n_added = 0;
    for (size_t i = 0; i < received_links.size(); ++i) {
        const diy::AMRLink& received_link = received_links[i];
        for (int j = 0; j < received_link.size(); ++j) {

            int candidate_gid = received_link.target(j).gid;

            if (link_contains_gid(l, candidate_gid))
                continue;

            if (not b->gid_must_be_in_link(candidate_gid))
                continue;

            n_added++;
            l->add_neighbor(received_link.target(j));
            l->add_bounds(received_link.level(j), received_link.refinement(j), received_link.core(j),
                          received_link.bounds(j));
            //                l->add_wrap(received_link.wrap(j));
        }
    }
    if (debug)
        fmt::print(
                "In expand_link for block = {}, b->done_ = {}, n_added = {}, new link size = {}, new link size_unqie = {}\n",
                b->gid, b->done_, n_added, l->size(), l->size_unique());
    cp.master()->add_expected(n_added);

}

template<unsigned D>
void send_to_neighbors_main(FabTmtBlock<Real, D>* b, const diy::Master::ProxyWithLink& cp)
{
    bool debug = true;
    if (debug) fmt::print("Called send_to_neighbors_main for block = {}\n", b->gid);

    auto* l = static_cast<diy::AMRLink*>(cp.link());

    cp.collectives()->clear();

    // get unique receivers
    std::set<diy::BlockID> receivers;
    for (int i = 0; i < l->size(); ++i) {
        if (l->target(i).gid != b->gid) {
            receivers.insert(l->target(i));
        }
    }

    std::map<int, int> componenents_per_gid;

    for (const Component& c : b->components_) {
        for (int component_receiver : c.current_neighbors_) {
            componenents_per_gid[component_receiver]++;
        }
    }

    //    if (debug) fmt::print("In send_to_neighbors for block = {}, link size = {}, unique = {}\n", b->gid, l->size(), receivers.size());
    for (const diy::BlockID& receiver : receivers) {
        int receiver_gid = receiver.gid;

        // if we have sent our tree to this receiver before, only send n_trees = 0
        // else send the tree and all outgoing edges
        int n_trees = componenents_per_gid[receiver_gid];

        cp.enqueue(receiver, n_trees);

        if (n_trees) {
            // send components that must communicate with the receiver
            for (Component& c : b->components_) {
                if (c.current_neighbors_.count(receiver_gid) == 1 and c.processed_neighbors_.count(receiver_gid) == 0) {
                    cp.enqueue(receiver, c.merge_tree_);
                    cp.enqueue(receiver, c.root_);
                    cp.enqueue(receiver, c.outgoing_edges_);
                    cp.enqueue(receiver, c.current_neighbors_);
                    c.processed_neighbors_.insert(receiver_gid);
                }
            }

            cp.enqueue(receiver, b->original_vertex_to_deepest_);

            diy::MemoryBuffer& out = cp.outgoing(receiver);
            diy::LinkFactory::save(out, l);
        }
    }

    int done = b->are_all_components_done();
    b->round_++;
    if (debug) fmt::print("In send_to_neighbors for block = {}, done = {}, b->round = {}\n", b->gid, done, b->round_);
    cp.all_reduce(done, std::logical_and<int>());
}

template<unsigned D>
void receive_main(FabTmtBlock<Real, D>* b, const diy::Master::ProxyWithLink& cp)
{
    bool debug = true;
    //    if (debug) fmt::print("Called receive_main for block = {}\n", b->gid);
    using Block = FabTmtBlock<Real, D>;
    using AmrTripletMergeTree = typename Block::TripletMergeTree;
    using AmrEdgeVector = typename Block::AmrEdgeContainer;
    using VertexVertexMap = typename Block::VertexVertexMap;
    using LinkVector = std::vector<diy::AMRLink>;


    auto* l = static_cast<diy::AMRLink*>(cp.link());

    std::set<diy::BlockID> senders;
    for (int i = 0; i < l->size(); ++i) {
        if (l->target(i).gid != b->gid) {
            senders.insert(l->target(i));
        }
    }

    AmrEdgeContainer all_received_edges;
    LinkVector received_links;

    if (debug) fmt::print("In receive_main for block = {}, # senders = {}\n", b->gid, senders.size());

    std::map<AmrVertexId, std::set<int>> component_to_neighbors;

    for (const diy::BlockID& sender : senders) {
        int n_trees;
        cp.dequeue(sender, n_trees);
        if (n_trees > 0) {
            for (int tree_idx = 0; tree_idx < n_trees; ++tree_idx) {
                AmrTripletMergeTree received_tree;
                cp.dequeue(sender, received_tree);

                AmrVertexId received_root;
                cp.dequeue(sender, received_root);

                AmrEdgeVector received_edges;
                cp.dequeue(sender, received_edges);
                all_received_edges.insert(all_received_edges.end(), received_edges.begin(), received_edges.end());

                GidContainer received_current_gids;
                cp.dequeue(sender, received_current_gids);

                // just for debug - check all received edges
                for (const AmrEdge& e : received_edges) {
                    AmrVertexId my_vertex = std::get<1>(e);
                    assert(my_vertex.gid != b->gid or
                           b->local_.mask_by_index(my_vertex) == MaskedBox::ACTIVE or
                           b->local_.mask_by_index(my_vertex) == MaskedBox::LOW);
                }

                // merge received trees
                r::merge(b->mt_, received_tree, received_edges, true);
                r::repair(b->mt_);
                if (debug) fmt::print("In receive_main for block = {}, repair OK\n", b->gid);

                // make_set in disjoint-sets of components
                b->add_component_to_disjoint_sets(received_root);

                // add gids from sender
                component_to_neighbors[received_root].insert(received_current_gids.begin(),
                                                             received_current_gids.end());


            } // loop over all trees received from current sender

            VertexVertexMap received_vertex_to_deepest;
            cp.dequeue(sender, received_vertex_to_deepest);
            // save information about vertex-component relation and component merging in block
            b->original_vertex_to_deepest_.insert(received_vertex_to_deepest.begin(),
                                         received_vertex_to_deepest.end());

            diy::MemoryBuffer& in = cp.incoming(sender.gid);
            diy::AMRLink* l = static_cast<diy::AMRLink*>(diy::LinkFactory::load(in));
            received_links.push_back(*l);
            delete l;
        }
    }

    // all tree from all senders were processed, now we can connect components in the disjoint sets data structure
    for (const AmrEdge& e : all_received_edges) {
        if (b->edge_exists(e)) {
            // edge e connects two vertices that we have, connect their components
            AmrVertexId deepest_a = b->original_deepest(std::get<0>(e));
            AmrVertexId deepest_b = b->original_deepest(std::get<1>(e));
            b->connect_components(deepest_a, deepest_b);
        } else {
            assert(b->edge_goes_out(e));
        }
    }
    all_received_edges.clear();

    // expand current neighbourhood of each component
    for (Component& c : b->components_) {
        for (const auto& sent_component_gids_pair : component_to_neighbors) {
            AmrVertexId sent_root = sent_component_gids_pair.first;
            const std::set<int>& sent_gids = sent_component_gids_pair.second;
            if (b->are_components_connected(c.root_, sent_root)) {
                c.current_neighbors_.insert(sent_gids.begin(), sent_gids.end());
            }
        }
    }

    expand_link(b, cp, l, received_links);
}

void read_from_file(std::string infn,
                    diy::mpi::communicator& world,
                    diy::Master& master_reader,
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

    bool absolute = ops >> Present('a', "absolute", "use absolute values for thresholds (instead of multiples of mean)");
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


    diy::Master master_reader(world, 1, -1, FabBlockR::create);
    diy::Master master(world, 1, -1);
    diy::ContiguousAssigner assigner(world.size(), 0);
    diy::MemoryBuffer header;
    diy::DiscreteBounds domain;

    read_from_file(input_filename, world, master_reader, assigner, header, domain);

    // copy FabBlocks to FabTmtBlocks
    // in FabTmtConstructor mask will be set and local trees will be computed
    // FabBlock can be safely discarded afterwards

    master_reader.foreach(
            [&master, domain, rho, negate, absolute](FabBlockR* b, const diy::Master::ProxyWithLink& cp) {
                auto* l = static_cast<diy::AMRLink*>(cp.link());
                diy::AMRLink* new_link = new diy::AMRLink(*l);

                // prepare neighbor box info to save in MaskedBox
                int local_ref = l->refinement();
                int local_lev = l->level();

                master.add(cp.gid(),
                           new Block(b->fab, local_ref, local_lev, domain, l->bounds(), l->core(), cp.gid(),
                                     new_link, rho, negate, absolute),
                           new_link);
            });

    fmt::print("FabBlocks copied\n");

    int global_n_undone = 1;
    int rounds = 0;

    // symmetrize edges
    master.foreach(&send_edges_to_neighbors<3>);
    master.exchange();
    master.foreach(&delete_low_edges<3>);

    while (true) {
        rounds++;
        master.foreach(&send_to_neighbors_main<3>);
        master.exchange();
        master.foreach(&receive_main<3>);
        master.exchange();
        global_n_undone = master.proxy(master.loaded_block()).read<int>();
        if (0 == global_n_undone)
            break;
        if (master.communicator().rank() == 0) {
            fmt::print("MASTER round {}, global_n_undone = {}\n", rounds, global_n_undone);
        }
    }
    return 0;
}
