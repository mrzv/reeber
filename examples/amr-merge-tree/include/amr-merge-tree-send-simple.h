#pragma once

#include <vector>

#include <diy/link.hpp>
#include <reeber/box.h>
#include "fab-block.h"
#include "fab-tmt-block.h"
#include "amr-merge-tree-helper.h"
using AMRLink = diy::AMRLink;

template<class Real, unsigned D>
void expand_link(FabTmtBlock<Real, D>* b,
                 const diy::Master::ProxyWithLink& cp,
                 AMRLink* l,
                 std::vector<AMRLink>& received_links,
                 std::vector<std::vector<int>>& received_original_gids)
{
//    bool debug = (b->gid == 3) || (b->gid == 11) || (b->gid == 0) || (b->gid == 1);
    bool debug = false;
    if (debug) fmt::print("in expand_link for block = {}, round = {}, started updating link\n", b->gid, b->round_);
    int n_added = 0;
    assert(received_links.size() == received_original_gids.size());
    std::set<int> added_gids;

    for(size_t i = 0; i < received_links.size(); ++i)
    {
        const AMRLink& received_link = received_links[i];
        assert(not received_original_gids[i].empty());
        //if (debug) fmt::print("in expand_link for block = {}, i = {}, received_links[i].size = {}\n", b->gid, i, received_links[i].size());

        for(int j = 0; j < received_link.size(); ++j)
        {
            // if we are already sending to this block, skip it
            int candidate_gid = received_link.target(j).gid;
            if (link_contains_gid(l, candidate_gid))
                continue;

            //if (debug) fmt::print("in expand_link for block = {}, candidate_gid = {}, received_original_links[{}] = {}\n", b->gid, candidate_gid, i, container_to_string(received_original_gids[i]));
            // skip non-original gids (we only include the original link)
            if (std::find(received_original_gids[i].begin(), received_original_gids[i].end(), candidate_gid) ==
                received_original_gids[i].end()) {
                //if (debug) fmt::print("in expand_link for block = {}, gid = {} not in original gids, skipping\n", b->gid, candidate_gid);
                continue;
            }
            n_added++;
            l->add_neighbor(received_link.target(j));
            added_gids.insert(candidate_gid);
            l->add_bounds(received_link.level(j), received_link.refinement(j), received_link.core(j),
                          received_link.bounds(j));
            //if (debug) fmt::print("in expand_link for block = {}, added gid = {}\n", b->gid, candidate_gid);
        }
    }

    //if (debug) fmt::print( "In expand_link for block = {}, round = {}, b->done_ = {}, n_added = {}, new link size = {}, new link size_unqie = {}, added_gids = {}\n", b->gid, b->round_, b->done_, n_added, l->size(), l->size_unique(), container_to_string(added_gids));
    if (debug)
        fmt::print(
                "In expand_link for block = {}, round = {}, b->done_ = {}, n_added = {}, new link size = {}, new link size_unqie = {}\n",
                b->gid, b->round_, b->done_, n_added, l->size(), l->size_unique());
    cp.master()->add_expected(n_added);
}

template<class Real, unsigned D>
void amr_tmt_send(FabTmtBlock<Real, D>* b, const diy::Master::ProxyWithLink& cp)
{

//    bool debug = (b->gid == 1 or b->gid == 100);
    bool debug = false;
    if (debug) fmt::print("Called send_simple for block = {}\n", b->gid);

    auto* l = static_cast<AMRLink*>(cp.link());


    auto receivers = link_unique(l, b->gid);

//    if (debug) fmt::print("In send_simple for block = {}, link size = {}, unique = {}\n", b->gid, l->size(), receivers.size());


    for (const diy::BlockID& receiver : receivers)
    {
        int receiver_gid = receiver.gid;

        cp.enqueue(receiver, b->get_original_link_gids());

        //if (debug) fmt::print("In send_simple for block = {}, receiver = {}, enqueued original_link_gids = {}\n", b->gid, receiver.gid, container_to_string(b->get_original_link_gids()));

        diy::MemoryBuffer& out = cp.outgoing(receiver);
        diy::LinkFactory::save(out, l);

        // if we have sent our tree to this receiver before, only send n_trees = 0
        // else send the tree and all outgoing edges
        int n_trees = (b->processed_receivers_.count(receiver_gid) == 0 and
                       b->new_receivers_.count(receiver_gid) == 1
                       and b->done_ == 0);


//        if (debug) fmt::print("In send_simple for block = {}, sending to {}, n_trees = {}\n", b->gid, receiver_gid, n_trees);

        cp.enqueue(receiver, n_trees);

        if (n_trees)
        {
            // send local tree and all outgoing edges that end in receiver
            cp.enqueue(receiver, b->original_tree_);
            cp.enqueue(receiver, b->original_vertex_to_deepest_);
            cp.enqueue(receiver, b->get_original_deepest_vertices());
            cp.enqueue(receiver, b->get_all_outgoing_edges());

            if (debug) fmt::print("In send_simple for block = {}, sent data to {}, n_trees = {}, tree_size = {}, deepest_vertices = {}, n_edges = {}\n",
                    b->gid, receiver_gid, n_trees, b->original_tree_.size(), b->get_original_deepest_vertices().size(), b->get_all_outgoing_edges().size());

            // mark receiver_gid as processed
            b->new_receivers_.erase(receiver_gid);
            b->processed_receivers_.insert(receiver_gid);
        }
    }

    int done = b->done_;
    b->round_++;
    if (debug) fmt::print("Exit send_simple for block = {}, b->done = {}, b->round = {}\n", b->gid, done, b->round_);
}

template<class Real, unsigned D>
void amr_tmt_receive(FabTmtBlock<Real, D>* b, const diy::Master::ProxyWithLink& cp)
{
#ifdef DO_DETAILED_TIMING
    dlog::Timer timer;
    dlog::Timer rl_loop_timer;
    dlog::Timer merge_timer;
    dlog::Timer union_find_timer;

    // detailed timings
    using DurationType = decltype(timer.elapsed());

//    DurationType time_to_receive_trees_and_gids;
//    DurationType rl_loop_time;
//    DurationType repair_time;
//    DurationType merge_call_time;
//    DurationType union_find_time;
//    DurationType sparsify_time;
//    DurationType expand_link_time;

#endif

    bool debug = false; // (b->gid == 1 or b->gid == 100);

    if (debug) fmt::print("Called receive_simple for block = {}\n", b->gid);

    using Block = FabTmtBlock<Real, D>;
    using AmrTripletMergeTree = typename Block::TripletMergeTree;
    using AmrEdgeVector = typename Block::AmrEdgeContainer;
    using VertexVertexMap = typename Block::VertexVertexMap;
    //    using VertexSizeMap = typename Block::VertexSizeMap;
    using LinkVector = std::vector<AMRLink>;


    auto* l = static_cast<AMRLink*>(cp.link());


    std::vector<AmrTripletMergeTree> received_trees;
    std::vector<VertexVertexMap> received_vertex_to_deepest;
    std::vector<AmrEdgeVector> received_edges;
    std::vector<std::vector<AmrVertexId>> received_deepest_vertices;
    std::vector<std::vector<int>> received_original_gids;

    LinkVector received_links;

#ifdef DO_DETAILED_TIMING
    timer.restart();
#endif

    auto senders = link_unique(l, b->gid);
    std::vector<int> sender_gids_debug;
    for(const auto& sender : senders)
    {
        sender_gids_debug.push_back(sender.gid);
    }

    if (debug) fmt::print("In receive_simple for block = {}, # senders = {}\n", b->gid, senders.size());

    for (const diy::BlockID& sender : senders)
    {
        int n_trees;

        received_original_gids.emplace_back();
        cp.dequeue(sender, received_original_gids.back());
        diy::MemoryBuffer& in = cp.incoming(sender.gid);
        AMRLink* l = static_cast<AMRLink*>(diy::LinkFactory::load(in));
        received_links.push_back(*l);
        delete l;

        cp.dequeue(sender, n_trees);

//        if (debug) fmt::print("In receive_simple for block = {}, dequeued from sender {} n_trees = {} \n", b->gid, sender.gid, n_trees);

        if (n_trees > 0)
        {
            assert(n_trees == 1);

            received_trees.emplace_back();
            received_vertex_to_deepest.emplace_back();
            received_deepest_vertices.emplace_back();
            received_edges.emplace_back();

            cp.dequeue(sender, received_trees.back());
            cp.dequeue(sender, received_vertex_to_deepest.back());
            cp.dequeue(sender, received_deepest_vertices.back());
            cp.dequeue(sender, received_edges.back());

//            if (debug) fmt::print( "In receive_simple for block = {}, dequeued from sender {} original link gids = {}\n", b->gid, sender.gid, container_to_string(received_original_gids.back()));

        }
    }

#ifdef DO_DETAILED_TIMING
    b->receive_trees_and_gids_time += timer.elapsed();
    timer.restart();
#endif


    assert(received_trees.size() == received_vertex_to_deepest.size() and
           received_trees.size() == received_edges.size() and
           received_trees.size() == received_deepest_vertices.size());

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

    //if (debug) fmt::print("In receive_simple for block = {}, dequeed all, edges checked OK\n", b->gid);

    // vertices in our processed neighbourhood, from which there is an edge going out to blocks we have not communicated with
    // if any of these vertices is in a connected component of original tree, we are not done
    std::vector<AmrVertexId> vertices_to_check;

    // merge all received trees
    for (size_t i = 0; i < received_trees.size(); ++i)
    {
        if (received_edges[i].empty())
            continue;
        AmrTripletMergeTree& rt = received_trees[i];

#ifdef DO_DETAILED_TIMING
        merge_timer.restart();
#endif

        r::merge(b->current_merge_tree_, rt, received_edges[i], true);

#ifdef DO_DETAILED_TIMING
        b->merge_call_time += merge_timer.elapsed();
#endif

        if (debug) fmt::print( "In receive_simple for block = {}, merge and repair OK for sender = {}, tree size = {}\n", b->gid, sender_gids_debug[i], b->get_merge_tree().size());

        // save information about vertex-component relation and component merging in block
        b->original_vertex_to_deepest_.insert(received_vertex_to_deepest[i].begin(),
                                              received_vertex_to_deepest[i].end());

#ifdef DO_DETAILED_TIMING
        union_find_timer.restart();
#endif

        // make_set in disjoint-sets of components
//        for (const AmrVertexId& new_deepest_vertex : received_deepest_vertices[i])
//        {
//            b->add_component_to_disjoint_sets(new_deepest_vertex);
//        }

#ifdef DO_DETAILED_TIMING
        b->union_find_time += union_find_timer.elapsed();
#endif

        //if (debug) fmt::print("In receive_simple for block = {}, add_component_to_disjoint_sets OK\n", b->gid);
        // update receivers - if a block from the 1-neighbourhood of sender has not received a tree from us
        // we must send it to this block in the next round

#ifdef DO_DETAILED_TIMING
        rl_loop_timer.restart();
#endif

        for (const AMRLink& rl : received_links)
        {
            for (int k = 0; k < rl.size(); ++k)
            {
                int sender_neighbor_gid = rl.target(k).gid;

                auto& original_gids = received_original_gids[i];

                bool is_in_original_gids = std::find(original_gids.begin(), original_gids.end(), sender_neighbor_gid) !=
                                           original_gids.end();
                bool is_not_processed = b->processed_receivers_.count(sender_neighbor_gid) == 0;

                //if (debug) fmt::print( "In receive_simple for block = {}, round = {}, sender_neighbor_gid = {}, is_in_original_gids = {}, is_not_processed = {}\n", b->gid, b->round_, sender_neighbor_gid, is_in_original_gids, is_not_processed);

                if (is_not_processed and is_in_original_gids)
                    b->new_receivers_.insert(sender_neighbor_gid);
            }
        }

#ifdef DO_DETAILED_TIMING
        b->rl_loop_time += rl_loop_timer.elapsed();
#endif

    }

#ifdef DO_DETAILED_TIMING
    b->whole_merge_tree_time += timer.elapsed();
    timer.restart();
#endif

    r::repair(b->current_merge_tree_);

#ifdef DO_DETAILED_TIMING
    b->repair_time += timer.elapsed();
    timer.restart();
#endif

//    if (debug) fmt::print("In receive_simple for block = {}, processed_receiveres_ OK\n", b->gid);

    // update disjoint sets data structure (some components are now connected to each other)
    for (size_t i = 0; i < received_trees.size(); ++i)
    {
        for (const AmrEdge& e : received_edges[i])
        {
            if (b->edge_exists(e))
            {
                // edge e connects two vertices that we have, connect their components
//                AmrVertexId deepest_a = b->original_deepest(std::get<0>(e));
//                AmrVertexId deepest_b = b->original_deepest(std::get<1>(e));
//                b->connect_components(deepest_a, deepest_b);
            } else
            {
                vertices_to_check.push_back(std::get<0>(e));
            }
        }
    }

#ifdef DO_DETAILED_TIMING
    b->union_find_time += timer.elapsed();
    timer.restart();
#endif

    b->sparsify_local_tree();

#ifdef DO_DETAILED_TIMING
    b->sparsify_time += timer.elapsed();
    timer.restart();
#endif

    if (debug) fmt::print("In receive_simple for block = {}, disjoint sets updated OK, tree size = {}\n", b->gid, b->get_merge_tree().size());

    expand_link(b, cp, l, received_links, received_original_gids);

#ifdef DO_DETAILED_TIMING
    b->expand_link_time += timer.elapsed();
    timer.restart();
#endif

    if (debug) fmt::print("Exit receive_simple for block = {}, expand_link OK\n", b->gid);

    b->done_ = b->is_done_simple(vertices_to_check);




    int n_undone = 1 - b->done_;

    cp.collectives()->clear();
    cp.all_reduce(n_undone, std::plus<int>());

#ifdef DO_DETAILED_TIMING
    b->is_done_time += timer.elapsed();
    timer.restart();
#endif

    if (debug)
        fmt::print("In receive_simple for block = {}, is_done_simple OK, vertices_to_check.size = {}\n", b->gid,
                   vertices_to_check.size());

    int old_size_unique = l->size_unique();
    int old_size = l->size();

    if (debug) fmt::print( "In receive_simple for block = {}, b->done_ = {}, old link size = {}, old link size_unqie = {}\n", b->gid, b->done_, old_size, old_size_unique);

}
