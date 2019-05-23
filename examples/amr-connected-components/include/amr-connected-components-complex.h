#pragma once

#include <vector>

#include <diy/link.hpp>
#include <reeber/box.h>
#include <reeber/amr_helper.h>
#include "fab-cc-block.h"

using AMRLink = diy::AMRLink;
using AmrVertexId = reeber::AmrVertexId;
using AmrEdgeContainer = reeber::AmrEdgeContainer;

template<class R, unsigned D>
void send_edges_to_neighbors_cc(FabComponentBlock<R, D>* b, const diy::Master::ProxyWithLink& cp)
{
    bool debug = false;

    if(debug) fmt::print("Called send_edges_to_neighbors for block = {}\n", b->gid);

    auto* l = static_cast<diy::AMRLink*>(cp.link());

    for(const diy::BlockID& receiver : link_unique(l, b->gid))
    {
        int receiver_gid = receiver.gid;
        cp.enqueue(receiver, b->gid_to_outgoing_edges_[receiver_gid]);
        cp.enqueue(receiver, b->vertex_to_deepest_);
        if(debug) fmt::print("In send_edges_to_neighbors for block = {}, sending to receiver= {}, cardinality = {}\n", b->gid, receiver_gid, b->gid_to_outgoing_edges_[receiver_gid].size());
    }
}

template<class R, unsigned D>
void delete_low_edges_cc(FabComponentBlock<R, D>* b, const diy::Master::ProxyWithLink& cp)
{
    using VertexVertexMap = typename FabComponentBlock<R, D>::VertexVertexMap;
    auto* l = static_cast<diy::AMRLink*>(cp.link());

    for(const diy::BlockID& sender : link_unique(l, b->gid))
    {
        AmrEdgeContainer edges_from_neighbor;
        VertexVertexMap received_vertex_to_deepest;
        cp.dequeue(sender, edges_from_neighbor);
        cp.dequeue(sender, received_vertex_to_deepest);
        b->delete_low_edges(sender.gid, edges_from_neighbor, received_vertex_to_deepest);
    }

    b->adjust_outgoing_edges();
//    b->sparsify_prune_original_tree();
}

template<class Real, unsigned D>
void expand_link(FabComponentBlock<Real, D>* b,
                 const diy::Master::ProxyWithLink& cp,
                 AMRLink* l,
                 std::vector<AMRLink>& received_links,
                 std::unordered_set<int> needed_gids)
{
//    bool debug = (b->gid == 3) || (b->gid == 11) || (b->gid == 0) || (b->gid == 1);
    bool debug = false;
    if (debug) fmt::print("in expand_link for block = {}, round = {}, started updating link\n", b->gid, b->round_);
    int n_added = 0;
    std::set<int> added_gids;

    for(size_t i = 0; i < received_links.size(); ++i)
    {
        const AMRLink& received_link = received_links[i];
        //if (debug) fmt::print("in expand_link for block = {}, i = {}, received_links[i].size = {}\n", b->gid, i, received_links[i].size());

        for(int j = 0; j < received_link.size(); ++j)
        {
            int candidate_gid = received_link.target(j).gid;

            // skip yourself
            if (candidate_gid == b->gid)
                continue;

            // if we are already sending to this block, skip it
            if (link_contains_gid(l, candidate_gid))
                continue;

            if (needed_gids.count(candidate_gid) == 0)
                continue;

            //if (debug) fmt::print("in expand_link for block = {}, candidate_gid = {}, received_original_links[{}] = {}\n", b->gid, candidate_gid, i, container_to_string(received_original_gids[i]));
            n_added++;
            l->add_neighbor(received_link.target(j));
            added_gids.insert(candidate_gid);
            l->add_bounds(received_link.level(j), received_link.refinement(j), received_link.core(j),
                          received_link.bounds(j));
            //if (debug) fmt::print("in expand_link for block = {}, added gid = {}\n", b->gid, candidate_gid);
        }
    }

    //if (debug) fmt::print( "In expand_link for block = {}, round = {}, b->done_ = {}, n_added = {}, new link size = {}, new link size_unique = {}, added_gids = {}\n", b->gid, b->round_, b->done_, n_added, l->size(), l->size_unique(), container_to_string(added_gids));
    if (debug)
        fmt::print( "In expand_link for block = {}, round = {}, b->done_ = {}, n_added = {}, new link size = {}, new link size_unique = {}\n", b->gid, b->round_, b->done_, n_added, l->size(), l->size_unique());

    if (debug)
    {
        auto new_link = link_unique(l, b->gid);
        std::set<int> new_gids;
        for(diy::BlockID bid : new_link)
            new_gids.insert(bid.gid);
        fmt::print("links[{}] = {} #PYTHON round = {}\n", b->gid, container_to_string(new_gids), b->round_);
    }

    cp.master()->add_expected(n_added);
}

template<class Real, unsigned D>
void amr_cc_send(FabComponentBlock<Real, D>* b, const diy::Master::ProxyWithLink& cp)
{
    using Component = typename FabComponentBlock<Real, D>::Component;

    b->round_++;

    bool debug = false; //b->gid == 23938 or b->gid == 23945;

    auto* l = static_cast<AMRLink*>(cp.link());
    auto receivers = link_unique(l, b->gid);

    for (const diy::BlockID& receiver : receivers)
    {
        int receiver_gid = receiver.gid;

        diy::MemoryBuffer& out = cp.outgoing(receiver);
        diy::LinkFactory::save(out, l);

        int n_sent = 0;
        for(Component& c : b->components_)
        {
            if (debug) fmt::print("in amr_cc_send, c = {}, receiver gid = {} START\n", c.original_deepest(), receiver_gid);
            if (c.is_done_sending())
                continue;

            if (not c.must_send_to_gid(receiver_gid))
                continue;

//            debug = ( debug_v == c.original_deepest() or debug_v_1 == c.original_deepest());

//            if (!c.must_send_tree_to_gid(receiver_gid))
//            {
//                if (debug) fmt::print("in amr_cc_send, not sending {} to gid {}\n", c.original_deepest(), receiver_gid);
//                continue;
//            }

            n_sent++;

            if (debug) fmt::print("in amr_cc_send, c = {}, receiver gid = {}, #current_neighbors = {}, #edges = {}\n",
                    c.original_deepest(), receiver_gid, c.current_neighbors().size(), c.edges().size());

            cp.enqueue(receiver, c.original_deepest());
            cp.enqueue(receiver, c.current_neighbors());
            int n_trees = c.must_send_tree_to_gid(receiver_gid);
            if (debug) fmt::print("in amr_cc_send, c = {}, receiver gid = {}, sent original deepest, n_trees = {}\n", c.original_deepest(), receiver_gid, n_trees);
            if (debug) fmt::print("in amr_cc_send, cfull = {}\n", c);
            cp.enqueue(receiver, n_trees);
            if (n_trees)
            {
                cp.enqueue(receiver, c.tree_);
                cp.enqueue(receiver, c.edges());
#ifdef EXTRA_INTEGRAL
                cp.enqueue(receiver, c.extra_values());
#endif
            }
//            if (debug) fmt::print("in amr_cc_send, sent {} to gid = {}, current_neighbors = {}\n", c.original_deepest(), receiver_gid, container_to_string(c.current_neighbors()));

        }
//        assert(n_sent == n_components);
    }

    for(Component& c : b->components_)
    {
        // check that all neighbors are in link
        assert(std::all_of(c.processed_gids().begin(), c.processed_gids().end(),
                [&receivers, b](const auto gid) { return b->gid == gid or std::any_of(receivers.begin(), receivers.end(), [gid](const diy::BlockID& receiver) {
            return receiver.gid == gid; }); } ));

        c.mark_all_gids_processed();
    }

    if (debug) fmt::print("Exit send_components for block = {}, b->done = {}, b->round = {}\n", b->gid, b->done_, b->round_);
}

template<class Real, unsigned D>
void amr_cc_receive(FabComponentBlock<Real, D>* b, const diy::Master::ProxyWithLink& cp)
{
#ifdef DO_DETAILED_TIMING
    dlog::Timer timer;
    dlog::Timer global_timer;
    dlog::Timer merge_call_timer;
    // detailed timings
    using DurationType = decltype(timer.elapsed());

//    DurationType time_to_receive_trees_and_gids;
//    DurationType rl_loop_time;
//    DurationType repair_time;
//    DurationType merge_call_time;
//    DurationType union_find_time;
//    DurationType sparsify_time;
//    DurationType expand_link_time;

    timer.restart();
    global_timer.restart();

#endif

    using Block = FabComponentBlock<Real, D>;
    using Component = typename FabComponentBlock<Real, D>::Component;
    using LinkVector = std::vector<AMRLink>;
    using AmrVertexSet = typename Block::AmrVertexSet;
    using AmrVertexContainer = typename Block::AmrVertexContainer;
    using TripletMergeTree = typename Block::TripletMergeTree;
    using GidSet = typename Block::GidSet;
    using ExtraValues = typename Block::ExtraValues;

    bool debug = false; //b->gid == 23938 or b->gid == 23945;
    if (debug) fmt::print("Called amr_cc_receive for block = {}, round = {}\n", b->gid, b->round_);

    //if (debug)
    //    for(const Component& c : b->components_)
    //    {
    //        fmt::print("Component: {}, #current_neighbors = {}, #processed_neighbors = {}\n",
    //                c.original_deepest(), c.current_neighbors().size(), c.processed_neighbors().size());
    //    }

    auto* l = static_cast<AMRLink*>(cp.link());


    LinkVector received_links;
    std::vector<TripletMergeTree> received_trees;

    auto senders = link_unique(l, b->gid);

    int total_received_trees = 0;

//    if (debug) fmt::print("In receive_simple for block = {}, # senders = {}\n", b->gid, senders.size());

    AmrVertexContainer received_deepest_vertices;
    std::unordered_map<AmrVertexId, AmrVertexSet> received_root_to_components;

#ifdef DO_DETAILED_TIMING
        timer.restart();
#endif


    for (const diy::BlockID& sender : senders)
    {
        //if (debug) fmt::print("amr_cc_receive, round = {}, block {}, receiving from {}\n", b->round_, b->gid, sender.gid);
        diy::MemoryBuffer& in = cp.incoming(sender.gid);
        AMRLink* l = static_cast<AMRLink*>(diy::LinkFactory::load(in));
        received_links.push_back(*l);
        delete l;

        while(cp.incoming(sender.gid))
        {
            AmrVertexSet received_current_neighbors;
            AmrVertexId received_original_deepest;
            int received_n_trees;
            TripletMergeTree received_tree;
            AmrEdgeContainer received_edges;
            ExtraValues received_extra_values;

            cp.dequeue(sender, received_original_deepest);
            cp.dequeue(sender, received_current_neighbors);
            cp.dequeue(sender, received_n_trees);
            total_received_trees += received_n_trees;
            if (received_n_trees)
            {
                cp.dequeue(sender, received_tree);
                cp.dequeue(sender, received_edges);
                cp.dequeue(sender, received_extra_values);

#ifdef DO_DETAILED_TIMING
                merge_call_timer.restart();
#endif

                r::merge(b->merge_tree_, received_tree, received_edges, true);

#ifdef DO_DETAILED_TIMING
                b->merge_call_time += merge_call_timer.elapsed();
#endif
            }

            received_deepest_vertices.push_back(received_original_deepest);
            received_root_to_components[received_original_deepest] = received_current_neighbors;

#ifdef EXTRA_INTEGRAL
            if (!b->local_integral_.count(received_original_deepest))
                b->local_integral_[received_original_deepest] = received_extra_values;
            else
                assert(b->local_integral_[received_original_deepest] == received_extra_values);
#endif

        }
    }

#ifdef DO_DETAILED_TIMING
        b->process_senders_time += timer.elapsed();
        timer.restart();
#endif


    if (total_received_trees)
    {
#ifdef DO_DETAILED_TIMING
        timer.restart();
#endif
        r::repair(b->merge_tree_);

#ifdef DO_DETAILED_TIMING
        b->repair_time += timer.elapsed();
        timer.restart();
#endif

        b->update_connectivity(received_deepest_vertices);
#ifdef DO_DETAILED_TIMING
        b->uc_time += timer.elapsed();
        timer.restart();
#endif
    }

    GidSet needed_gids;

    // process internal components that are united after merging
    {
        std::unordered_map<AmrVertexId, AmrVertexSet> global_deepest_to_neighbors;  // to accumulate current_neighbors of all local components that are in one global component
        std::unordered_map<AmrVertexId, AmrVertexSet> global_deepest_to_original_deepests;

         // collect all neighbors in merged components and corresponding roots
        for(const Component& c : b->components_)
        {
            AmrVertexId original_deepest = c.original_deepest();
            AmrVertexId global_deepest = b->vertex_to_deepest_.at(original_deepest);
            global_deepest_to_neighbors[global_deepest].insert(c.current_neighbors().begin(), c.current_neighbors().end());
            global_deepest_to_original_deepests[global_deepest].insert(original_deepest);
       }

#ifdef DO_DETAILED_TIMING
        b->comps_loop_time += timer.elapsed();
        timer.restart();
#endif
        for(const auto& root_deepest_set_pair : received_root_to_components)
       {
           const AmrVertexId& original_deepest = root_deepest_set_pair.first;
           AmrVertexId global_deepest = b->vertex_to_deepest_.at(original_deepest);
           const AmrVertexSet& cn = root_deepest_set_pair.second;
           global_deepest_to_neighbors[global_deepest].insert(cn.begin(), cn.end());
       }

#ifdef DO_DETAILED_TIMING
        b->rrtc_loop_time += timer.elapsed();
        timer.restart();
#endif

       // update current neighbors
       for(const auto& deepest_roots_pair : global_deepest_to_original_deepests)
       {
           AmrVertexId global_deepest = deepest_roots_pair.first;
           for(AmrVertexId original_deepest : deepest_roots_pair.second)
           {
               Component& c = b->get_component_by_deepest(original_deepest);
               auto& cn = global_deepest_to_neighbors.at(global_deepest);
               c.set_current_neighbors(cn);

               for(const auto& deepest : cn){
                   needed_gids.insert(deepest.gid);
               }
           }
       }

#ifdef DO_DETAILED_TIMING
        b->ucn_loop_time += timer.elapsed();
        timer.restart();
#endif
    }

//    int old_size_unique = l->size_unique();
//    int old_size = l->size();

    expand_link(b, cp, l, received_links, needed_gids);

#ifdef DO_DETAILED_TIMING
        b->expand_link_time += timer.elapsed();
        timer.restart();
#endif

    b->done_ = b->are_all_components_done();
    int undone = 1 - b->done_;

#ifdef DO_DETAILED_TIMING
        b->is_done_time += timer.elapsed();
        timer.restart();
#endif

    if (debug)
        for(const Component& c : b->components_)
        {
            fmt::print("amr_cc_receive END: Component: {}, #current_gids = {}, #processed_gids = {}, must_send_neighbors = {}\n",
                    c.original_deepest(), c.current_gids().size(), c.processed_gids().size(), c.must_send_neighbors());
        }

    cp.collectives()->clear();
    cp.all_reduce(undone, std::plus<int>());

#ifdef DO_DETAILED_TIMING
        b->collectives_time += timer.elapsed();
        timer.restart();

        b->global_receive_time += global_timer.elapsed();
#endif


    if (debug) fmt::print("Exit amr_cc_receive, gid = {}, round = {}, mt.size = {}, done = {}\n", b->gid, b->round_, b->get_merge_tree().size(), b->done_);

}
