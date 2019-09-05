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
    auto* l = static_cast<diy::AMRLink*>(cp.link());

    for(const diy::BlockID& receiver : link_unique(l, b->gid))
    {
        int receiver_gid = receiver.gid;
        cp.enqueue(receiver, b->gid_to_outgoing_edges_[receiver_gid]);
        cp.enqueue(receiver, b->vertex_to_deepest_);
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
}

template<class Real, unsigned D>
void expand_link(FabComponentBlock<Real, D>* b,
        const diy::Master::ProxyWithLink& cp,
        AMRLink* l,
        std::vector<AMRLink>& received_links,
        std::unordered_set<int> needed_gids)
{
    int n_added = 0;
    std::set<int> added_gids;

    for(size_t i = 0; i < received_links.size(); ++i)
    {
        const AMRLink& received_link = received_links[i];

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

            n_added++;
            l->add_neighbor(received_link.target(j));
            added_gids.insert(candidate_gid);
            l->add_bounds(received_link.level(j), received_link.refinement(j), received_link.core(j),
                    received_link.bounds(j));
        }
    }

    cp.master()->add_expected(n_added);
}

template<class Real, unsigned D>
void amr_cc_send(FabComponentBlock<Real, D>* b, const diy::Master::ProxyWithLink& cp)
{
    using Component = typename FabComponentBlock<Real, D>::Component;

    b->round_++;

    auto* l = static_cast<AMRLink*>(cp.link());
    auto receivers = link_unique(l, b->gid);

    for(const diy::BlockID& receiver : receivers)
    {
        diy::MemoryBuffer& out = cp.outgoing(receiver);
        diy::LinkFactory::save(out, l);

        int n_sent = 0;
        int receiver_gid = receiver.gid;

        for(Component& c : b->components_)
        {
            if (c.is_done_sending())
                continue;

            if (not c.must_send_to_gid(receiver_gid))
                continue;

            n_sent++;

            cp.enqueue(receiver, c.original_deepest());
            cp.enqueue(receiver, c.current_neighbors());
            int n_trees = c.must_send_tree_to_gid(receiver_gid);
            cp.enqueue(receiver, n_trees);
            if (n_trees)
            {
                cp.enqueue(receiver, c.tree_);
                cp.enqueue(receiver, c.edges());
#ifdef REEBER_EXTRA_INTEGRAL
                cp.enqueue(receiver, c.extra_values());
#endif
            }
        } // loop over components
    } // loop over receivers

    for(Component& c : b->components_)
    {
        // check that all neighbors are in link
//        if (not (std::all_of(c.processed_gids().begin(), c.processed_gids().end(),
//                [&receivers, b](const auto gid) {
//                    return b->gid == gid
//                            or std::any_of(receivers.begin(), receivers.end(), [gid](const diy::BlockID& receiver) {
//                                return receiver.gid == gid;
//                            });
//                }))) throw std::runtime_error("not all neighbors are in link!");

        c.mark_all_gids_processed();
    }
}

template<class Real, unsigned D>
void amr_cc_receive(FabComponentBlock<Real, D>* b, const diy::Master::ProxyWithLink& cp)
{

#ifdef REEBER_DO_DETAILED_TIMING
    dlog::Timer timer;
    dlog::Timer global_timer;
    dlog::Timer merge_call_timer;
    using DurationType = decltype(timer.elapsed());
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

    auto* l = static_cast<AMRLink*>(cp.link());

    LinkVector received_links;
    std::vector<TripletMergeTree> received_trees;

    auto senders = link_unique(l, b->gid);

    int total_received_trees = 0;

    AmrVertexContainer received_deepest_vertices;
    std::unordered_map<AmrVertexId, AmrVertexSet> received_root_to_components;

    AmrVertexSet keep;

#ifdef REEBER_DO_DETAILED_TIMING
    timer.restart();
#endif

    for(const diy::BlockID& sender : senders)
    {
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

#ifdef REEBER_EXTRA_INTEGRAL
            ExtraValues received_extra_values;
#endif
            cp.dequeue(sender, received_original_deepest);
            cp.dequeue(sender, received_current_neighbors);
            cp.dequeue(sender, received_n_trees);
            total_received_trees += received_n_trees;
            if (received_n_trees)
            {
                cp.dequeue(sender, received_tree);
                cp.dequeue(sender, received_edges);
#ifdef REEBER_EXTRA_INTEGRAL
                cp.dequeue(sender, received_extra_values);
#endif

#ifdef REEBER_DO_DETAILED_TIMING
                merge_call_timer.restart();
#endif

                r::merge(b->merge_tree_, received_tree, received_edges, true);

                for(auto e : received_edges)
                {
                    keep.insert(std::get<0>(e));
                    keep.insert(std::get<1>(e));
                }

#ifdef REEBER_DO_DETAILED_TIMING
                b->merge_call_time += merge_call_timer.elapsed();
                b->merge_calls += 1;
                b->edges_in_merge += received_edges.size();
#endif
            }

            received_deepest_vertices.push_back(received_original_deepest);
            received_root_to_components[received_original_deepest] = received_current_neighbors;

#ifdef REEBER_EXTRA_INTEGRAL
            if (!b->local_integral_.count(received_original_deepest))
                b->local_integral_[received_original_deepest] = received_extra_values;
            else
                assert(b->local_integral_[received_original_deepest] == received_extra_values);
#endif
        } // loop over all components sent from sender
    } // loop over senders

#ifdef REEBER_DO_DETAILED_TIMING
    b->process_senders_time += timer.elapsed();
    timer.restart();
#endif

    if (total_received_trees)
    {
#ifdef REEBER_DO_DETAILED_TIMING
        timer.restart();
#endif
        r::repair(b->merge_tree_);

#ifdef REEBER_DO_DETAILED_TIMING
        b->repair_time += timer.elapsed();
        timer.restart();
#endif

        b->update_connectivity(received_deepest_vertices);
#ifdef REEBER_DO_DETAILED_TIMING
        b->uc_time += timer.elapsed();
        timer.restart();
#endif
        b->sparsify_local_tree(keep);
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
            global_deepest_to_neighbors[global_deepest].insert(c.current_neighbors().begin(),
                    c.current_neighbors().end());
            global_deepest_to_original_deepests[global_deepest].insert(original_deepest);
        }

#ifdef REEBER_DO_DETAILED_TIMING
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

#ifdef REEBER_DO_DETAILED_TIMING
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

                for(const auto& deepest : cn)
                {
                    needed_gids.insert(deepest.gid);
                }
            }
        }

#ifdef REEBER_DO_DETAILED_TIMING
        b->ucn_loop_time += timer.elapsed();
        timer.restart();
#endif
    }

    expand_link(b, cp, l, received_links, needed_gids);

#ifdef REEBER_DO_DETAILED_TIMING
    b->expand_link_time += timer.elapsed();
    timer.restart();
#endif

    b->done_ = b->are_all_components_done();

    int undone = 1 - b->done_;

#ifdef REEBER_DO_DETAILED_TIMING
    b->is_done_time += timer.elapsed();
    timer.restart();
#endif

    cp.collectives()->clear();
    cp.all_reduce(undone, std::plus<int>());

#ifdef REEBER_DO_DETAILED_TIMING
    b->collectives_time += timer.elapsed();
    timer.restart();

    b->global_receive_time += global_timer.elapsed();
#endif
}
