#pragma once

#include <vector>

#include <diy/link.hpp>
#include <reeber/box.h>
#include "fab-cc-block.h"
#include "amr-merge-tree-helper.h"
using AMRLink = diy::AMRLink;

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
    b->sparsify_prune_original_tree();
}

template<class Real, unsigned D>
void expand_link(FabComponentBlock<Real, D>* b,
                 const diy::Master::ProxyWithLink& cp,
                 AMRLink* l,
                 std::vector<AMRLink>& received_links)
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
            // if we are already sending to this block, skip it
            int candidate_gid = received_link.target(j).gid;
            if (link_contains_gid(l, candidate_gid))
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

    //if (debug) fmt::print( "In expand_link for block = {}, round = {}, b->done_ = {}, n_added = {}, new link size = {}, new link size_unqie = {}, added_gids = {}\n", b->gid, b->round_, b->done_, n_added, l->size(), l->size_unique(), container_to_string(added_gids));
    if (debug)
        fmt::print( "In expand_link for block = {}, round = {}, b->done_ = {}, n_added = {}, new link size = {}, new link size_unqie = {}\n", b->gid, b->round_, b->done_, n_added, l->size(), l->size_unique());
    cp.master()->add_expected(n_added);
}

template<class Real, unsigned D>
void amr_cc_send(FabComponentBlock<Real, D>* b, const diy::Master::ProxyWithLink& cp)
{
    using Component = typename FabComponentBlock<Real, D>::Component;
    bool debug = false;
    if (debug) fmt::print("Called send_simple for block = {}\n", b->gid);

    auto* l = static_cast<AMRLink*>(cp.link());
    auto receivers = link_unique(l, b->gid);

    if (debug) fmt::print("In send_simple for block = {}, link size = {}, unique = {}\n", b->gid, l->size(), receivers.size());

    for (const diy::BlockID& receiver : receivers)
    {
        int receiver_gid = receiver.gid;

        diy::MemoryBuffer& out = cp.outgoing(receiver);
        diy::LinkFactory::save(out, l);

        int n_components = b->get_n_components_for_gid(receiver_gid);
        //if (debug) fmt::print("In send_simple for block = {}, sending to {}, n_components = {}\n", b->gid, receiver_gid, n_components);
        cp.enqueue(receiver, n_components);
        int n_sent = 0;
        for(Component& c : b->components_)
        {
            if (!c.must_send_to_gid(receiver_gid))
                continue;

            n_sent++;
            cp.enqueue(receiver, c.original_deepest());
            cp.enqueue(receiver, c.global_deepest());
            cp.enqueue(receiver, c.current_neighbors());
            cp.enqueue(receiver, c.processed_neighbors());
            cp.enqueue(receiver, b->original_integral_values_.at(c.original_deepest()));
            cp.enqueue(receiver, c.global_deepest_value());

            c.mark_gid_processed(receiver_gid);
        }
        assert(n_sent == n_components);
    }

    int done = b->done_;
    b->round_++;
    if (debug) fmt::print("Exit send_simple for block = {}, b->done = {}, b->round = {}\n", b->gid, done, b->round_);
}

template<class Real, unsigned D>
void amr_cc_receive(FabComponentBlock<Real, D>* b, const diy::Master::ProxyWithLink& cp)
{
    using Block = FabComponentBlock<Real, D>;
    using Component = typename FabComponentBlock<Real, D>::Component;
    using AmrEdgeVector = typename Block::AmrEdgeContainer;
    using LinkVector = std::vector<AMRLink>;
    using VertexValue = typename Block::VertexValue;
    using AmrVertexSet = typename Block::AmrVertexSet;

    bool debug = false;
    if (debug) fmt::print("Called amr_cc_receive for block = {}\n", b->gid);

    auto* l = static_cast<AMRLink*>(cp.link());

    cp.collectives()->clear();

    LinkVector received_links;

    auto senders = link_unique(l, b->gid);

    if (debug) fmt::print("In receive_simple for block = {}, # senders = {}\n", b->gid, senders.size());

    for (const diy::BlockID& sender : senders)
    {
        int n_components;

        diy::MemoryBuffer& in = cp.incoming(sender.gid);
        AMRLink* l = static_cast<AMRLink*>(diy::LinkFactory::load(in));
        received_links.push_back(*l);
        delete l;

        cp.dequeue(sender, n_components);

        for(int i = 0; i < n_components; ++i)
        {
            AmrVertexSet received_current_neighbors;
            AmrVertexSet received_processed_neighbors;
            AmrVertexId received_original_deepest;
            AmrVertexId received_global_deepest;
            Real received_original_integral_value;
            Real received_global_deepest_value;

            cp.dequeue(sender, received_original_deepest);
            cp.dequeue(sender, received_global_deepest);
            cp.dequeue(sender, received_current_neighbors);
            cp.dequeue(sender, received_processed_neighbors);
            cp.dequeue(sender, received_original_integral_value);
            cp.dequeue(sender, received_global_deepest_value);

            b->disjoint_sets_.make_component_if_not_exists(received_original_deepest);


//            {
//                auto x = b->original_integral_values_.count(received_original_deepest);
//                if (x != 0)
//                {
//                    auto v1 = b->original_integral_values_.at(received_original_deepest);
//                    auto diff = v1 - received_original_integral_value;
//                    assert(diff == 0.0);
//                }
//            }

            assert(b->original_integral_values_.count(received_original_deepest) == 0 or
                   b->original_integral_values_.at(received_original_deepest) == received_original_integral_value);
            b->original_integral_values_[received_original_deepest] = received_original_integral_value;

            for(const AmrVertexId& rcn : received_current_neighbors)
            {
                b->disjoint_sets_.make_component_if_not_exists(rcn);
                b->disjoint_sets_.unite_components(rcn, received_original_deepest);

                if (rcn.gid != b->gid)
                    continue;

                Component& c = b->get_component_by_deepest(rcn);
                assert(rcn == c.original_deepest_);
                if (c.processed_neighbors().count(received_original_deepest) == 0)
                {
                    for(const auto& cn : received_current_neighbors)
                        c.add_current_neighbor(cn);
                    c.mark_neighbor_processed(received_original_deepest);
                    b->disjoint_sets_.unite_components(c.original_deepest_, received_original_deepest);
                    if (b->cmp(received_global_deepest_value, c.global_deepest_value_)) {
                        c.set_global_deepest({received_global_deepest, received_global_deepest_value});
                    }
                }
            }
        }
    }


    // for debug only
    for(const Component& c : b->components_)
    {
        for(const AmrVertexId cn : c.current_neighbors())
        {
            assert(b->disjoint_sets_.are_connected(c.original_deepest_, cn));
        }
    }

    {
        std::unordered_map<AmrVertexId, AmrVertexSet> root_to_neighbors;  // to accumulate current_neighbors of all local components that are in one global component
        std::unordered_map<AmrVertexId, AmrVertexSet> root_to_deepest;    // to accumulate deepest vertices of all local components that are in one global component
        std::unordered_map<AmrVertexId, VertexValue> root_to_global_deepest;  // to update global_deepest in all local components that are in one global component
        for(const Component& c : b->components_)
        {
            auto root = b->disjoint_sets_.find_component(c.original_deepest_);
            root_to_neighbors[root].insert(c.current_neighbors().begin(), c.current_neighbors_.end());
            root_to_deepest[root].insert(c.global_deepest_);
            auto iter = root_to_global_deepest.find(root);
            if (iter != root_to_global_deepest.end())
            {
                if (b->cmp(c.global_deepest_value_, iter->second.value))
                    iter->second = { c.global_deepest_, c.global_deepest_value_ };
            }
            else
            {
                root_to_global_deepest[root] = { c.global_deepest_, c.global_deepest_value_ };
            };
        }

        for(auto root_deepest_set_pair : root_to_deepest)
        {
            AmrVertexId root = root_deepest_set_pair.first;
            for(auto deepest : root_deepest_set_pair.second)
            {
                if (deepest.gid != b->gid)
                    continue;

                Component& c = b->get_component_by_deepest(deepest);

                if (b->gid == 1)
                {
                    fmt::print("Before setting global deepest , c = {}\n", c);
                }

                c.set_current_neighbors(root_to_neighbors.at(root));
                c.set_global_deepest(root_to_global_deepest.at(root));

                if (b->gid == 1)
                {
                    fmt::print("After setting global deepest , c = {}\n", c);
                }
            }
        }
    }

    // exchange information between nodes in one block
    // only update integral values
    for(Component& c : b->components_)
    {
        for(const AmrVertexId& local_deepest : c.current_neighbors())
        {
            if (local_deepest.gid != b->gid || c.processed_neighbors().count(local_deepest))
                continue;
            Component& other_component = b->get_component_by_deepest(local_deepest);
            b->disjoint_sets_.unite_components(c.original_deepest_, other_component.original_deepest_);
            c.mark_neighbor_processed(local_deepest);
        }
    }

//    if (debug) fmt::print("In receive_simple for block = {}, processed_receiveres_ OK\n", b->gid);

    b->done_ = b->are_all_components_done();
    int n_undone = 1 - b->done_;

//    debug = b->gid == 0;
    if (debug) fmt::print("b->done = {}, n_undone = {}\n", b->done_, n_undone);

    cp.all_reduce(n_undone, std::plus<int>());

    int old_size_unique = l->size_unique();
    int old_size = l->size();

    //if (debug) fmt::print( "In receive_simple for block = {}, b->done_ = {}, old link size = {}, old link size_unqie = {}\n", b->gid, b->done_, old_size, old_size_unique);

    expand_link(b, cp, l, received_links);

    if (debug) fmt::print("Exit receive_simple for block = {}, expand_link OK\n", b->gid);

    std::vector<AmrVertexId> vv_check { {7, 34947}, {1, 32005}, {2, 12113}, {7, 34947}, {1, 32005},
                                        {2, 12113}, {3, 4773}, {1, 16980}, {1, 32005}, {3, 4773},
                                        {1, 16980}, {0, 33913}, {1, 34350} };

}
