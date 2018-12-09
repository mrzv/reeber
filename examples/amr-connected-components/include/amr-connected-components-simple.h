#pragma once

#include <vector>

#include <diy/link.hpp>
#include <reeber/box.h>
#include "fab-cc-block.h"
#include "amr-merge-tree-helper.h"
using AMRLink = diy::AMRLink;

template<class Real, unsigned D>
void expand_link(FabComponentBlock<Real, D>* b,
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
        fmt::print( "In expand_link for block = {}, round = {}, b->done_ = {}, n_added = {}, new link size = {}, new link size_unqie = {}\n", b->gid, b->round_, b->done_, n_added, l->size(), l->size_unique());
    cp.master()->add_expected(n_added);
}

template<class Real, unsigned D>
void amr_tmt_send(FabComponentBlock<Real, D>* b, const diy::Master::ProxyWithLink& cp)
{

    //bool debug = (b->gid == 0);
    bool debug = false;
    if (debug) fmt::print("Called send_simple for block = {}\n", b->gid);

    auto* l = static_cast<AMRLink*>(cp.link());


    auto receivers = link_unique(l, b->gid);

    if (debug)
        fmt::print("In send_simple for block = {}, link size = {}, unique = {}\n", b->gid, l->size(), receivers.size());


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

        //if (debug) fmt::print("In send_simple for block = {}, sending to {}, n_trees = {}\n", b->gid, receiver_gid, n_trees);

        cp.enqueue(receiver, n_trees);

        if (n_trees)
        {
            // send local tree and all outgoing edges that end in receiver
            cp.enqueue(receiver, b->disjoint_sets_);
            cp.enqueue(receiver, b->get_all_outgoing_edges());
            cp.enqueue(receiver, b->vertex_values_);

            if (debug) fmt::print("In amr_tmt_send, enqued parent {}, size {}, edges {}", b->disjoint_sets_.parent_.size(), b->disjoint_sets_.size_.size(), b->get_all_outgoing_edges().size());

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
void amr_tmt_receive(FabComponentBlock<Real, D>* b, const diy::Master::ProxyWithLink& cp)
{
//    bool debug = (b->gid == 3) || (b->gid == 11) || (b->gid == 0) || (b->gid == 1);
    bool debug = false;

    //    if (debug) fmt::print("Called receive_simple for block = {}\n", b->gid);

    using Block = FabComponentBlock<Real, D>;
    using AmrEdgeVector = typename Block::AmrEdgeContainer;
    using VertexVertexMap = typename Block::VertexVertexMap;
    using VertexSizeMap = typename Block::VertexSizeMap;
    using LinkVector = std::vector<AMRLink>;
    using VertexValue = typename Block::VertexValue;
    using VertexDeepestMap = typename Block::VertexDeepestMap;
    using DjSets = typename DisjointSets<AmrVertexId>;

    auto* l = static_cast<AMRLink*>(cp.link());

    cp.collectives()->clear();

    std::vector<AmrEdgeVector> received_edges;
    std::vector<std::vector<int>> received_original_gids;

    LinkVector received_links;

    auto senders = link_unique(l, b->gid);

    if (debug) fmt::print("In receive_simple for block = {}, # senders = {}\n", b->gid, senders.size());

    VertexDeepestMap new_root_to_deepest { b->root_to_deepest_};


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

        //if (debug) fmt::print("In receive_simple for block = {}, dequeued from sender {} n_trees = {} \n", b->gid, sender.gid, n_trees);

        if (n_trees > 0)
        {
            assert(n_trees == 1);

            DjSets received_disjoint_sets;
            received_edges.emplace_back();
            VertexDeepestMap received_root_to_deepest;

            cp.dequeue(sender, received_disjoint_sets);
            cp.dequeue(sender, received_edges.back());
            cp.dequeue(sender, received_root_to_deepest);

            int old_size = b->disjoint_sets_.parent_.size();

            b->disjoint_sets_.disjoint_union(received_disjoint_sets);
            new_root_to_deepest.insert(received_root_to_deepest.begin(), received_root_to_deepest.end());

            if (debug) fmt::print("In amr_tmt_receive, old parent size {}, new parent size {}, edges received {}\n", old_size, b->disjoint_sets_.parent_.size(), received_edges.back().size());
//            if (debug) fmt::print( "In receive_simple for block = {}, dequeued from sender {} original link gids = {}\n", b->gid, sender.gid, container_to_string(received_original_gids.back()));

        }
    }

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
    for (size_t i = 0; i < received_edges.size(); ++i)
    {
        for(const AmrEdge& e : received_edges[i])
        {
            if (b->edge_exists(e))
            {
                auto x = std::get<0>(e);
                auto y = std::get<1>(e);
                auto x_root = b->disjoint_sets_.find_component(x);
                auto y_root = b->disjoint_sets_.find_component(y);
                if (x_root != y_root)
                {
                    auto survived_root = b->disjoint_sets_.unite_components_by_roots(x_root, y_root);
                    auto iter_x = new_root_to_deepest.find(x_root);
                    auto iter_y = new_root_to_deepest.find(y_root);
                    if (iter_x != new_root_to_deepest.end() and iter_y != new_root_to_deepest.end())
                    {
                        VertexValue vv = (b->cmp(iter_x->second.value, iter_y->second.value)) ? iter_x->second : iter_y->second;
                        if (survived_root == x_root)
                        {
                            iter_x->second = vv;
                            new_root_to_deepest.erase(iter_y);
                        } else
                        {
                            iter_y->second = vv;
                            new_root_to_deepest.erase(iter_x);
                        }
                    }
                }
            } else
            {
                vertices_to_check.push_back(std::get<0>(e));
            }
        }

        b->root_to_deepest_ = new_root_to_deepest;

        //if (debug) fmt::print("In receive_simple for block = {}, add_component_to_disjoint_sets OK\n", b->gid);
        // update receivers - if a block from the 1-neighbourhood of sender has not received a tree from us
        // we must send it to this block in the next round

        for (const AMRLink& rl : received_links)
        {
            for (int k = 0; k < rl.size(); ++k)
            {
                int sender_neighbor_gid = rl.target(k).gid;

                auto& original_gids = received_original_gids[i];

                bool is_in_original_gids = std::find(original_gids.begin(), original_gids.end(), sender_neighbor_gid) !=
                                           original_gids.end();
                bool is_not_processed = b->processed_receivers_.count(sender_neighbor_gid) == 0;

                if (debug) fmt::print( "In receive_simple for block = {}, round = {}, sender_neighbor_gid = {}, is_in_original_gids = {}, is_not_processed = {}\n", b->gid, b->round_, sender_neighbor_gid, is_in_original_gids, is_not_processed);

                if (is_not_processed and is_in_original_gids)
                    b->new_receivers_.insert(sender_neighbor_gid);
            }
        }
    }

    if (debug) fmt::print("In receive_simple for block = {}, processed_receiveres_ OK\n", b->gid);

    b->done_ = b->is_done_simple(vertices_to_check);
    int n_undone = 1 - b->done_;

    cp.all_reduce(n_undone, std::plus<int>());

    if (debug) fmt::print("In receive_simple for block = {}, is_done_simple OK, vertices_to_check.size = {}\n", b->gid, vertices_to_check.size());

    int old_size_unique = l->size_unique();
    int old_size = l->size();

    //if (debug) fmt::print( "In receive_simple for block = {}, b->done_ = {}, old link size = {}, old link size_unqie = {}\n", b->gid, b->done_, old_size, old_size_unique);

    expand_link(b, cp, l, received_links, received_original_gids);

    if (debug) fmt::print("Exit receive_simple for block = {}, expand_link OK\n", b->gid);
}
