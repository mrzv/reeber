#pragma once

#include <vector>

#include "amr-merge-tree-helper.h"
#include <diy/link.hpp>
#include "fab-tmt-block.h"
//
//#include "../../local-global/output-persistence.h"

template<class Real, unsigned D>
void expand_link(FabTmtBlock<Real, D>* b,
                 const diy::Master::ProxyWithLink& cp,
                 diy::AMRLink* l,
                 std::vector<diy::AMRLink>& received_links)
{
//    bool debug = (b->gid == 3) || (b->gid == 11) || (b->gid == 0) || (b->gid == 1);
    bool debug = false;
    if (debug) fmt::print("in expand_link for block = {}, round = {}, started updating link\n", b->gid, b->round_);
    int n_added = 0;
    assert(received_links.size() == received_original_gids.size());
    std::set<int> added_gids;

    for(size_t i = 0; i < received_links.size(); ++i)
    {
        const diy::AMRLink& received_link = received_links[i];
        assert(not received_original_gids[i].empty());
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
        fmt::print(
                "In expand_link for block = {}, round = {}, b->done_ = {}, n_added = {}, new link size = {}, new link size_unqie = {}\n",
                b->gid, b->round_, b->done_, n_added, l->size(), l->size_unique());
    cp.master()->add_expected(n_added);
}

template<class Real, unsigned D>
void amr_tmt_send(FabTmtBlock<Real, D>* b, const diy::Master::ProxyWithLink& cp)
{
    using Component = typename FabTmtBlock<Real, D>::Component;
    bool debug = true;
    if (debug) fmt::print("Called send_componentwise for block = {}\n", b->gid);
    auto* l = static_cast<diy::AMRLink*>(cp.link());
    cp.collectives()->clear();
    auto receivers = link_unique(l, b->gid);

    // find out how many components we send to each gid
    std::map<int, int> componenents_per_gid;
    for(const Component& c : b->components_)
    {
        for(const auto& receiver : receivers)
        {
            if (c.must_send_to_gid(receiver.gid))
                componenents_per_gid[receiver.gid]++;
        }
    }

    //    if (debug) fmt::print("In send_simple for block = {}, link size = {}, unique = {}\n", b->gid, l->size(), receivers.size());
    for(const diy::BlockID& receiver : receivers)
    {
        int n_trees = componenents_per_gid[receiver.gid];
        cp.enqueue(receiver, n_trees);
        if (n_trees == 0)
            continue;
        // send components that must communicate with the receiver
        for(Component& c : b->components_)
        {
            if (not c.must_send_to_gid(receiver.gid))
                continue;
            cp.enqueue(receiver, c.merge_tree_);
            cp.enqueue(receiver, c.root_);
            cp.enqueue(receiver, c.outgoing_edges_);
            cp.enqueue(receiver, c.current_neighbors_);
            c.processed_neighbors_.insert(receiver.gid);
        }
        cp.enqueue(receiver, b->original_vertex_to_deepest_);

        diy::MemoryBuffer& out = cp.outgoing(receiver);
        diy::LinkFactory::save(out, l);
    }
}

template<class Real, unsigned D>
void amr_tmt_receive(FabTmtBlock<Real, D>* b, const diy::Master::ProxyWithLink& cp)
{
    bool debug = true;
    //    if (debug) fmt::print("Called receive_componentwise for block = {}\n", b->gid);

    using Block = FabTmtBlock<Real, D>;
    using TripletMergeTree = typename Block::TripletMergeTree;
    using GidContainer = typename Block::GidContainer;
    using Component = typename Block::Component;
    using LinkVector = std::vector<diy::AMRLink>;


    auto* l = static_cast<diy::AMRLink*>(cp.link());

    AmrEdgeContainer all_received_edges;
    LinkVector received_links;

    auto senders = link_unique(l, b->gid);

    if (debug) fmt::print("In receive_componentwise for block = {}, # senders = {}\n", b->gid, senders.size());

    std::map<AmrVertexId, std::set<int>> component_to_neighbors;

    for(const diy::BlockID& sender : senders)
    {
        int n_trees;
        cp.dequeue(sender, n_trees);
        if (n_trees == 0)
            continue;
        for(int tree_idx = 0; tree_idx < n_trees; ++tree_idx)
        {
            TripletMergeTree received_tree;
            cp.dequeue(sender, received_tree);

            AmrVertexId received_root;
            cp.dequeue(sender, received_root);

            AmrEdgeContainer received_edges;
            cp.dequeue(sender, received_edges);
            all_received_edges.insert(all_received_edges.end(), received_edges.begin(), received_edges.end());

            GidContainer received_current_gids;
            cp.dequeue(sender, received_current_gids);

            // just for debug - check all received edges
            for(const AmrEdge& e : received_edges)
            {
                AmrVertexId my_vertex = std::get<1>(e);
                assert(my_vertex.gid != b->gid or
                       b->local_.mask_by_index(my_vertex) == MaskedBox::ACTIVE or
                       b->local_.mask_by_index(my_vertex) == MaskedBox::LOW);
            }

            // merge received trees
            r::merge(b->mt_, received_tree, received_edges, true);
            r::repair(b->mt_);
            if (debug) fmt::print("In receive_componentwise for block = {}, repair OK\n", b->gid);

            // make_set in disjoint-sets of components
            b->add_component_to_disjoint_sets(received_root);

            // add gids from sender
            component_to_neighbors[received_root].insert(received_current_gids.begin(),
                                                         received_current_gids.end());


        } // loop over all trees received from current sender

        typename Block::VertexVertexMap received_vertex_to_deepest;
        cp.dequeue(sender, received_vertex_to_deepest);
        // save information about vertex-component relation and component merging in block
        b->add_received_original_vertices(received_vertex_to_deepest);

        diy::MemoryBuffer& in = cp.incoming(sender.gid);
        diy::AMRLink* l = static_cast<diy::AMRLink*>(diy::LinkFactory::load(in));
        received_links.push_back(*l);
        delete l;
    }

    // all tree from all senders were processed, now we can connect components in the disjoint sets data structure
    for(const AmrEdge& e : all_received_edges)
    {
        if (b->edge_exists(e))
        {
            // edge e connects two vertices that we have, connect their components
            AmrVertexId deepest_a = b->original_deepest(std::get<0>(e));
            AmrVertexId deepest_b = b->original_deepest(std::get<1>(e));
            b->connect_components(deepest_a, deepest_b);
        } else
        {
            assert(b->edge_goes_out(e));
        }
    }

    all_received_edges.clear();

    // expand current neighbourhood of each component
    // first collect all incoming neighbourhoods
    for(Component& c : b->components_)
    {
        for(const auto& deepest : b->current_deepest_)
        {
            if (deepest.gid == b->gid)
                continue; // skip local roots for now
            if (b->are_components_connected(c.root_, deepest))
            {
                c.current_neighbors_.insert(component_to_neighbors[deepest].begin(),
                                            component_to_neighbors[deepest].end());
            }
        }
    }

    std::map<AmrVertexId, std::set<Component*>> m;
    std::map<AmrVertexId, std::set<int>> new_neighborhoods;
    std::set<AmrVertexId> s;
    // now neighbourhoods of merged components must be united
    for(Component& c : b->components_)
    {
        if (s.find(c.root_) != s.end())
            continue;

        if (m.find(c.root_) == m.end() and s.find(c.root_) == s.end())
        {
            m[c.root_].insert(&c);
            new_neighborhoods[c.root_] = c.current_neighbors_;
        }

        for(Component& other_c : b->components_)
        {
            if (other_c.root_ <= c.root_)
                continue; // skip yourself and already processed pairs
            s.insert(other_c.root_);
            if (b->are_components_connected(c.root_, other_c.root_))
            {
                m[c.root_].insert(&other_c);
                new_neighborhoods[c.root_].insert(other_c.current_neighbors_.begin(), other_c.current_neighbors_.end());
            }
        }
    }

    for(auto& root_component_set_pair : m)
    {
        for(Component* c : root_component_set_pair.second)
        {
            c->current_neighbors_ = new_neighborhoods.at(root_component_set_pair.first);
        }
    }
    expand_link(b, cp, l, received_links);
}
