#pragma once

#include <tuple>
#include <unordered_map>
#include <unordered_set>

#include <diy/master.hpp>
#include <diy/reduce.hpp>
#include <diy/partners/swap.hpp>

#include "triplet-merge-tree.h"
#include "edges.h"

namespace reeber
{

namespace detail
{
    template<class Block, class Vertex, class Value, class Partners, class GidGenerator>
    struct MergeSparsify;

    template<class Topology, class LocalTest>
    struct LocalTopology;

    template<class Topology, class LocalTest>
    LocalTopology<Topology, LocalTest>
    make_local_topology(const Topology& t, const LocalTest& lt)     { return LocalTopology<Topology, LocalTest>(t, lt); }
}

template<class Block, class Vertex, class Value>
void
resolve_edges(diy::Master&                            master,
              TripletMergeTree<Vertex,Value> Block::* tmt_,
              EdgeMaps<Vertex,Value> Block::*         edge_maps_);

template<class Block, class Vertex, class Value, class GidGenerator, class Partners>
void resolve_and_merge(diy::Master&                            master,
                       diy::Assigner&                          assigner,
                       TripletMergeTree<Vertex,Value> Block::* tmt_,
                       EdgeMaps<Vertex,Value> Block::*         edge_maps_,
                       const GidGenerator&                     gid_generator,
                       const Partners&                         partners)
{
    resolve_edges<Block, Vertex, Value>(master, tmt_, edge_maps_);

    // remove (simplify into vertices vectors) degree-2 nodes
    master.foreach([&](Block* b, const diy::Master::ProxyWithLink& cp)
    {
        auto& tmt       = b->*tmt_;
        auto& edge_maps = b->*edge_maps_;

        std::unordered_set<Vertex> special;
        for (auto& kv_em : edge_maps)
        {
            auto& edges = kv_em.second;
            for (auto &kv : edges)
            {
                Vertex s = std::get<1>(kv.second);
                if (tmt.contains(s)) special.insert(s);
            }
        }
        remove_degree_two(tmt, [&special](Vertex u) { return special.find(u) != special.end(); });
        LOG_SEV(debug) << "[" << b->gid << "] " << "Tree size after pruning degree-2: " << tmt.size();
    });

    // perform the global swap-reduce
    diy::reduce(master, assigner, partners,
                detail::MergeSparsify<Block, Vertex, Value, Partners, GidGenerator>(tmt_, edge_maps_, gid_generator));
}

template<class Block, class Vertex, class Value, class GidGenerator>
void resolve_and_merge(diy::Master&                            master,
                       diy::Assigner&                          assigner,
                       TripletMergeTree<Vertex,Value> Block::* tmt_,
                       EdgeMaps<Vertex,Value> Block::*         edge_maps_,
                       const GidGenerator&                     gid_generator)
{
    // By default, use 1-D domain decomposition. Clearly inefficient, but the
    // best we can hope for in absence of other assumptions.
    int nblocks = assigner.nblocks();
    diy::RegularDecomposer<diy::DiscreteBounds>  decomposer(1, diy::interval(0, nblocks - 1), nblocks);
    resolve_and_merge(master, assigner, tmt_, edge_maps_, gid_generator, diy::RegularSwapPartners(decomposer, 2, true));
}

template<class Block, class Vertex, class Value,
         class TopologyGenerator, class FunctionGenerator, class GidGenerator,
         class Partners>
void compute_merge_tree(diy::Master&                            master,
                        diy::Assigner&                          assigner,
                        TripletMergeTree<Vertex,Value> Block::* tmt_,
                        EdgeMaps<Vertex,Value> Block::*         edge_maps_,
                        const TopologyGenerator&                topology_generator,
                        const FunctionGenerator&                function_generator,
                        const GidGenerator&                     gid_generator,
                        const Partners&                         partners)
{
    // This makes two passes over the topology: first to compute the local
    // tree; second to figure out outgoing edges. Not super-efficient, but good
    // enough for now.

    // compute local tree
    master.foreach([&](Block* b, const diy::Master::ProxyWithLink& cp)
    {
        using Topology      = decltype(topology_generator(b));
        using Function      = decltype(function_generator(b));
        using GidGen        = decltype(gid_generator(b));

        Topology topology   = topology_generator(b);
        Function function   = function_generator(b);
        GidGen   gid_gen    = gid_generator(b);

        LOG_SEV(debug) << "Got topology: " << topology;

        int   gid = cp.gid();
        auto& tmt = b->*tmt_;
        compute_merge_tree2(tmt, detail::make_local_topology(topology, [&](Vertex v) { return gid_gen(v) == gid; }), function);
        LOG_SEV(debug) << "[" << b->gid << "] " << "Initial tree size: " << tmt.size();
    });

    // fill edges
    master.foreach([&](Block* b, const diy::Master::ProxyWithLink& cp)
    {
        auto    topology    = topology_generator(b);
        auto    gid_gen     = gid_generator(b);
        int     gid         = cp.gid();
        auto&   edge_maps   = b->*edge_maps_;

        for (auto u : topology.vertices())
        {
            if (gid_gen(u) != gid)
                continue;

            for (auto v : topology.link(u))
            {
                int v_gid = gid_gen(v);
                if (v_gid == gid)      // need non-local neighbors
                    continue;

                using ValueVertex   = std::tuple<Value, Vertex>;
                edge_maps[v_gid].emplace(std::make_tuple(u,v), ValueVertex());     // just keep the edges
            }
        }
    });

    resolve_and_merge(master, assigner, tmt_, edge_maps_, gid_generator, partners);
}

template<class Block, class Vertex, class Value, class TopologyGenerator, class FunctionGenerator, class GidGenerator>
void compute_merge_tree(diy::Master&                            master,
                        diy::Assigner&                          assigner,
                        TripletMergeTree<Vertex,Value> Block::* tmt_,
                        EdgeMaps<Vertex,Value> Block::*         edge_maps_,
                        const TopologyGenerator&                topology_generator,
                        const FunctionGenerator&                function_generator,
                        const GidGenerator&                     gid_generator)
{
    // By default, use 1-D domain decomposition. Clearly inefficient, but the
    // best we can hope for in absence of other assumptions.
    int nblocks = assigner.nblocks();
    diy::RegularDecomposer<diy::DiscreteBounds>  decomposer(1, diy::interval(0, nblocks - 1), nblocks);
    compute_merge_tree(master, assigner, tmt_, edge_maps_,
                       topology_generator, function_generator, gid_generator,
                       diy::RegularSwapPartners(decomposer, 2, true));
}

}

#include "distributed-tmt.hpp"
