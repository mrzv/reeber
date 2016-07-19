#ifndef REEBER_TRIPLET_MERGE_TREE_BLOCK_H
#define REEBER_TRIPLET_MERGE_TREE_BLOCK_H

#include <boost/range/adaptor/map.hpp>
namespace ba = boost::adaptors;

#include <diy/serialization.hpp>

#include <reeber/triplet-merge-tree.h>
#include <reeber/grid.h>
#include <reeber/grid-serialization.h>
#include <reeber/box.h>
#include <reeber/triplet-merge-tree-serialization.h>
namespace r = reeber;

#include "reeber-real.h"

struct TripletMergeTreeBlock
{
    typedef     r::Grid<Real, 3>                  Grid;
    typedef     r::OffsetGrid<Real, 3>            OffsetGrid;
    typedef     r::GridRestriction<Real, 3>       GridRestriction;
    typedef     Grid::Index                       Index;
    typedef     Grid::Vertex                      Vertex;
    typedef     Grid::Value                       Value;
    typedef     r::Box<3>                         Box;
    typedef     r::TripletMergeTree<Index, Value> TripletMergeTree;
    typedef     std::vector<Real>                 Size;

    typedef     std::tuple<Index, Index>          Edge;

    struct edge_hash : public std::unary_function<Edge, std::size_t>
    {
        std::size_t operator()(const Edge& k) const
        {
            size_t x = std::hash<Index>()(std::get<0>(k));
            size_t y = std::hash<Index>()(std::get<1>(k));
            return x ^ (y + 0x9e3779b9 + (x<<6) + (x>>2));
        }
    };

    struct edge_equal : public std::binary_function<Edge, Edge, bool>
    {
        bool operator()(const Edge& v0, const Edge& v1) const
        {
            return std::get<0>(v0) == std::get<0>(v1) && std::get<1>(v0) == std::get<1>(v1);
        }
    };
   
    typedef     std::unordered_map<Edge, std::tuple<Value, Index>, edge_hash, edge_equal>
                                                  EdgeMap;

    static void*            create()                                        { return new TripletMergeTreeBlock; }
    static void             destroy(void* b)                                { delete static_cast<TripletMergeTreeBlock*>(b); }
    static void             save(const void* b, diy::BinaryBuffer& bb)      { diy::save(bb, *static_cast<const TripletMergeTreeBlock*>(b)); }
    static void             load(      void* b, diy::BinaryBuffer& bb)      { diy::load(bb, *static_cast<TripletMergeTreeBlock*>(b)); }

    inline void             compute_average(const diy::Master::ProxyWithLink& cp, void*);

    int                     gid;
    Box                     local;
    Box                     global;
    TripletMergeTree        mt;
    OffsetGrid              grid;
    Size                    cell_size;
    EdgeMap                 edges;
    std::vector<EdgeMap>    edge_maps;
};

namespace diy
{
    template<>
    struct Serialization<TripletMergeTreeBlock>
    {
        static void             save(diy::BinaryBuffer& bb, const TripletMergeTreeBlock& b)
        {
            diy::save(bb, b.gid);
            diy::save(bb, b.local);
            diy::save(bb, b.global);
            diy::save(bb, b.mt);
            diy::save(bb, b.grid);
            diy::save(bb, b.cell_size);
            diy::save(bb, b.edges);
            diy::save(bb, b.edge_maps);
        }
        static void             load(diy::BinaryBuffer& bb, TripletMergeTreeBlock& b)
        {
            diy::load(bb, b.gid);
            diy::load(bb, b.local);
            diy::load(bb, b.global);
            diy::load(bb, b.mt);
            diy::load(bb, b.grid);
            diy::load(bb, b.cell_size);
            diy::load(bb, b.edges);
            diy::load(bb, b.edge_maps);
        }
    };
}

#endif
