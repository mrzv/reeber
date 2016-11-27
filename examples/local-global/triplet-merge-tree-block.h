#ifndef REEBER_TRIPLET_MERGE_TREE_BLOCK_H
#define REEBER_TRIPLET_MERGE_TREE_BLOCK_H

#include <diy/serialization.hpp>

#include <reeber/triplet-merge-tree.h>
#include <reeber/grid.h>
#include <reeber/grid-serialization.h>
#include <reeber/box.h>
#include <reeber/triplet-merge-tree-serialization.h>
#include <reeber/edges.h>
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

    using EdgeMap  = reeber::EdgeMap<Index, Value>;
    using EdgeMaps = reeber::EdgeMaps<Index, Value>;

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
    EdgeMaps                edge_maps;
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
