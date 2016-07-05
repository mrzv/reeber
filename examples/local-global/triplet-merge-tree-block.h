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

    static void*            create()                                        { return new TripletMergeTreeBlock; }
    static void             destroy(void* b)                                { delete static_cast<TripletMergeTreeBlock*>(b); }
    static void             save(const void* b, diy::BinaryBuffer& bb)      { diy::save(bb, *static_cast<const TripletMergeTreeBlock*>(b)); }
    static void             load(      void* b, diy::BinaryBuffer& bb)      { diy::load(bb, *static_cast<TripletMergeTreeBlock*>(b)); }

    inline void             compute_average(const diy::Master::ProxyWithLink& cp, void*);

    int                     gid;
    Box                     core;
    Box                     local;
    Box                     global;
    TripletMergeTree        mt;
    OffsetGrid              grid;
    Size                    cell_size;
};

namespace diy
{
    template<>
    struct Serialization<TripletMergeTreeBlock>
    {
        static void             save(diy::BinaryBuffer& bb, const TripletMergeTreeBlock& b)
        {
            diy::save(bb, b.gid);
            diy::save(bb, b.core);
            diy::save(bb, b.local);
            diy::save(bb, b.global);
            diy::save(bb, b.mt);
            diy::save(bb, b.grid);
            diy::save(bb, b.cell_size);
        }
        static void             load(diy::BinaryBuffer& bb, TripletMergeTreeBlock& b)
        {
            diy::load(bb, b.gid);
            diy::load(bb, b.core);
            diy::load(bb, b.local);
            diy::load(bb, b.global);
            diy::load(bb, b.mt);
            diy::load(bb, b.grid);
            diy::load(bb, b.cell_size);
        }
    };
}

void
TripletMergeTreeBlock::
compute_average(const diy::Master::ProxyWithLink& cp, void*)
{
    double                value = 0;
    size_t                count = 0;

    BOOST_FOREACH(const TripletMergeTree::Neighbor node, const_cast<const TripletMergeTree&>(mt).nodes() | ba::map_values)
    {
        if (core.contains(node->vertex))
        {
            value += node->value;
            count += 1;
        }
    }

    cp.all_reduce(value, std::plus<double>());
    cp.all_reduce(count, std::plus<size_t>());
}

#endif
