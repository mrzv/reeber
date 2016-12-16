#ifndef REEBER_MERGE_TREE_BLOCK_H
#define REEBER_MERGE_TREE_BLOCK_H

#include <diy/serialization.hpp>

#include <reeber/merge-tree.h>
#include <reeber/grid.h>
#include <reeber/grid-serialization.h>
#include <reeber/box.h>
#include <reeber/merge-tree-serialization.h>
#include <reeber/range/map.h>
namespace r = reeber;

#include "reeber-real.h"

struct MergeTreeBlock
{
    typedef     r::Grid<Real, 3>                Grid;
    typedef     r::OffsetGrid<Real, 3>          OffsetGrid;
    typedef     r::GridRestriction<Real, 3>     GridRestriction;
    typedef     Grid::Index                     Index;
    typedef     Grid::Vertex                    Vertex;
    typedef     Grid::Value                     Value;
    typedef     r::Box<3>                       Box;
    typedef     r::MergeTree<Index, Value>      MergeTree;
    typedef     std::vector<Real>               Size;

    static void*            create()                                        { return new MergeTreeBlock; }
    static void             destroy(void* b)                                { delete static_cast<MergeTreeBlock*>(b); }
    static void             save(const void* b, diy::BinaryBuffer& bb)      { diy::save(bb, *static_cast<const MergeTreeBlock*>(b)); }
    static void             load(      void* b, diy::BinaryBuffer& bb)      { diy::load(bb, *static_cast<MergeTreeBlock*>(b)); }

    inline void             compute_average(const diy::Master::ProxyWithLink& cp);

    int                     gid;
    Box                     core;
    Box                     local;
    Box                     global;
    MergeTree               mt;
    OffsetGrid              grid;
    Size                    cell_size;
};

namespace diy
{
    template<>
    struct Serialization<MergeTreeBlock>
    {
        static void             save(diy::BinaryBuffer& bb, const MergeTreeBlock& b)
        {
            diy::save(bb, b.gid);
            diy::save(bb, b.core);
            diy::save(bb, b.local);
            diy::save(bb, b.global);
            diy::save(bb, b.mt);
            diy::save(bb, b.grid);
            diy::save(bb, b.cell_size);
        }
        static void             load(diy::BinaryBuffer& bb, MergeTreeBlock& b)
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
MergeTreeBlock::
compute_average(const diy::Master::ProxyWithLink& cp)
{
    double                value = 0;
    size_t                count = 0;

    for(const MergeTree::Neighbor node : const_cast<const MergeTree&>(mt).nodes() | r::range::map_values)
    {
        if (core.contains(node->vertex))
        {
            value += node->value;
            count += 1;
        }

        for(const MergeTree::Node::ValueVertex& x : node->vertices)
            if (core.contains(x.second))
            {
                value += x.first;
                count += 1;
            }
    }

    cp.all_reduce(value, std::plus<double>());
    cp.all_reduce(count, std::plus<size_t>());
}

#endif
