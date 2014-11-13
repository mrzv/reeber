#ifndef REEBER_MERGE_TREE_BLOCK_H
#define REEBER_MERGE_TREE_BLOCK_H

#include <diy/serialization.hpp>

#include <reeber/merge-tree.h>
#include <reeber/grid.h>
#include <reeber/box.h>
#include <reeber/merge-tree-serialization.h>
namespace r = reeber;


typedef     REEBER_REAL                     Real;

struct MergeTreeBlock
{
    typedef     r::Grid<Real, 3>                Grid;
    typedef     r::OffsetGrid<Real, 3>          OffsetGrid;
    typedef     Grid::Index                     Index;
    typedef     Grid::Vertex                    Vertex;
    typedef     Grid::Value                     Value;
    typedef     r::Box<3>                       Box;
    typedef     r::MergeTree<Index, Value>      MergeTree;

    static void*            create()                                        { return new MergeTreeBlock; }
    static void             destroy(void* b)                                { delete static_cast<MergeTreeBlock*>(b); }
    static void             save(const void* b, diy::BinaryBuffer& bb)      { diy::save(bb, *static_cast<const MergeTreeBlock*>(b)); }
    static void             load(      void* b, diy::BinaryBuffer& bb)      { diy::load(bb, *static_cast<MergeTreeBlock*>(b)); }

    int                     gid;
    Box                     local;
    Box                     global;
    MergeTree               mt;
};

namespace diy
{
    template<>
    struct Serialization<MergeTreeBlock>
    {
        static void             save(diy::BinaryBuffer& bb, const MergeTreeBlock& b)
        {
            diy::save(bb, b.gid);
            diy::save(bb, b.local);
            diy::save(bb, b.global);
            diy::save(bb, b.mt);
        }
        static void             load(diy::BinaryBuffer& bb, MergeTreeBlock& b)
        {
            diy::load(bb, b.gid);
            diy::load(bb, b.local);
            diy::load(bb, b.global);
            diy::load(bb, b.mt);
        }
    };
}

#endif
