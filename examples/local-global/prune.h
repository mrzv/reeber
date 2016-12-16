#ifndef REEBER_PRUNE_H
#define REEBER_PRUNE_H

#include <array>
#include "merge-tree-block.h"

// Prune internal vertices + vertices that are not minima in the boundary
// (more accurately in the pairwise intersection)
struct PruneInitial
{
    typedef MergeTreeBlock::OffsetGrid          OffsetGrid;

            PruneInitial(MergeTreeBlock* b_, const OffsetGrid& g_):
                b(b_), g(g_)        {}

    bool    operator()(MergeTreeBlock::Index v) const
    {
        if (b->local.internal_test()(v))
            return true;

        MergeTreeBlock::Vertex  p    = b->local.position(v);
        std::array<int, 3> side {{ 0, 0, 0 }};
        for (int i = 0; i < 3; ++i)
            if (p[i] == b->local.from()[i])
                side[i] = -1;
            else if (p[i] == b->local.to()[i])
                side[i] = 1;

        int zeroes = 0;
        MergeTreeBlock::Box     side_box = b->local;
        for (int i = 0; i < 3; ++i)
        {
            if (side[i] == -1)
                side_box.to()[i] = side_box.from()[i];
            else if (side[i] == 1)
                side_box.from()[i] = side_box.to()[i];
            else // (side[i] == 0)
                ++zeroes;
        }
        if (zeroes < 2)     // corner
            return false;

        typedef     MergeTreeBlock::MergeTree::Node::ValueVertex    ValueVertex;
        ValueVertex vval(g(v), v);
        for(MergeTreeBlock::Index u : side_box.link(v))
        {
            ValueVertex uval(g(u), u);
            if (b->mt.cmp(uval, vval))      // v is not a minimum
                return true;
        }

        return false;
    }

    MergeTreeBlock*     b;
    const OffsetGrid&   g;
};

#endif
