#ifndef REEBER_GHOSTS_H
#define REEBER_GHOSTS_H

#include <diy/link.hpp>

#include <reeber/box.h>
#include <reeber/grid.h>

unsigned spread_bits(unsigned x, unsigned factor)
{
    unsigned res = 0;
    for (unsigned j = 0; (1 << j) <= x; ++j)
        if (x & (1 << j))
            res |= (1 << factor*j);     // spread the bits into even positions
    return res;
}

// send ghosts to the lower neighbors
template<class Block_>
struct EnqueueGhosts
{
    typedef     Block_                              Block;
    typedef     typename Block::OffsetGrid          Grid;
    typedef     typename Block::Box                 Box;
    typedef     Grid                                Block::*GridPtr;
    typedef     Box                                 Block::*BoxPtr;

                EnqueueGhosts(GridPtr grid_, BoxPtr local_):
                    grid(grid_), local(local_)          {}

    void        operator()(void* b_, const diy::Master::ProxyWithLink& cp, void*) const
    {
        typedef     diy::RegularGridLink                        RGLink;
        typedef     typename reeber::RestrictGrid<Grid>::type   GridRestriction;

        Block*      b = static_cast<Block*>(b_);
        RGLink*     l = static_cast<RGLink*>(cp.link());

        // enqueue to lower side
        for (unsigned i = 0; i < 8; ++i)
        {
            unsigned side = spread_bits(i, 2);      // spread the bits into even positions

            int nbr = l->direction(diy::Direction(side));
            if (nbr == -1)
                continue;

            GridRestriction grid_side = GridRestriction::side(b->*grid, side);
            for (unsigned i = 0; i < 3; ++i)
                if (grid_side.to()[i] != grid_side.from()[i] &&
                    (b->*local).to()[i]  != (b->*local).grid_shape()[i]-1)       // reduce the grid sides by one
                    grid_side.to()[i]--;
            cp.enqueue(l->target(nbr), grid_side);
        }
    }

    GridPtr grid;
    BoxPtr  local;
};

// receive ghosts from the upper neighbors
template<class Block_>
struct DequeueGhosts
{
    typedef     Block_                              Block;
    typedef     typename Block::OffsetGrid          Grid;
    typedef     typename Block::Box                 Box;
    typedef     Grid                                Block::*GridPtr;
    typedef     Box                                 Block::*BoxPtr;


                DequeueGhosts(GridPtr grid_, BoxPtr local_):
                    grid(grid_), local(local_)          {}

    void        operator()(void* b_, const diy::Master::ProxyWithLink& cp, void*) const
    {
        typedef     diy::RegularGridLink                        RGLink;
        typedef     typename reeber::RestrictGrid<Grid>::type   GridRestriction;

        Block*          b = static_cast<Block*>(b_);
        RGLink*         l = static_cast<RGLink*>(cp.link());

        // dequeue from upper sides
        for (unsigned i = 0; i < 8; ++i)
        {
            unsigned side = spread_bits(i, 2) << 1;      // spread the bits into odd positions

            int nbr = l->direction(diy::Direction(side));
            if (nbr == -1)
                continue;

            GridRestriction grid_side = GridRestriction::side(b->*grid, side);
            for (unsigned i = 0; i < 3; ++i)
                if (grid_side.to()[i] != grid_side.from()[i] &&
                    (b->*local).to()[i]  != (b->*local).grid_shape()[i]-1)       // reduce the grid sides by one
                    grid_side.to()[i]--;
            cp.dequeue(l->target(nbr).gid, grid_side);
        }
    }

    GridPtr grid;
    BoxPtr  local;
};

#endif
