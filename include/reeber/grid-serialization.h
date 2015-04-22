#ifndef REEBER_GRID_SERIALIZATION_H
#define REEBER_GRID_SERIALIZATION_H

#include "grid.h"
#include "vertices.h"
#include <diy/serialization.hpp>

namespace diy
{
    template<class C, unsigned D>
    struct Serialization< ::reeber::Grid<C,D> >
    {
        typedef     ::reeber::Grid<C,D>     Grid;
        typedef     typename Grid::Vertex   Vertex;

        static void save(BinaryBuffer& bb, const Grid& g)
        {
            diy::save(bb, g.shape());
            diy::save(bb, g.data(), g.size());
        }

        static void load(BinaryBuffer& bb, Grid& g)
        {
            Vertex shape;
            diy::load(bb, shape);
            Grid tmp(shape);
            diy::load(bb, tmp.data(), tmp.size());
            g.swap(tmp);
        }
    };

    template<class C, unsigned D>
    struct Serialization< ::reeber::OffsetGrid<C,D> >
    {
        typedef     ::reeber::OffsetGrid<C,D>       OffsetGrid;
        typedef     typename OffsetGrid::Grid       Grid;
        typedef     typename OffsetGrid::GridProxy  GridProxy;
        typedef     typename OffsetGrid::Vertex     Vertex;

        static void save(BinaryBuffer& bb, const OffsetGrid& g)
        {
            diy::save(bb, g.g_.shape());
            diy::save(bb, g.offset);
            diy::save(bb, static_cast<const Grid&>(g));
        }

        static void load(BinaryBuffer& bb, OffsetGrid& g)
        {
            Vertex full_shape;
            diy::load(bb, full_shape);
            GridProxy(0, full_shape).swap(g.g_);
            diy::load(bb, g.offset);
            diy::load(bb, static_cast<Grid&>(g));
        }
    };

    // Deliberately doesn't save/load from and to, so that one side could be mapped into another side
    template<class C, unsigned D>
    struct Serialization< ::reeber::GridRestriction<C,D> >
    {
        typedef     ::reeber::GridRestriction<C,D>                  GridRestriction;
        typedef     typename GridRestriction::Grid                  Grid;
        typedef     typename Grid::Vertex                           Vertex;
        typedef     ::reeber::VerticesIterator<Vertex>              VerticesIterator;

        static void save(BinaryBuffer& bb, const GridRestriction& gr)
        {
            VerticesIterator vi  = VerticesIterator::begin(gr.from(), gr.to()),
                             end = VerticesIterator::end(gr.from(), gr.to());
            while (vi != end)
                diy::save(bb, gr.grid()(*vi++));
        }

        static void load(BinaryBuffer& bb, GridRestriction& gr)
        {
            VerticesIterator vi  = VerticesIterator::begin(gr.from(), gr.to()),
                             end = VerticesIterator::end(gr.from(), gr.to());
            while (vi != end)
                diy::load(bb, gr.grid()(*vi++));
        }
    };
}

#endif
