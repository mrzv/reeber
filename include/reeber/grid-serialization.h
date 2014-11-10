#ifndef REEBER_GRID_SERIALIZATION_H
#define REEBER_GRID_SERIALIZATION_H

#include "grid.h"
#include <diy/serialization.hpp>

namespace diy
{
    template<class C, unsigned D>
    struct Serialization< ::reeber::Grid<C,D> >
    {
        typedef     ::reeber::Grid<C,D>        Grid;
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
}

#endif
