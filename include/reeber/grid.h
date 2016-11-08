#ifndef REEBER_GRID_H
#define REEBER_GRID_H

#include "point.h"
#include <diy/grid.hpp>

namespace reeber
{

template<class C, unsigned D>
using Grid = ::diy::Grid<C,D>;

template<class C, unsigned D>
using GridRef = ::diy::GridRef<C,D>;

template<class C, unsigned D>
struct OffsetGrid: public Grid<C, D>
{
    typedef         ::reeber::Grid<C,D>                     Grid;
    typedef         typename Grid::Value                    Value;
    typedef         typename Grid::Vertex                   Vertex;
    typedef         typename Grid::Index                    Index;
    typedef         GridRef<void*, 3>                       GridProxy;      // used for translation operations on the full grid

                    OffsetGrid():
                        Grid(Vertex::zero()), g_(0, Vertex::zero()), offset(Vertex::zero()) {}

                    OffsetGrid(const Vertex& full_shape, const Vertex& from, const Vertex& to, bool c_order = true):
                        Grid(to - from + Vertex::one(), c_order),
                        g_(0, full_shape),
                        offset(from)                        {}

    // These operations take global indices as input and translate them into the local values
    template<class Int>
    Value           operator()(const Point<Int, D>& v) const            { return Grid::operator()(local(v)); }

    template<class Int>
    Value&          operator()(const Point<Int, D>& v)                  { return Grid::operator()(local(v)); }

    Value           operator()(Index i) const                           { return Grid::operator()(local(g_.vertex(i))); }
    Value&          operator()(Index i)                                 { return Grid::operator()(local(g_.vertex(i))); }

    Vertex          local(Vertex v) const                               { v -= offset; for (unsigned i = 0; i < D; ++i) if (v[i] < 0) v[i] += g_.shape()[i]; return v; }

    void            swap(OffsetGrid& other)                             { Grid::swap(other); std::swap(g_, other.g_); std::swap(offset, other.offset); }

    GridProxy       g_;
    Vertex          offset;
};

template<class C, unsigned D>
struct GridRestriction
{
    typedef         GridRef<C,D>                            Grid;
    typedef         typename Grid::Vertex                   Vertex;

                    GridRestriction(Grid& grid, const Vertex& from, const Vertex& to):
                        grid_(&grid), from_(from), to_(to)   {}

    static
    GridRestriction side(Grid& g, unsigned s)
    {
        Vertex from = Vertex::zero();
        Vertex to   = g.shape() - Vertex::one();

        unsigned dim = 0;
        while (s != 0)
        {
            if (s & 1)                      // lower side
                to[dim] = from[dim];
            else if (s & 2)                 // upper side
                from[dim] = to[dim];

            s >>= 2;
            ++dim;
        }

        return GridRestriction(g, from, to);
    }

    const Vertex&   from() const                            { return from_; }
    Vertex&         from()                                  { return from_; }
    const Vertex&   to() const                              { return to_; }
    Vertex&         to()                                    { return to_; }

    const Grid&     grid() const                            { return *grid_; }
    Grid&           grid()                                  { return *grid_; }

    Grid*           grid_;
    Vertex          from_, to_;
};

template<class Grid>
struct RestrictGrid;

template<class C, unsigned D>
struct RestrictGrid< Grid<C,D> >
{
    typedef     GridRestriction<C,D>        type;
};

template<class C, unsigned D>
struct RestrictGrid< OffsetGrid<C,D> >
{
    typedef     GridRestriction<C,D>        type;
};

}

#endif
