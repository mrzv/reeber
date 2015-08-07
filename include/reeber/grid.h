#ifndef REEBER_GRID_H
#define REEBER_GRID_H

#include "point.h"

namespace reeber
{

template<class C, unsigned D>
struct Grid;

template<class C, unsigned D>
struct GridRef
{
    public:
        typedef     C                                           Value;

        typedef     Point<int, D>                               Vertex;
        typedef     size_t                                      Index;

    public:
        template<class Int>
                GridRef(C* data, const Point<Int,D>& shape, bool c_order = true):
                    data_(data), shape_(shape), c_order_(c_order)   { set_stride(); }

                GridRef(Grid<C,D>& g):
                    data_(g.data()), shape_(g.shape()),
                    c_order_(g.c_order())                       { set_stride(); }

        template<class Int>
        C       operator()(const Point<Int, D>& v) const        { return (*this)(index(v)); }

        template<class Int>
        C&      operator()(const Point<Int, D>& v)              { return (*this)(index(v)); }

        C       operator()(Index i) const                       { return data_[i]; }
        C&      operator()(Index i)                             { return data_[i]; }

        const Vertex&
                shape() const                                   { return shape_; }

        const C*
                data() const                                    { return data_; }
        C*      data()                                          { return data_; }

        // Set every element to the given value
        GridRef&    operator=(C value)                          { Index s = size(); for (Index i = 0; i < s; ++i) data_[i] = value; return *this; }
        GridRef&    operator/=(C value)                         { Index s = size(); for (Index i = 0; i < s; ++i) data_[i] /= value; return *this; }

        Vertex      vertex(Index idx) const                     { Vertex v; for (unsigned i = 0; i < D; ++i) { v[i] = idx / stride_[i]; idx %= stride_[i]; } return v; }
        Index       index(const Vertex& v) const                { Index idx = 0; for (unsigned i = 0; i < D; ++i) { idx += ((Index) v[i]) * ((Index) stride_[i]); } return idx; }

        Index       size() const                                { return size(shape()); }
        void        swap(GridRef& other)                        { std::swap(data_, other.data_); std::swap(shape_, other.shape_); std::swap(stride_, other.stride_); std::swap(c_order_, other.c_order_); }

        bool        c_order() const                             { return c_order_; }

        static
        unsigned    dimension()                                 { return D; }

    protected:
        static Index
                size(const Vertex& v)                           { Index res = 1; for (unsigned i = 0; i < D; ++i) res *= v[i]; return res; }

        void    set_stride()
        {
            Index cur = 1;
            if (c_order_)
                for (unsigned i = D; i > 0; --i) { stride_[i-1] = cur; cur *= shape_[i-1]; }
            else
                for (unsigned i = 0; i < D; ++i) { stride_[i] = cur; cur *= shape_[i]; }

        }
        void    set_shape(const Vertex& v)                      { shape_ = v; set_stride(); }
        void    set_data(C* data)                               { data_ = data; }
        void    set_c_order(bool order)                         { c_order_ = order; }

    private:
        C*      data_;
        Vertex  shape_;
        Vertex  stride_;
        bool    c_order_;
};


template<class C, unsigned D>
struct Grid: public GridRef<C,D>
{
    public:
        typedef     GridRef<C,D>                                Parent;
        typedef     typename Parent::Value                      Value;
        typedef     typename Parent::Index                      Index;
        typedef     typename Parent::Vertex                     Vertex;
        typedef     Parent                                      Reference;

        template<class U>
        struct rebind { typedef Grid<U,D>                       type; };

    public:
                Grid():
                    Parent(new C[0], Vertex::zero())            {}
        template<class Int>
                Grid(const Point<Int, D>& shape, bool c_order = true):
                    Parent(new C[size(shape)], shape, c_order)
                {}

                Grid(const Parent& g):
                    Parent(new C[size(g.shape())], g.shape(),
                           g.c_order())                     { copy_data(g.data()); }

        template<class OtherGrid>
                Grid(const OtherGrid& g):
                    Parent(new C[size(g.shape())],
                           g.shape(),
                           g.c_order())   { copy_data(g.data()); }

                ~Grid()                                         { delete[] Parent::data(); }

        template<class OC>
        Grid&   operator=(const GridRef<OC, D>& other)
        {
            delete[] Parent::data();
            Parent::set_c_order(other.c_order());       // NB: order needs to be set before the shape, to set the stride correctly
            Parent::set_shape(other.shape());
            Index s = size(shape());
            Parent::set_data(new C[s]);
            copy_data(other.data());
            return *this;
        }

        using Parent::data;
        using Parent::shape;
        using Parent::operator();
        using Parent::operator=;
        using Parent::size;

    private:
        template<class OC>
        void    copy_data(const OC* data)
        {
            Index s = size(shape());
            for (Index i = 0; i < s; ++i)
                Parent::data()[i] = data[i];
        }
};


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
