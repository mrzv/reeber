#ifndef REEBER_BOX_H
#define REEBER_BOX_H

#include "range/filtered.h"
#include "range/transformed.h"
#include "range/utility.h"

#include "grid.h"
#include "vertices.h"

namespace reeber
{

template<unsigned D>
class Box
{
    public:
        typedef             GridRef<void*, D>                                       GridProxy;
        typedef             typename GridProxy::Vertex                              Position;

        struct              InternalTest;
        struct              BoundaryTest;
        struct              BoundsTest;
        struct              PositionToVertex;

        class               FreudenthalLinkIterator;
        typedef             range::iterator_range<FreudenthalLinkIterator>          FreudenthalLinkRange;

        typedef             VerticesIterator<Position>                              VI;
        typedef             range::transformed_range
                                <range::iterator_range<VI>, PositionToVertex>       VertexRange;

        // Topology interface
        typedef             typename GridProxy::Index                               Vertex;
        typedef             range::filtered_range
                                     <FreudenthalLinkRange, BoundsTest>             PositionLink;
        typedef             range::transformed_range
                                <PositionLink, PositionToVertex>                    Link;


                            Box(): g_(0, Position())                                {}
                            Box(const Position& shape):
                                g_(0, shape), to_(shape - Position::one())          {}
                            Box(const Position& shape,
                                const Position& from,
                                const Position& to):
                                g_(0, shape), from_(from), to_(to)                  {}
                            Box(const Position& from,
                                const Position& to):
                                g_(0, to - from + Position::one()),
                                from_(from), to_(to)                                {}


        const Position&     from() const                                            { return from_; }
        const Position&     to() const                                              { return to_; }
        Position&           from()                                                  { return from_; }
        Position&           to()                                                    { return to_; }
        Position            shape() const                                           { return to_ - from_ + Position::one(); }
        const Position&     grid_shape() const                                      { return g_.shape(); }
        static unsigned     dimension()                                             { return D; }

        size_t              size() const                                            { size_t c = 1; for (unsigned i = 0; i < D; ++i) c *= (to_[i] - from_[i] + 1); return c; }

        VertexRange         vertices() const                                        { return range::iterator_range<VI>(VI::begin(from_, to_), VI::end(from_, to_))
                                                                                                | range::transformed(position_to_vertex()); }

        decltype(auto)      positions() const                                       { return range::iterator_range<VI>(VI::begin(from_, to_), VI::end(from_, to_));}

        Link                link(const Position& p) const                           { return position_link(p)
                                                                                                | range::transformed(position_to_vertex()); }
        Link                link(const Vertex& v) const                             { return link(position(v)); }

        PositionLink        position_link(const Position& p) const                  { return FreudenthalLinkRange(FreudenthalLinkIterator::begin(p), FreudenthalLinkIterator::end(p))
                                                                                                | range::filtered(bounds_test()); }
        PositionLink        position_link(const Vertex& v) const                    { return position_link(position(v)); }

        Box                 intersect(const Box& other) const;
        bool                intersects(const Box& other) const;
        void                merge(const Box& other);

        bool                contains(const Position& p) const;
        bool                contains(const Vertex& v) const                         { return contains(position(v)); }

        bool                boundary(const Position& p, bool degenerate = false) const;
        bool                boundary(const Vertex& v, bool deg = false) const       { return boundary(position(v), deg); }
        Box                 side(unsigned axis, bool upper) const;

        InternalTest        internal_test() const                                   { return InternalTest(*this); }
        BoundaryTest        boundary_test() const                                   { return BoundaryTest(*this); }
        BoundsTest          bounds_test() const                                     { return BoundsTest(*this); }
        PositionToVertex    position_to_vertex() const                              { return PositionToVertex(*this); }


        void                swap(Box& other)                                        { g_.swap(other.g_); std::swap(from_, other.from_); std::swap(to_, other.to_); }

        bool                operator==(const Box& other) const                      { return from_ == other.from_ && to_ == other.to_; }

        template<class C_, class T_>
        friend std::basic_ostream<C_,T_>&
        operator<<(std::basic_ostream<C_,T_>& out, const Box& b)                    { out << "Box: " << b.from_ << " - " << b.to_ << " inside " << b.g_.shape(); return out; }

        struct InternalTest
        {
                            InternalTest(const Box& box): box_(box)                 {}
            bool            operator()(const Vertex& v) const                       { return !box_.boundary(v); }
            const Box&      box_;
        };

        struct BoundaryTest
        {
                            BoundaryTest(const Box& box): box_(box)                 {}
            bool            operator()(const Vertex& v) const                       { return box_.boundary(v); }
            bool            operator()(const Position& p) const                     { return box_.boundary(box_.position_to_vertex()(p)); }
            const Box&      box_;
        };

        struct BoundsTest
        {
                            BoundsTest(const Box& box): box_(box)                   {}
            bool            operator()(const Position& p) const                     { return box_.contains(p); }
            bool            operator()(const Vertex& v) const                       { return box_.contains(v); }
            const Box&      box_;
        };

        struct PositionToVertex
        {
            typedef         Vertex                                                  result_type;
                            PositionToVertex(const Box& box): box_(box)             {}
            Vertex          operator()(Position p) const                            { for (unsigned i = 0; i < D; ++i) p[i] %= box_.grid_shape()[i]; return box_.g_.index(p); }
            const Box&      box_;
        };

        // computes position inside the box (adjusted for the wrap-around, if need be)
        Position            position(const Vertex& v) const                         { Position p = g_.vertex(v); for (unsigned i = 0; i < D; ++i) if (p[i] < from()[i]) p[i] += grid_shape()[i]; return p; }

        Position            positive_position(Position p) const                     { for (unsigned i = 0; i < D; ++i) if (p[i] < 0) p[i] += grid_shape()[i]; return p; }

    private:
        GridProxy           g_;
        Position            from_, to_;
};

}

#include "box.hpp"

#endif
