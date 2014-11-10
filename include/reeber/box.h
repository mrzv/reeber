#ifndef REEBER_BOX_H
#define REEBER_BOX_H

#include <boost/range/adaptor/filtered.hpp>
#include <boost/range/adaptor/transformed.hpp>

#include "grid.h"
#include "vertices.h"

namespace reeber
{

namespace ba = boost::adaptors;

template<unsigned D>
class Box
{
    public:
        typedef             GridRef<void*, D>                                       GridProxy;
        typedef             typename GridProxy::Vertex                              Position;

        struct              InternalTest;
        struct              BoundsTest;
        struct              PositionToVertex;

        class               FreudenthalLinkIterator;
        typedef             boost::iterator_range<FreudenthalLinkIterator>          FreudenthalLinkRange;

        typedef             VerticesIterator<Position>                              VI;
        typedef             boost::transformed_range
                                <PositionToVertex, const boost::iterator_range<VI> >    VertexRange;

        // Topology interface
        typedef             typename GridProxy::Index                               Vertex;
        typedef             boost::transformed_range
                                <PositionToVertex,
                                 const boost::filtered_range
                                     <BoundsTest,
                                      const FreudenthalLinkRange> >                 Link;


                            Box()                                                   {}
                            Box(const Position& shape):
                                g_(0, shape), to_(shape - Position::one())          {}
                            Box(const Position& shape,
                                const Position& from,
                                const Position& to):
                                g_(0, shape), from_(from), to_(to)                  {}


        const Position&     from() const                                            { return from_; }
        const Position&     to() const                                              { return to_; }
        Position&           from()                                                  { return from_; }
        Position&           to()                                                    { return to_; }
        Position            shape() const                                           { return to_ - from_ + Position::one(); }
        const Position&     grid_shape() const                                      { return g_.shape(); }
        static unsigned     dimension()                                             { return D; }

        size_t              size() const                                            { size_t c = 1; for (unsigned i = 0; i < D; ++i) c *= (to_[i] - from_[i] + 1); return c; }

        VertexRange         vertices() const                                        { return boost::iterator_range<VI>(VI::begin(from_, to_), VI::end(from_, to_))
                                                                                                | ba::transformed(PositionToVertex(*this)); }
        Link                link(const Position& p) const                           { return FreudenthalLinkRange(FreudenthalLinkIterator::begin(p), FreudenthalLinkIterator::end(p))
                                                                                                | ba::filtered(BoundsTest(*this))
                                                                                                | ba::transformed(PositionToVertex(*this)); }
        Link                link(const Vertex& v) const                             { return link(g_.vertex(v)); }

        Box                 intersect(const Box& other) const;
        bool                intersects(const Box& other) const;
        void                merge(const Box& other);

        bool                contains(const Position& p) const;
        bool                contains(const Vertex& v) const                         { return contains(g_.vertex(v)); }

        bool                boundary(const Position& p, bool degenerate = false) const;
        bool                boundary(const Vertex& v, bool deg = false) const       { return boundary(g_.vertex(v), deg); }
        Box                 side(unsigned axis, bool upper) const;

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
            Vertex          operator()(const Position& p) const                     { return box_.g_.index(p); }
            const Box&      box_;
        };

        // debug
        Position            position(const Vertex& v) const                         { return g_.vertex(v); }

    private:
        GridProxy           g_;
        Position            from_, to_;
};

}

#include "box.hpp"

#endif
