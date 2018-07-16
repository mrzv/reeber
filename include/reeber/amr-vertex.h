#pragma once

#include <assert.h>
#include <stdexcept>
#include <boost/functional/hash.hpp>

#include <diy/types.hpp>
#include <diy/point.hpp>
#include <diy/serialization.hpp>

namespace reeber {

    struct AmrVertexId
    {
        int gid;
        size_t vertex;

        AmrVertexId() :
                gid(-1), vertex(-1)
        {}

        AmrVertexId(int _gid, size_t _vertex) :
                gid(_gid), vertex(_vertex)
        {}

        // comparison - lexicographic

        bool operator==(const AmrVertexId& other) const
        {
            return std::tie(gid, vertex) == std::tie(other.gid, other.vertex);
        }

        bool operator!=(const AmrVertexId& other) const
        {
            return !(*this == other);
        }

        bool operator<(const AmrVertexId& other) const
        {
            return std::tie(vertex, gid) < std::tie(other.vertex, other.gid);
        }

        bool operator>(const AmrVertexId& other) const
        {
            return other < *this;
        }

        bool operator<=(const AmrVertexId& other) const
        {
            return !(*this > other);
        }


        bool operator>=(const AmrVertexId& other) const
        {
            return !(*this < other);
        }

        friend std::ostream& operator<<(std::ostream& os, const AmrVertexId& a)
        {
            os << "(gid = " << a.gid << ", index = " << a.vertex << ")";
            return os;
        }

        // implicit conversion
        operator size_t() const
        {
            return vertex;
        }
    };

    using AmrEdge = std::tuple<AmrVertexId, AmrVertexId>;

    AmrEdge reverse_amr_edge(const AmrEdge& e)
    {
        return AmrEdge { std::get<1>(e), std::get<0>(e) };
    }

    using AmrEdgeContainer = std::vector<AmrEdge>;

} // namespace reeber


namespace std {

    //    template<>
    //    struct hash<diy::BlockID>
    //    {
    //        std::size_t operator()(const diy::BlockID gid)
    //        {
    //            std::size_t seed = 0;
    //            boost::hash_combine(seed, gid);
    //            return seed;
    //        }
    //    };

    template<>
    struct hash<reeber::AmrVertexId>
    {
        std::size_t operator()(const reeber::AmrVertexId& id) const
        {
            std::size_t seed = 0;
            boost::hash_combine(seed, id.gid);
            boost::hash_combine(seed, id.vertex);
            return seed;
        }
    };

    template<>
    struct hash<reeber::AmrEdge>
    {
        std::size_t operator()(const reeber::AmrEdge& e) const
        {
            std::size_t seed = 0;

            boost::hash_combine(seed, std::get<0>(e).gid);
            boost::hash_combine(seed, std::get<0>(e).vertex);

            boost::hash_combine(seed, std::get<1>(e).gid);
            boost::hash_combine(seed, std::get<1>(e).vertex);

            return seed;
        }
    };

    std::ostream& operator<<(std::ostream& os, const reeber::AmrEdge& e)
    {
        os << "(" << std::get<0>(e) << " <-> " << std::get<1>(e) << ")";
        return os;
    }
}


