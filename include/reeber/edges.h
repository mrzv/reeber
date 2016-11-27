#pragma once

namespace reeber
{

namespace detail
{
    template<class Vertex>
    using Edge  = std::tuple<Vertex, Vertex>;

    template<class Vertex>
    struct edge_hash
    {
        std::size_t operator()(const Edge<Vertex>& k) const
        {
            size_t x = std::hash<Vertex>()(std::get<0>(k));
            size_t y = std::hash<Vertex>()(std::get<1>(k));
            return x ^ (y + 0x9e3779b9 + (x<<6) + (x>>2));
        }
    };

    template<class Vertex>
    struct edge_equal
    {
        bool operator()(const Edge<Vertex>& v0, const Edge<Vertex>& v1) const
        {
            return std::get<0>(v0) == std::get<0>(v1) && std::get<1>(v0) == std::get<1>(v1);
        }
    };
}

template<class Vertex, class Value>
using EdgeMap = std::unordered_map<detail::Edge<Vertex>, std::tuple<Value, Vertex>, detail::edge_hash<Vertex>, detail::edge_equal<Vertex>>;

template<class Vertex, class Value>
using EdgeMaps = std::unordered_map<int, EdgeMap<Vertex,Value>>;

}
