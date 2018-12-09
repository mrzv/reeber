#pragma once

#include <vector>
#include <unordered_map>

#include <format.h>

#include <diy/serialization.hpp>

template<class Vertex_>
struct DisjointSets
{
    using Vertex = Vertex_;
    using VertexVertexMap = std::unordered_map<Vertex, Vertex>;
    using VertexSizeMap = std::unordered_map<Vertex, int>;

    VertexVertexMap parent_;
    VertexSizeMap size_;

    DisjointSets() {}

    template<class VertexContainer>
    DisjointSets(const VertexContainer& vc)
    {
        for(const auto &v : vc)
        {
            create_component(v);
        }
    }

    void make_component(const Vertex& v)
    {
        assert(parent_.count(v) == 0);
        parent_[v] = v;
        size_[v] = 1;
    }

    bool has_component(const Vertex& v) const
    {
        return parent_.count(v) == 1;
    }

    bool are_connected(const Vertex& a, const Vertex_& b)
    {
        return find_component(a) == find_component(a);
    }

    Vertex find_component(const Vertex& a)
    {
        Vertex x = a;
        bool debug = false;
        while (parent_.at(x) != x)
        {
            auto next = parent_[x];
            if (debug) fmt::print("in find_component_in_disjoint_sets, x = {}, next = {}\n", x, next);
            parent_[x] = parent_.at(next);
            x = next;
        }
        if (debug) fmt::print("Exiting find_component_in_disjoint_sets, x = {}\n", x);
        return x;
    }

    Vertex unite_components_by_roots(Vertex a_root, Vertex b_root)
    {
        bool debug = false;
        if (a_root == b_root)
            return a_root;
        if (size_[a_root] < size_[b_root])
            std::swap(a_root, b_root);
        parent_[b_root] = a_root;
        size_[a_root] += size_[b_root];
        if (debug) fmt::print("Successfully united roots {} and {}\n", a_root, b_root);
        return a_root;
    }

    void unite_components(const Vertex& a, const Vertex& b)
    {
        auto a_root = find_component(a);
        auto b_root = find_component(b);
        unite_components_by_roots(a_root, b_root);
    }

    void disjoint_union(const DisjointSets<Vertex>& other)
    {
        for(const auto& key_val_pair : other.parent_)
        {
            assert(parent_.count(key_val_pair.first) == 0);
        }
        parent_.insert(other.parent_.begin(), other.parent_.end());
        size_.insert(other.size_.begin(), other.size_.end());
    }

};


template <>
template<class Vertex>
struct diy::Serialization<DisjointSets<Vertex>>
{
    using DS = DisjointSets<Vertex>;

    static void save(BinaryBuffer& bb, const DS& disjoint_sets)
    {
        diy::save(bb, disjoint_sets.parent_);
        diy::save(bb, disjoint_sets.size_);

    }

    static void load(BinaryBuffer& bb, const DS& disjoint_sets)
    {
        diy::load(bb, disjoint_sets.parent_);
        diy::load(bb, disjoint_sets.size_);

    }
};
