#ifndef REEBER_MERGE_TREE_H
#define REEBER_MERGE_TREE_H

#include <vector>
#include <unordered_map>
#include <tuple>

#include <boost/range/adaptor/map.hpp>
#include <boost/range/algorithm.hpp>

#include "serialization.h"
#include "format.h"

namespace reeber
{

namespace ba = boost::adaptors;

template<class Vertex_, class Value_>
struct MergeTreeNode
{
    typedef                     Vertex_                         Vertex;
    typedef                     Value_                          Value;

    typedef                     MergeTreeNode*                  Neighbor;

    typedef                     std::tuple<Neighbor, Neighbor> Parent;

    bool                        operator<(const MergeTreeNode& other) const     { return value < other.value || (value == other.value && vertex < other.vertex); }
    bool                        operator>(const MergeTreeNode& other) const     { return value > other.value || (value == other.value && vertex > other.vertex); }

    Vertex                      vertex;
    Value                       value;
    Parent                      parent;
    Neighbor                    cur_deepest;
};

template<class Vertex_, class Value_>
class MergeTree
{
    public:
        typedef     Vertex_                             Vertex;
        typedef     Value_                              Value;

        typedef     MergeTreeNode<Vertex,Value>         Node;
        typedef     typename Node::Neighbor             Neighbor;

        typedef     std::unordered_map<Vertex, Neighbor>
                                                        VertexNeighborMap;

    public:
                    MergeTree(bool negate = false):
                        negate_(negate)                 {}

        void        repair(const Neighbor u);
        Neighbor    add(const Vertex& x, Value v);
        void        link(const Neighbor u, const Neighbor s, const Neighbor v)
                                                        { u->parent = std::make_tuple(s, v); }

        Neighbor    find_deepest(const Neighbor u);

        size_t      size() const                        { return nodes_.size(); }

        bool        contains(const Vertex& x) const     { return nodes_.find(x) != nodes_.end(); }

        void        swap(MergeTree& other)              { std::swap(negate_, other.negate_); nodes_.swap(other.nodes_); }

        bool        negate() const                      { return negate_; }
        void        set_negate(bool negate)             { negate_ = negate; }

        template<class T>
        bool        cmp(const T& x, const T& y) const   { return negate_ ? x > y : x < y; }
        bool        cmp(Neighbor x, Neighbor y) const   { return cmp(*x, *y); }

        Neighbor    node(const Vertex& x) const         { return nodes_.find(x)->second; }

        const VertexNeighborMap& nodes() const          { return nodes_; }
        const std::vector<Vertex> vertices() const      { std::vector<Vertex> vs; boost::copy(nodes_ | ba::map_keys, std::back_inserter(vs)); return vs; }

        friend struct ::reeber::Serialization<MergeTree>;

    private:
        VertexNeighborMap& nodes()                      { return nodes_; }

        template<class MT, class T, class F>
        friend void
        compute_merge_tree(MT& mt, const T& t, const F& f);

        template<class MT, class T, class F>
        friend void
        compute_merge_tree2(MT& mt, const T& t, const F& f);

        template<class MT, class F>
        friend void
        traverse_persistence(const MT& mt, const F& f);

        template<class MT, class S>
        friend void
        sparsify(MT& out, const MT& in, const S& special);

        template<class MT>
        friend void
        merge(MT& mt, typename MT::Neighbor u, typename MT::Neighbor s, typename MT::Neighbor v);

        template<class MT, class E>
        friend void
        merge(MT& mt1, const MT& mt2, const E& edges);

    private:
        bool                        negate_;
        VertexNeighborMap           nodes_;
};

/**
 * Topology defines a range vertices() and a link(v) function;
 *          vertices should be allowed to repeat (will simplify uniting multiple trees).
 */
template<class MergeTree, class Topology, class Function>
void compute_merge_tree(MergeTree& mt, const Topology& topology, const Function& f);

template<class MergeTree, class Topology, class Function>
void compute_merge_tree2(MergeTree& mt, const Topology& topology, const Function& f);

template<class MergeTree, class Functor>
void traverse_persistence(const MergeTree& mt, const Functor& f);

template<class MergeTree, class Special>
void sparsify(MergeTree& out, const MergeTree& in, const Special& special);

template<class MergeTree>
void merge(MergeTree& mt, typename MergeTree::Neighbor u, typename MergeTree::Neighbor s, typename MergeTree::Neighbor v);

template<class MergeTree, class Edges>
void merge(MergeTree& mt1, const MergeTree& mt2, const Edges& edges);

}

#include "triplet-merge-tree.hpp"

#endif
