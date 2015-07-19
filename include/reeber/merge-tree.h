#ifndef REEBER_MERGE_TREE_H
#define REEBER_MERGE_TREE_H

#include <vector>
#include <map>

#include <boost/foreach.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/range/adaptor/map.hpp>

#include "serialization.h"

namespace reeber
{

namespace ba = boost::adaptors;

template<class Vertex_, class Value_>
struct MergeTreeNode
{
    typedef                     Vertex_                         Vertex;
    typedef                     Value_                          Value;

    typedef                     std::pair<Value, Vertex>        ValueVertex;
    typedef                     std::vector<ValueVertex>        VerticesVector;

    typedef                     MergeTreeNode*                  Neighbor;
    typedef                     std::vector<Neighbor>           Neighbors;

    bool                        operator<(const MergeTreeNode& other) const     { return value < other.value || (value == other.value && vertex < other.vertex); }
    bool                        operator>(const MergeTreeNode& other) const     { return value > other.value || (value == other.value && vertex > other.vertex); }

    template<class F>
    bool                        any_vertex(const F& f) const                    { if (f(vertex)) return true; BOOST_FOREACH(const ValueVertex& vv, vertices) if (f(vv.second)) return true; return false; }

    Neighbor                    parent;
    Neighbors                   children;

    Vertex                      vertex;
    Value                       value;
    VerticesVector              vertices;       // degree-2 vertices on the path from this node to its parent

    void*                       aux;
};

template<class Vertex_, class Value_>
class MergeTree
{
    public:
        typedef     Vertex_                             Vertex;
        typedef     Value_                              Value;

        typedef     MergeTreeNode<Vertex,Value>         Node;
        typedef     typename Node::Neighbor             Neighbor;

        typedef     std::map<Vertex, Neighbor>          VertexNeighborMap;               // TODO: change something like unordered_map; possibly a template parameter

    public:
                    MergeTree(bool negate = false):
                        negate_(negate)                 {}
                    ~MergeTree()                        { BOOST_FOREACH(typename VertexNeighborMap::value_type& x, nodes_) delete x.second; }

        Neighbor    operator[](const Vertex& x) const   { return nodes_.find(x)->second; }
        Neighbor    add(const Vertex& x, Value v);

        size_t      size() const                        { return nodes_.size(); }

        bool        contains(const Vertex& x) const     { return nodes_.find(x) != nodes_.end(); }

        Neighbor    find(const Vertex& x) const;        // finds root of the subtree containing x
        static void link(Neighbor xn, Neighbor yn);     // links yn as a child of xn

        void        swap(MergeTree& other)              { std::swap(negate_, other.negate_); nodes_.swap(other.nodes_); }

        bool        negate() const                      { return negate_; }
        void        set_negate(bool negate)             { negate_ = negate; }

        template<class T>
        bool        cmp(const T& x, const T& y) const   { return negate_ ? x > y : x < y; }

        const VertexNeighborMap&
                    nodes() const                       { return nodes_; }

        Neighbor    find_root() const                   { BOOST_FOREACH(Neighbor n, nodes() | ba::map_values) if (!n->parent) return n; return 0; }
        size_t      count_roots() const                 { size_t s = 0; BOOST_FOREACH(Neighbor n, nodes() | ba::map_values) if (!n->parent) ++s; return s; }
        void        reset_aux() const                   { BOOST_FOREACH(Neighbor n, nodes() | ba::map_values) aux_neighbor(n) = 0; }

        friend struct ::reeber::Serialization<MergeTree>;

    public:
        struct      CollapseEvent   { static const char* name() { return "collapsed"; } };
        struct      EraseEvent      { static const char* name() { return "erased"; } };
        struct      FindStepEvent   { static const char* name() { return "find-step"; } };

    private:
        VertexNeighborMap&
                    nodes()                             { return nodes_; }

        static
        Neighbor&   aux_neighbor(Neighbor n)            { return (Neighbor&) n->aux; }

        template<class MT, class T, class F, class C>
        friend void
        compute_merge_tree(MT& mt, const T& t, const F& f, const C& c);

        template<class MT, class F>
        friend void
        traverse_persistence(const MT& mt, const F& f);

        template<class MT, class S>
        friend void
        sparsify(MT& mt, const S& s);

        template<class MT, class S>
        friend void
        sparsify(MT& out, const MT& in, const S& s);

        template<class MT>
        friend void
        merge(MT& mt, const std::vector<MT>& trees);

        template<class MT, class P, class S>
        friend void
        remove_degree2(MT& mt, const P& p, const S& s);

    private:
        bool                        negate_;
        VertexNeighborMap           nodes_;
};

/**
 * Topology defines a range vertices() and a link(v) function;
 *          vertices should be allowed to repeat (will simplify uniting multiple trees).
 * Collapsible indicates whether a node can be simplified away (if it's a degree-2 node).
 */
template<class MergeTree, class Topology, class Function, class Collapsible>
void compute_merge_tree(MergeTree& mt, const Topology& topology, const Function& f, const Collapsible& collapsible);

// By default, everything is collapsible
template<class MergeTree, class Topology, class Function>
void compute_merge_tree(MergeTree& mt, const Topology& topology, const Function& f)
{
    compute_merge_tree(mt, topology, f, boost::lambda::constant(true));
}

template<class MergeTree, class Functor>
void traverse_persistence(const MergeTree& mt, const Functor& f);

template<class MergeTree, class Special>
void sparsify(MergeTree& mt, const Special& special);

template<class MergeTree, class Special>
void sparsify(MergeTree& out, const MergeTree& in, const Special& special);

template<class MergeTree>
void merge(MergeTree& mt, const std::vector<MergeTree>& trees);

template<class MergeTree, class Edges>
void merge(MergeTree& mt, const std::vector<MergeTree>& trees, const Edges& edges);

// preserve indicates whether the a degree-2 node should be collapsed or removed
// special indicates whether the node should be kept in place even if it's degree2
template<class MergeTree, class Preserve, class Special>
void remove_degree2(MergeTree& mt, const Preserve& preserve, const Special& special);

template<class MergeTree, class Preserve>
void remove_degree2(MergeTree& mt, const Preserve& preserve)
{
    remove_degree2(mt, preserve, boost::lambda::constant(false));
}


namespace detail
{
    template<class Vertex>
    struct EmptyEdges;
}

}

#include "merge-tree.hpp"

#endif
