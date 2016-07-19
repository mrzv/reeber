#ifndef REEBER_TRIPLET_MERGE_TREE_H
#define REEBER_TRIPLET_MERGE_TREE_H

#include <vector>
#include <unordered_map>
#include <tuple>

#include <boost/range/adaptor/map.hpp>
#include <boost/range/algorithm.hpp>

#include "parallel-tbb.h"

#include "serialization.h"
#include "format.h"

namespace reeber
{

namespace ba = boost::adaptors;

template<class Vertex_, class Value_>
struct TripletMergeTreeNode
{
    typedef                     Vertex_                         Vertex;
    typedef                     Value_                          Value;

    typedef                     std::pair<Value, Vertex>        ValueVertex;
    typedef                     std::vector<ValueVertex>        VerticesVector;

    typedef                     TripletMergeTreeNode*           Neighbor;
    struct Parent
    {
        Neighbor    through;
        Neighbor    to;
    };

    bool                        operator<(const TripletMergeTreeNode& other) const     { return value < other.value || (value == other.value && vertex < other.vertex); }
    bool                        operator>(const TripletMergeTreeNode& other) const     { return value > other.value || (value == other.value && vertex > other.vertex); }

    std::tuple<Neighbor, Neighbor>
                                parent() const                                          { Parent p = parent_; return {p.through, p.to}; }
    static Parent               make_parent(Neighbor s, Neighbor v)                     { return { s, v }; }

    Vertex                      vertex;
    Value                       value;
    atomic<Parent>              parent_;
    Neighbor                    cur_deepest;
    VerticesVector              vertices;
};

template<class Vertex_, class Value_>
class TripletMergeTree
{
    public:
        typedef     Vertex_                             Vertex;
        typedef     Value_                              Value;

        typedef     TripletMergeTreeNode<Vertex,Value>  Node;
        typedef     typename Node::Neighbor             Neighbor;

        typedef     map<Vertex, Neighbor>               VertexNeighborMap;

    public:
                    TripletMergeTree(bool negate = false):
                        negate_(negate)                 {}
                    ~TripletMergeTree()                 { for (auto n : nodes_) delete n.second; }

        // It's Ok to move the tree; it's not Ok to copy it (because of the dynamically allocated nodes)
                    TripletMergeTree(const TripletMergeTree&)   =delete;
        TripletMergeTree&
                    operator=(const TripletMergeTree&)          =delete;
                    TripletMergeTree(TripletMergeTree&&)        =default;
        TripletMergeTree&
                    operator=(TripletMergeTree&&)               =default;

        void        repair(const Neighbor u);
        Neighbor    add(const Vertex& x, Value v);
        Neighbor    find_or_add(const Vertex& x, Value v);
        Neighbor    add_or_update(const Vertex& x, Value v);
        void        link(const Neighbor u, const Neighbor s, const Neighbor v)
                                                        { u->parent_ = Node::make_parent(s, v); }
        bool        cas_link(const Neighbor u, const Neighbor os, const Neighbor ov, const Neighbor s, const Neighbor v)
                                                        { auto op = Node::make_parent(os,ov); auto p = Node::make_parent(s,v); return compare_exchange(u->parent_, op, p); }

        Neighbor    find_deepest(const Neighbor u);

        size_t      size() const                        { return nodes_.size(); }

        bool        contains(const Vertex& x) const     { return nodes_.find(x) != nodes_.end(); }

        void        swap(TripletMergeTree& other)       { std::swap(negate_, other.negate_); nodes_.swap(other.nodes_); }

        bool        negate() const                      { return negate_; }
        void        set_negate(bool negate)             { negate_ = negate; }

        template<class T>
        bool        cmp(const T& x, const T& y) const   { return negate_ ? x > y : x < y; }
        bool        cmp(Neighbor x, Neighbor y) const   { return cmp(*x, *y); }

        Neighbor    node(const Vertex& x) const         { return nodes_.find(x)->second; }
        Neighbor    operator[](const Vertex& x) const   { return nodes_.find(x)->second; }

        const VertexNeighborMap& nodes() const          { return nodes_; }

        friend struct ::reeber::Serialization<TripletMergeTree>;

    private:
        VertexNeighborMap& nodes()                      { return nodes_; }

        template<class Vert, class Val, class T, class F>
        friend void
        compute_merge_tree(TripletMergeTree<Vert, Val>& mt, const T& t, const F& f);

        template<class Vert, class Val, class S>
        friend void
        remove_degree_two(TripletMergeTree<Vert, Val>& mt, const S& s);

        template<class Vert, class Val, class T, class F>
        friend void
        compute_merge_tree2(TripletMergeTree<Vert, Val>& mt, const T& t, const F& f);

        template<class Vert, class Val, class F>
        friend void
        traverse_persistence(const TripletMergeTree<Vert, Val>& mt, const F& f);

        template<class Vert, class Val, class S>
        friend void
        sparsify(TripletMergeTree<Vert, Val>& out, TripletMergeTree<Vert, Val>& in, const S& s);

        template<class Vert, class Val, class S>
        friend void
        sparsify(TripletMergeTree<Vert, Val>& mt, const S& s);

        template<class Vert, class Val>
        friend void
        merge(TripletMergeTree<Vert, Val>& mt, typename TripletMergeTree<Vert, Val>::Neighbor u, typename TripletMergeTree<Vert, Val>::Neighbor s, typename TripletMergeTree<Vert, Val>::Neighbor v);

        template<class Vert, class Val, class E>
        friend void
        merge(TripletMergeTree<Vert, Val>& mt1, TripletMergeTree<Vert, Val>& mt2, const E& edges);

    private:
        bool                        negate_;
        VertexNeighborMap           nodes_;
};

/**
 * Topology defines a range vertices() and a link(v) function;
 *          vertices should be allowed to repeat (will simplify uniting multiple trees).
 */
template<class Vertex, class Value, class Topology, class Function>
void compute_merge_tree(TripletMergeTree<Vertex, Value>& mt, const Topology& topology, const Function& f);

template<class Vertex, class Value, class Special>
void remove_degree_two(TripletMergeTree<Vertex, Value>& mt, const Special& special);

template<class Vertex, class Value, class Topology, class Function>
void compute_merge_tree2(TripletMergeTree<Vertex, Value>& mt, const Topology& topology, const Function& f);

template<class Vertex, class Value, class Functor>
void traverse_persistence(const TripletMergeTree<Vertex, Value>& mt, const Functor& f);

template<class Vertex, class Value, class Special>
void sparsify(TripletMergeTree<Vertex, Value>& out, TripletMergeTree<Vertex, Value>& in, const Special& special);

template<class Vertex, class Value, class Special>
void sparsify(TripletMergeTree<Vertex, Value>& mt, const Special& special);

template<class Vertex, class Value>
void merge(TripletMergeTree<Vertex, Value>& mt, typename TripletMergeTree<Vertex, Value>::Neighbor u, typename TripletMergeTree<Vertex, Value>::Neighbor s, typename TripletMergeTree<Vertex, Value>::Neighbor v);

template<class Vertex, class Value, class Edges>
void merge(TripletMergeTree<Vertex, Value>& mt1, TripletMergeTree<Vertex, Value>& mt2, const Edges& edges);

}

#include "triplet-merge-tree.hpp"

#endif
