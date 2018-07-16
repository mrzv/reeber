#ifndef REEBER_TRIPLET_MERGE_TREE_H
#define REEBER_TRIPLET_MERGE_TREE_H

#include <vector>
#include <unordered_map>
#include <tuple>
#include <set>

#include "parallel-tbb.h"

#include "serialization.h"
#include "format.h"

namespace reeber
{
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

    bool                        operator< (const TripletMergeTreeNode& other) const     { return std::tie(value, vertex) <  std::tie(other.value, other.vertex); }
    bool                        operator<=(const TripletMergeTreeNode& other) const     { return std::tie(value, vertex) <= std::tie(other.value, other.vertex); }
    bool                        operator> (const TripletMergeTreeNode& other) const     { return std::tie(value, vertex) >  std::tie(other.value, other.vertex); }
    bool                        operator>=(const TripletMergeTreeNode& other) const     { return std::tie(value, vertex) >= std::tie(other.value, other.vertex); }

    bool                        operator==(const TripletMergeTreeNode& other) const     { return std::tie(vertex, value) == std::tie(other.vertex, other.value); }
    bool                        operator!=(const TripletMergeTreeNode& other) const     { return !(*this == other); }

    std::tuple<Neighbor, Neighbor>
                                parent() const                                          { Parent p = parent_; return std::make_tuple(p.through, p.to); }
    static Parent               make_parent(Neighbor s, Neighbor v)                     { return { s, v }; }

    Vertex                      vertex;
    Value                       value;
    atomic<Parent>              parent_;
    Neighbor                    cur_deepest;
    VerticesVector              vertices;

    friend std::ostream&        operator<<(std::ostream& os, const TripletMergeTreeNode& n) { os << "Node(vertex = " << n.vertex  << ", value = " << n.value << ")"; return os; }

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
                    ~TripletMergeTree()                 { if (is_node_owner_) for (auto n : nodes_) delete_node(n.second); }

        // It's Ok to move the tree; it's not Ok to copy it (because of the dynamically allocated nodes)
                    TripletMergeTree(const TripletMergeTree&)   =delete;
        TripletMergeTree&
                    operator=(const TripletMergeTree&)          =delete;
                    TripletMergeTree(TripletMergeTree&&)        =default;
        TripletMergeTree&
                    operator=(TripletMergeTree&&)               =default;

        std::tuple<Neighbor,Neighbor>
                    repair(const Neighbor u);
        Neighbor    add(const Vertex& x, Value v);
        Neighbor    find_or_add(const Vertex& x, Value v);
        Neighbor    add_or_update(const Vertex& x, Value v);
        void        make_shallow_copy(TripletMergeTree& other);
        void        make_deep_copy(TripletMergeTree& other);
        void        link(const Neighbor u, const Neighbor s, const Neighbor v)
                                                        { u->parent_ = Node::make_parent(s, v); }
        bool        cas_link(const Neighbor u, const Neighbor os, const Neighbor ov, const Neighbor s, const Neighbor v)
                                                        { auto op = Node::make_parent(os,ov); auto p = Node::make_parent(s,v); return compare_exchange(u->parent_, op, p); }

        void        merge(Neighbor u, Neighbor v);
        void        merge(Neighbor u, Neighbor s, Neighbor v);
        Neighbor    representative(Neighbor u, Neighbor a) const;

        Neighbor    find_deepest(const Neighbor u);

        size_t      size() const                        { return nodes_.size(); }

        bool        contains(const Vertex& x) const     { return nodes_.find(x) != nodes_.end(); }

        void        swap(TripletMergeTree& other)       { std::swap(negate_, other.negate_); nodes_.swap(other.nodes_); }

        bool        negate() const                      { return negate_; }
        void        set_negate(bool negate)             { negate_ = negate; }

        template<class T>
        bool        cmp(const T& x, const T& y) const   { return negate_ ? x > y : x < y; }
        bool        cmp(Neighbor x, Neighbor y) const   { return cmp(*x, *y); }

        Neighbor    operator[](const Vertex& x) const   { return nodes_.find(x)->second; }

        const VertexNeighborMap& nodes() const          { return nodes_; }

        friend struct ::reeber::Serialization<TripletMergeTree>;

        Neighbor    new_node()                          { Neighbor p = alloc_.allocate(1); alloc_.construct(p); return p; }
        void        delete_node(Neighbor p)             { alloc_.destroy(p); alloc_.deallocate(p,1); }

    private:
        VertexNeighborMap& nodes()                      { return nodes_; }

        template<class Vert, class Val, class T, class F>
        friend void
        compute_merge_tree(TripletMergeTree<Vert, Val>& mt, const T& t, const F& f);

        template<class Vert, class Val, class S>
        friend void
        remove_degree_two(TripletMergeTree<Vert, Val>& mt, const S& s);

        template<class Vertex, class Value>
        friend void
        repair(TripletMergeTree<Vertex, Value>& mt);

        template<class Vert, class Val, class T, class F>
        friend void
        compute_merge_tree2(TripletMergeTree<Vert, Val>& mt, const T& t, const F& f);

        template<class Vert, class Val, class F>
        friend void
        traverse_persistence(const TripletMergeTree<Vert, Val>& mt, const F& f);

        template<class Vert, class Val, class S>
        friend set<Vert>
        sparsify_keep(TripletMergeTree<Vert, Val>& mt, const S& s);

        template<class Vert, class Val, class S>
        friend void
        sparsify(TripletMergeTree<Vert, Val>& out, TripletMergeTree<Vert, Val>& in, const S& s);

        template<class Vert, class Val, class S>
        friend void
        sparsify(TripletMergeTree<Vert, Val>& mt, const S& s);

        template<class Vert, class Val>
        friend typename TripletMergeTree<Vert, Val>::Neighbor
        representative(TripletMergeTree<Vert, Val>& mt, typename TripletMergeTree<Vert, Val>::Neighbor u, typename TripletMergeTree<Vert, Val>::Neighbor a);

        template<class Vert, class Val, class E>
        friend void
        merge(TripletMergeTree<Vert, Val>& mt1, TripletMergeTree<Vert, Val>& mt2, const E& edges, bool ignore_missing_edges);

    private:
        bool                        negate_;
        VertexNeighborMap           nodes_;
        allocator<Node>             alloc_;
        bool                        is_node_owner_ { true };
};

/**
 * Topology defines a range vertices() and a link(v) function;
 *          vertices should be allowed to repeat (will simplify uniting multiple trees).
 */
template<class Vertex, class Value, class Topology, class Function>
void compute_merge_tree(TripletMergeTree<Vertex, Value>& mt, const Topology& topology, const Function& f);

template<class Vertex, class Value, class Special>
void remove_degree_two(TripletMergeTree<Vertex, Value>& mt, const Special& special);

template<class Vertex, class Value>
void repair(TripletMergeTree<Vertex, Value>& mt);

template<class Vertex, class Value, class Topology, class Function>
void compute_merge_tree2(TripletMergeTree<Vertex, Value>& mt, const Topology& topology, const Function& f);

template<class Vertex, class Value, class Functor>
void traverse_persistence(const TripletMergeTree<Vertex, Value>& mt, const Functor& f);

template<class Vertex, class Value, class Special>
void sparsify(TripletMergeTree<Vertex, Value>& out, TripletMergeTree<Vertex, Value>& in, const Special& special);

template<class Vertex, class Value, class Special>
void sparsify(TripletMergeTree<Vertex, Value>& mt, const Special& special);

template<class Vertex, class Value>
typename TripletMergeTree<Vertex, Value>::Neighbor
representative(TripletMergeTree<Vertex, Value>& mt, typename TripletMergeTree<Vertex, Value>::Neighbor u, typename TripletMergeTree<Vertex, Value>::Neighbor a);

template<class Vertex, class Value, class Edges>
void merge(TripletMergeTree<Vertex, Value>& mt1, TripletMergeTree<Vertex, Value>& mt2, const Edges& edges, bool ignore_missing_edges = false);

template<class Vertex, class Value, class Special>
set<Vertex>
sparsify_keep(TripletMergeTree<Vertex, Value>& mt, const Special& special);

}

#include "triplet-merge-tree.hpp"

#endif
