#ifndef REEBER_PATH_MERGE_TREE_H
#define REEBER_PATH_MERGE_TREE_H

#include "parallel-tbb.h"
#include "serialization.h"

/*
   Implements a simple merge tree that only keeps track of the parent at each
   node and supports construction only via merging paths to the root
*/

namespace reeber
{

template<class Vertex_, class Value_>
struct PathMergeTreeNode
{
    typedef                     Vertex_                         Vertex;
    typedef                     Value_                          Value;

    typedef                     PathMergeTreeNode*              Neighbor;

    bool                        operator<(const PathMergeTreeNode& other) const     { return value < other.value || (value == other.value && vertex < other.vertex); }
    bool                        operator>(const PathMergeTreeNode& other) const     { return value > other.value || (value == other.value && vertex > other.vertex); }

    Vertex                      vertex;
    Value                       value;
    atomic<Neighbor>            parent;
};

template<class Vertex_, class Value_>
class PathMergeTree
{
    public:
        typedef     Vertex_                             Vertex;
        typedef     Value_                              Value;

        typedef     PathMergeTreeNode<Vertex,Value>     Node;
        typedef     typename Node::Neighbor             Neighbor;

        typedef     map<Vertex, Neighbor>               VertexNeighborMap;

    public:
                    PathMergeTree(bool negate = false):
                        negate_(negate)                 {}
                    ~PathMergeTree()                    { for (auto n : nodes_) delete_node(n.second); }

        // It's Ok to move the tree; it's not Ok to copy it (because of the dynamically allocated nodes)
                    PathMergeTree(const PathMergeTree&)     =delete;
        PathMergeTree&
                    operator=(const PathMergeTree&)         =delete;
                    PathMergeTree(PathMergeTree&&)          =default;
        PathMergeTree&
                    operator=(PathMergeTree&&)              =default;

        Neighbor    add(const Vertex& x, Value v);
        void        merge(Neighbor u, Neighbor v);

        Neighbor    add_or_update(const Vertex& x, Value v);
        Neighbor    find_or_add(const Vertex& x, Value v);

        size_t      size() const                        { return nodes_.size(); }
        bool        contains(const Vertex& x) const     { return nodes_.find(x) != nodes_.end(); }
        void        swap(PathMergeTree& other)          { std::swap(negate_, other.negate_); nodes_.swap(other.nodes_); }

        bool        negate() const                      { return negate_; }
        void        set_negate(bool negate)             { negate_ = negate; }

        template<class T>
        bool        cmp(const T& x, const T& y) const   { return negate_ ? x > y : x < y; }
        bool        cmp(Neighbor x, Neighbor y) const   { return cmp(*x, *y); }

        Neighbor    operator[](const Vertex& x) const   { return nodes_.find(x)->second; }

        const VertexNeighborMap& nodes() const          { return nodes_; }

        friend struct ::reeber::Serialization<PathMergeTree>;

        Neighbor    new_node()                          { Neighbor p = alloc_.allocate(1); alloc_.construct(p); return p; }
        void        delete_node(Neighbor p)             { alloc_.destroy(p); alloc_.deallocate(p,1); }

    private:
        VertexNeighborMap& nodes()                      { return nodes_; }

        template<class Vert, class Val, class T, class F>
        friend void
        compute_merge_tree2(PathMergeTree<Vert, Val>& mt, const T& t, const F& f);

        template<class Vert, class Val, class E>
        friend void
        merge(PathMergeTree<Vert, Val>& mt1, PathMergeTree<Vert, Val>& mt2, const E& edges);

    private:
        bool                        negate_;
        VertexNeighborMap           nodes_;
        allocator<Node>             alloc_;
};

/**
 * Topology defines a range vertices() and a link(v) function;
 *          vertices should be allowed to repeat (will simplify uniting multiple trees).
 */
template<class Vertex, class Value, class Topology, class Function>
void compute_merge_tree2(PathMergeTree<Vertex, Value>& mt, const Topology& topology, const Function& f);

template<class Vertex, class Value, class Edges>
void merge(PathMergeTree<Vertex, Value>& mt1, PathMergeTree<Vertex, Value>& mt2, const Edges& edges);

}

#include "path-merge-tree.hpp"

#endif
