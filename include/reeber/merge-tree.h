#ifndef REEBER_MERGE_TREE_H
#define REEBER_MERGE_TREE_H

#include <vector>
#include <map>

#include <boost/foreach.hpp>
#include <boost/lambda/lambda.hpp>

namespace reeber
{

template<class Vertex_, class Value_>
struct MergeTreeNode
{
    typedef                     Vertex_                         Vertex;
    typedef                     Value_                          Value;

    typedef                     std::pair<Value, Vertex>        ValueVertex;

    typedef                     MergeTreeNode*                  Neighbor;
    typedef                     std::vector<Neighbor>           Neighbors;

    // TODO?: add operator<(const MergeTreeNode& other)

    Neighbor                    parent;
    Neighbors                   children;

    Vertex                      vertex;
    Value                       value;
    std::vector<ValueVertex>    vertices;       // degree-2 vertices on the path from this node to its parent

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

        template<class T>
        bool        cmp(const T& x, const T& y) const   { return negate_ ? x > y : x < y; }

        const VertexNeighborMap&
                    nodes() const                       { return nodes_; }

    public:
        struct      CollapseEvent   { static const char* name() { return "collapsed"; } };
        struct      EraseEvent      { static const char* name() { return "erased"; } };

    private:
        VertexNeighborMap&
                    nodes()                             { return nodes_; }

        static
        Neighbor&   compressed_parent(Neighbor n)       { return (Neighbor&) n->aux; }

        template<class MT, class T, class F, class C>
        friend void
        compute_merge_tree(MT& mt, const T& t, const F& f, const C& c);

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


// TODO: set up union topology from the trees and feed it into compute_merge_tree()
template<class MergeTree>
void    merge(std::vector<MergeTree>& trees);


}

#include "merge-tree.hpp"

#endif
