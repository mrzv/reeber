#include <stack>

#include <dlog/stats.h>
#include <dlog/log.h>
#include <dlog/counters.h>

#include "tree-union.h"

template<class Vertex, class Value>
typename reeber::MergeTree<Vertex, Value>::Neighbor
reeber::MergeTree<Vertex, Value>::
add(const Vertex& x, Value v)
{
    Neighbor    n = new Node;
    nodes_[x] = n;

    n->vertex = x;
    n->value  = v;
    n->parent = 0;
    aux_neighbor(n) = 0;

    return n;
}

template<class Vertex, class Value>
typename reeber::MergeTree<Vertex, Value>::Neighbor
reeber::MergeTree<Vertex, Value>::
find(const Vertex& x) const
{
    // aux_neighbor functions as compressed parent

    Neighbor xn = (*this)[x];
    Neighbor res = xn;
    while (aux_neighbor(res) != 0)
    {
        COUNTER(FindStepEvent)++;
        res = aux_neighbor(res);
    }

    // compress the path to res
    Neighbor up = aux_neighbor(xn);
    while (up != 0)
    {
        aux_neighbor(xn) = res;
        xn = up;
        up = aux_neighbor(xn);
    }

    return res;
}

template<class Vertex, class Value>
void
reeber::MergeTree<Vertex, Value>::
link(Neighbor xn, Neighbor yn)
{
    yn->parent = xn;
    aux_neighbor(yn) = xn;
    xn->children.push_back(yn);
}


/* Compute merge tree */
template<class MergeTree, class Topology, class Function, class Collapsible>
void
reeber::compute_merge_tree(MergeTree& mt, const Topology& topology, const Function& f, const Collapsible& collapsible)
{
    dlog::prof << "compute-merge-tree";
    typedef     typename Topology::Vertex       Vertex;
    typedef     typename Topology::Link         Link;
    typedef     typename Function::Value        Value;
    typedef     std::pair<Value, Vertex>        ValueVertexPair;
    typedef     typename MergeTree::Neighbor    Neighbor;

    std::vector<ValueVertexPair>     vertices;
    vertices.reserve(topology.size());
    BOOST_FOREACH(const Vertex& v, topology.vertices())
        vertices.push_back(std::make_pair(f(v), v));

    if (mt.negate())
        std::sort(vertices.begin(), vertices.end(), std::greater<ValueVertexPair>());
    else
        std::sort(vertices.begin(), vertices.end(), std::less<ValueVertexPair>());


    BOOST_FOREACH(const ValueVertexPair& fu, vertices)
    {
        Value val; Vertex u;
        boost::tie(val, u) = fu;

        std::set<Neighbor>  roots;
        BOOST_FOREACH(const Vertex& v, topology.link(u))
        {
            if (mt.contains(v))
            {
                Neighbor v_root = mt.find(v);
                roots.insert(v_root);
            }
        }

        if (roots.size() == 1 && collapsible(u))
        {
            Neighbor n = *roots.begin();
            n->vertices.push_back(fu);
            mt.nodes()[u] = n;
            COUNTER(typename MergeTree::CollapseEvent)++;
        } else
        {
            Neighbor u_root = mt.add(u, val);
            BOOST_FOREACH(Neighbor n, roots)
                mt.link(u_root, n);
        }
    }

    // clean up
    typename MergeTree::VertexNeighborMap::iterator it = mt.nodes().begin();
    while(it != mt.nodes().end())
    {
        if (it->first != it->second->vertex)
        {
            COUNTER(typename MergeTree::EraseEvent)++;
            mt.nodes().erase(it++);
        }
        else
        {
            MergeTree::aux_neighbor(it->second) = 0;      // reset aux
            ++it;
        }
    }
    dlog::prof >> "compute-merge-tree";
}

template<class MergeTree, class Functor>
void
reeber::traverse_persistence(const MergeTree& mt, const Functor& f)
{
    dlog::prof << "traverse-persistence";

    typedef     typename MergeTree::Neighbor        Neighbor;

    // find root
    Neighbor root = mt.find_root();

    std::stack<Neighbor> s;
    s.push(root);
    while(!s.empty())
    {
        Neighbor n = s.top();
        if (n->children.empty())
        {
            MergeTree::aux_neighbor(n) = n;
            s.pop();
        } else if (!MergeTree::aux_neighbor(n->children[0]))
        {
            BOOST_FOREACH(Neighbor child, n->children)
                s.push(child);
        } else
        {
            // find the deepest subtree
            Neighbor deepest = MergeTree::aux_neighbor(n->children[0]);
            for (unsigned i = 1; i < n->children.size(); ++i)
            {
                Neighbor child = n->children[i];
                if (mt.cmp(*MergeTree::aux_neighbor(child), *MergeTree::aux_neighbor(deepest)))
                    deepest = child;
            }

            MergeTree::aux_neighbor(n) = MergeTree::aux_neighbor(deepest);

            // report the rest of the pairs
            BOOST_FOREACH(Neighbor child, n->children)
            {
                if (child == deepest)
                    continue;
                f(MergeTree::aux_neighbor(child), n, MergeTree::aux_neighbor(deepest));
            }
            s.pop();
        }
    }

    f(MergeTree::aux_neighbor(root), root, MergeTree::aux_neighbor(root));

    mt.reset_aux();

    dlog::prof >> "traverse-persistence";
}

template<class MergeTree, class Special>
void
reeber::sparsify(MergeTree& mt, const Special& special)
{
    // Traverses the tree in a post-order fashion;
    // aux_neighbor for internal nodes keeps track of the deepest leaf
    // aux_neighbor for the leaves marks whether the subtree needs to be preserved
    // (this information may change in the course of the traversal)

    dlog::prof << "sparsify";

    typedef     typename MergeTree::Neighbor        Neighbor;

    Neighbor root = mt.find_root();

    std::stack<Neighbor> s;
    s.push(root);
    while(!s.empty())
    {
        Neighbor n = s.top();
        if (n->children.empty())
        {
            MergeTree::aux_neighbor(n) = 0;
            if (special(n->vertex))
                MergeTree::aux_neighbor(n) = n;
            AssertMsg(n->parent, "Parent of " << n->vertex << " must be present");
            if (mt.cmp(*n, *MergeTree::aux_neighbor(n->parent)))
                MergeTree::aux_neighbor(n->parent) = n;
            s.pop();
        } else if (!MergeTree::aux_neighbor(n))
        {
            BOOST_FOREACH(Neighbor child, n->children)
                s.push(child);
            MergeTree::aux_neighbor(n) = n;
        } else
        {
            Neighbor deepest = MergeTree::aux_neighbor(n);

            unsigned end = n->children.size();
            for (unsigned i = 0; i < end; )
            {
                Neighbor child = n->children[i];
                if (child == deepest)
                {
                    ++i; continue;
                }
                if (MergeTree::aux_neighbor(child))     // needs to be preserved
                {
                    MergeTree::aux_neighbor(deepest) = deepest;     // mark this subtree as in need of preserving
                    ++i; continue;
                }

                std::swap(n->children[i], n->children[end-1]);
                --end;
            }

            // remove subtrees at n->children[end..]
            std::stack<Neighbor> rms;
            for (unsigned i = end; i < n->children.size(); ++i)
                rms.push(n->children[i]);
            n->children.resize(end);        // this won't actually shrink the space, but it probably doesn't matter

            while (!rms.empty())
            {
                Neighbor rm = rms.top();
                rms.pop();
                BOOST_FOREACH(Neighbor child, rm->children)
                    rms.push(child);
                mt.nodes().erase(rm->vertex);
                delete rm;
            }

            s.pop();
        }
    }

    mt.reset_aux();

    dlog::prof >> "sparsify";
}

#include "merge-tree-serialization.h"       // TODO: remove once the code below is fixed
template<class MergeTree, class Special>
void
reeber::sparsify(MergeTree& out, const MergeTree& in, const Special& special)
{
    // TODO: a hideous hack for now; implement properly later
    diy::BinaryBuffer bb;
    LOG_SEV(debug) << "Saving the tree: " << in.size();
    diy::save(bb, in);
    bb.reset();
    diy::load(bb, out);
    LOG_SEV(debug) << "Loaded the tree: " << out.size();
    sparsify(out, special);
}

template<class MergeTree>
void
reeber::merge(MergeTree& mt, const std::vector<MergeTree>& trees)
{
    dlog::prof << "merge";
    TreeUnionTopology<MergeTree> union_topology(trees);
    TreeUnionFunction<MergeTree> union_function(trees);
    compute_merge_tree(mt, union_topology, union_function, boost::lambda::constant(false));
    dlog::prof >> "merge";
}

template<class MergeTree, class Preserve>
void
reeber::remove_degree2(MergeTree& mt, const Preserve& preserve)
{
    dlog::prof << "remove-degree2";
    typedef     typename MergeTree::Neighbor        Neighbor;

    Neighbor root = mt.find_root();

    std::stack<Neighbor> s;
    s.push(root);
    while(!s.empty())
    {
        Neighbor n = s.top();
        s.pop();
        for (unsigned i = 0; i < n->children.size(); ++i)
        {
            Neighbor child = n->children[i];
            if (child->children.size() == 1)
            {
                Neighbor descendant = child->children[0];
                while (descendant->children.size() == 1)
                    descendant = descendant->children[0];

                // save the nodes that need to be preserved
                Neighbor cur = descendant->parent;
                while (cur != n)
                {
                    if (preserve(cur->vertex))
                    {
                        typedef     typename MergeTree::Node::ValueVertex       ValueVertex;
                        descendant->vertices.push_back(ValueVertex(cur->value, cur->vertex));
                        BOOST_FOREACH(const ValueVertex& vv, cur->vertices)
                            descendant->vertices.push_back(vv);
                    }
                }

                // remove the path
                cur = descendant->parent;
                while (cur != n)
                {
                    Neighbor rm = cur;
                    cur = cur->parent;
                    mt.nodes().erase(rm->vertex);
                    delete rm;
                }

                // link descendant and n
                descendant->parent = n;
                n->children[i] = descendant;
                s.push(descendant);
            } else
                s.push(child);
        }
    }
    dlog::prof >> "remove-degree2";
}
