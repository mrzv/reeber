#include <stack>

#include <dlog/stats.h>
#include <dlog/log.h>
#include <dlog/counters.h>

#include <boost/bind.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/range/algorithm_ext/push_back.hpp>

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
find(Neighbor xn) const
{
    // aux_neighbor functions as compressed parent

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
    typedef     typename Function::Value        Value;
    typedef     std::pair<Value, Vertex>        ValueVertexPair;
    typedef     typename MergeTree::Neighbor    Neighbor;

    std::vector<ValueVertexPair>     vertices;
    vertices.reserve(topology.size());
    BOOST_FOREACH(Vertex v, topology.vertices())
        vertices.push_back(std::make_pair(f(v), v));

    if (mt.negate())
        std::sort(vertices.begin(), vertices.end(), std::greater<ValueVertexPair>());
    else
        std::sort(vertices.begin(), vertices.end(), std::less<ValueVertexPair>());

    LOG_SEV(debug) << "Computing merge tree out of " << vertices.size() << " vertices";


    BOOST_FOREACH(const ValueVertexPair& fu, vertices)
    {
        Value val; Vertex u;
        boost::tie(val, u) = fu;

        std::set<Neighbor>  roots;
        BOOST_FOREACH(Vertex v, topology.link(u))
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
    Neighbor root = 0;
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
            if (!it->second->parent) root = it->second;
            ++it;
        }
    }

    // pull out the correct root
    if (!root->vertices.empty())
    {
        typedef     typename MergeTree::Node::VerticesVector        VerticesVector;
        typename VerticesVector::iterator max_elem = std::max_element(root->vertices.begin(), root->vertices.end(),
                                                                      boost::bind(&MergeTree::template cmp<typename VerticesVector::value_type>,
                                                                                  &mt, _1, _2));
        Neighbor new_root = mt.add(max_elem->second, max_elem->first);
        std::swap(*max_elem, root->vertices.back());
        root->vertices.pop_back();
        root->parent = new_root;
        new_root->children.push_back(root);
        LOG_SEV(debug) << "Pulled out new root: " << new_root->vertex << " " << new_root->value;
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
            Neighbor deepest = n->children[0];
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
        if (n->children.empty())                // leaf
        {
            MergeTree::aux_neighbor(n) = 0;
            if (n->any_vertex(special))
                MergeTree::aux_neighbor(n) = n;
            AssertMsg(n->parent, "Parent of " << n->vertex << " must be present");
            if (mt.cmp(*n, *MergeTree::aux_neighbor(n->parent)))
                MergeTree::aux_neighbor(n->parent) = n;
            s.pop();
        } else if (!MergeTree::aux_neighbor(n)) // on the way down
        {
            BOOST_FOREACH(Neighbor child, n->children)
                s.push(child);
            MergeTree::aux_neighbor(n) = n;
        } else                                  // on the way up
        {
            Neighbor deepest = MergeTree::aux_neighbor(n);
            AssertMsg(mt.contains(deepest->vertex), "deepest must be in the tree");
            AssertMsg(deepest->children.empty(),    "deepest must be a leaf");

            // propagate deepest up
            if (n->parent && mt.cmp(*deepest, *MergeTree::aux_neighbor(n->parent)))
                MergeTree::aux_neighbor(n->parent) = deepest;

            bool preserve = n->any_vertex(special) || !n->parent;
            unsigned end = n->children.size();
            for (unsigned i = 0; i < end; )
            {
                Neighbor child         = n->children[i];
                Neighbor child_deepest = MergeTree::aux_neighbor(child);
                if (!child_deepest) child_deepest = child;
                if (child_deepest && MergeTree::aux_neighbor(child_deepest))     // needs to be preserved
                {
                    preserve = true;
                    ++i; continue;
                }
                if (child_deepest == deepest)
                {
                    ++i; continue;
                }

                --end;
                std::swap(n->children[i], n->children[end]);
            }

            if (!preserve)
                MergeTree::aux_neighbor(deepest) = 0;         // this subtree can be removed
            else
                MergeTree::aux_neighbor(deepest) = deepest;

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

template<class MergeTree, class Special>
void
reeber::sparsify(MergeTree& out, const MergeTree& in, const Special& special)
{
    // Puts the sparsified tree into out
    // Traverses the tree in a post-order fashion;
    // aux_neighbor for internal nodes keeps track of the deepest leaf
    // aux_neighbor for the leaves marks whether the subtree needs to be preserved
    // (this information may change in the course of the traversal)

    dlog::prof << "sparsify";

    typedef     typename MergeTree::Neighbor        Neighbor;
    typedef     typename MergeTree::Vertex          Vertex;
    typedef     std::map<Vertex, Neighbor>          VertexNeighborMap;

    VertexNeighborMap   deepest_root;       // map from deepest node to its current root in the new tree

    Neighbor root = in.find_root();
    std::stack<Neighbor> s;
    s.push(root);
    while(!s.empty())
    {
        Neighbor n = s.top();
        if (n->children.empty())                // leaf
        {
            MergeTree::aux_neighbor(n) = 0;
            if (n->any_vertex(special))
            {
                MergeTree::aux_neighbor(n) = n;
                Neighbor new_n = out.add(n->vertex, n->value);
                deepest_root[new_n->vertex] = new_n;
            }
            AssertMsg(n->parent, "Parent of " << n->vertex << " must be present");
            if (in.cmp(*n, *MergeTree::aux_neighbor(n->parent)))
                MergeTree::aux_neighbor(n->parent) = n;
            s.pop();
        } else if (!MergeTree::aux_neighbor(n)) // on the way down
        {
            BOOST_FOREACH(Neighbor child, n->children)
                s.push(child);
            MergeTree::aux_neighbor(n) = n;
        } else                                  // on the way up
        {
            Neighbor deepest = MergeTree::aux_neighbor(n);
            AssertMsg(in.contains(deepest->vertex), "deepest must be in the tree");
            AssertMsg(deepest->children.empty(),    "deepest must be a leaf");

            // propagate deepest up
            if (n->parent && in.cmp(*deepest, *MergeTree::aux_neighbor(n->parent)))
                MergeTree::aux_neighbor(n->parent) = deepest;

            bool preserve = n->any_vertex(special) || !n->parent;
            unsigned end = n->children.size();
            for (unsigned i = 0; i < end; )
            {
                Neighbor child         = n->children[i];
                Neighbor child_deepest = MergeTree::aux_neighbor(child);
                if (!child_deepest) child_deepest = child;
                if (child_deepest && MergeTree::aux_neighbor(child_deepest))     // needs to be preserved
                    ++i;
                else if (child_deepest == deepest)
                    ++i;
                else
                {
                    --end;
                    std::swap(n->children[i], n->children[end]);
                }
            }

            if (preserve || end > 1)
            {
                Neighbor new_n = out.add(n->vertex, n->value);
                for (unsigned i = 0; i < end; ++i)
                {
                    Neighbor        child_deepest = MergeTree::aux_neighbor(n->children[i]);
                    if (!child_deepest)
                        child_deepest = n->children[i];
                    const Vertex&   v             = child_deepest->vertex;
                    typename VertexNeighborMap::iterator it = deepest_root.find(v);
                    if (it != deepest_root.end())
                    {
                        Neighbor y = it->second;
                        MergeTree::link(new_n, y);
                        it->second = new_n;
                    } else
                    {
                        Neighbor x = out.add(child_deepest->vertex, child_deepest->value);
                        MergeTree::link(new_n, x);
                        deepest_root[v] = new_n;
                    }
                }
                MergeTree::aux_neighbor(deepest) = deepest;
            }

            s.pop();
        }
    }

    in.reset_aux();
    out.reset_aux();

    dlog::prof >> "sparsify";
}

template<class MergeTree>
void
reeber::merge(MergeTree& mt, const std::vector<MergeTree>& trees)
{
    merge(mt, trees, detail::EmptyEdges<typename MergeTree::Vertex>());
}

template<class MergeTree, class Edges>
void
reeber::merge(MergeTree& mt, const std::vector<MergeTree>& trees, const Edges& edges)
{
    dlog::prof << "merge";

    // Fill and sort the nodes
    typedef     typename MergeTree::Neighbor        Neighbor;
    typedef     typename MergeTree::Vertex          Vertex;

    std::vector<Neighbor> nodes;
    for (unsigned i = 0; i < trees.size(); ++i)
        boost::push_back(nodes, trees[i].nodes() | ba::map_values);
    boost::sort(nodes,
                boost::bind(static_cast< bool (MergeTree::*) (Neighbor,Neighbor) const >(&MergeTree::cmp), mt, _1, _2));

    BOOST_FOREACH(Neighbor n, nodes)
    {
        Neighbor nn;
        if (!mt.contains(n->vertex))
        {
            nn = mt.add(n->vertex, n->value);

            // deal with the edges (only once per vertex)
            BOOST_FOREACH(Vertex v, edges(n->vertex))
                if (mt.contains(v))
                {
                    Neighbor cn = mt.find(v);
                    if (cn != nn)
                        mt.link(nn, cn);
                }
        }
        else
            nn = mt[n->vertex];

        for(size_t i = 0; i < n->children.size(); ++i)
        {
            assert(mt.contains(n->children[i]->vertex));
            Neighbor cn = mt.find(n->children[i]->vertex);
            if (cn != nn)
                mt.link(nn, cn);
        }
    }

    // reset aux
    BOOST_FOREACH(Neighbor n, mt.nodes() | ba::map_values)
        mt.aux_neighbor(n) = 0;

    dlog::prof >> "merge";
}

template<class Vertex_>
struct reeber::detail::
EmptyEdges
{
    typedef         Vertex_                                 Vertex;
    typedef         std::pair<Vertex*, Vertex*>             VertexRange;
    VertexRange     operator()(const Vertex& u) const       { return VertexRange(0,0); }
};


template<class MergeTree, class Preserve, class Special>
void
reeber::remove_degree2(MergeTree& mt, const Preserve& preserve, const Special& special)
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
            if (child->children.size() == 1 && !special(child->vertex))
            {
                Neighbor descendant = child->children[0];
                while (descendant->children.size() == 1 && !special(descendant->vertex))
                    descendant = descendant->children[0];

                // save the nodes that need to be preserved
                Neighbor cur = descendant->parent;
                while (cur != n)
                {
                    typedef     typename MergeTree::Node::ValueVertex       ValueVertex;
                    if (preserve(cur->vertex))
                        descendant->vertices.push_back(ValueVertex(cur->value, cur->vertex));

                    BOOST_FOREACH(const ValueVertex& vv, cur->vertices)
                        if (preserve(vv.second))
                            descendant->vertices.push_back(vv);
                    cur = cur->parent;
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

template<class MergeTree>
void
reeber::redistribute_vertices(MergeTree& mt)
{
    dlog::prof << "redistribute-vertices";
    typedef     typename MergeTree::Neighbor            Neighbor;
    typedef     typename MergeTree::Node::ValueVertex   ValueVertex;

    Neighbor root = mt.find_root();

    std::stack< std::pair<Neighbor,int> > s;
    s.push(std::make_pair(root,0));
    while(!s.empty())
    {
        Neighbor n = s.top().first;
        int&     c = s.top().second;

        if (c == n->children.size())
        {
            s.pop();

            // distribute the vertices
            std::set<ValueVertex>                       vertices(n->vertices.begin(), n->vertices.end());
            typename std::set<ValueVertex>::iterator    it = vertices.begin();
            n->vertices.clear();
            while(it != vertices.end())
            {
                // move vertex up
                ValueVertex parent_vv(n->parent->value, n->parent->vertex);
                if (n->parent && !mt.cmp(*it, parent_vv))
                {
                    if (mt.cmp(parent_vv, *it))
                        n->parent->vertices.push_back(*it);
                } else
                    n->vertices.push_back(*it);
                ++it;
            }
        } else
        {
            c += 1;
            s.push(std::make_pair(n->children[c-1], 0));
        }
    }
    dlog::prof >> "redistribute-vertices";
}
