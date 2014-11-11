#include <stack>

#include <boost/range/adaptor/map.hpp>

#include <dlog/stats.h>
#include <dlog/log.h>
#include <dlog/counters.h>

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
    typedef     typename MergeTree::Neighbor        Neighbor;
    std::stack<Neighbor> s;
    namespace ba = boost::adaptors;

    // find root
    Neighbor root;
    BOOST_FOREACH(Neighbor n, mt.nodes() | ba::map_values)
        if (!n->parent)
        {
            root = n;
            break;
        }

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

    // reset aux
    BOOST_FOREACH(Neighbor n, mt.nodes() | ba::map_values)
        MergeTree::aux_neighbor(n) = 0;
}

