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
    compressed_parent(n) = 0;

    return n;
}

template<class Vertex, class Value>
typename reeber::MergeTree<Vertex, Value>::Neighbor
reeber::MergeTree<Vertex, Value>::
find(const Vertex& x) const
{
    Neighbor xn = (*this)[x];
    Neighbor res = xn;
    while (compressed_parent(res) != 0)
    {
        COUNTER(FindStepEvent)++;
        res = compressed_parent(res);
    }

    // compress the path to res
    Neighbor up = compressed_parent(xn);
    while (up != 0)
    {
        compressed_parent(xn) = res;
        xn = up;
        up = compressed_parent(xn);
    }

    return res;
}

template<class Vertex, class Value>
void
reeber::MergeTree<Vertex, Value>::
link(Neighbor xn, Neighbor yn)
{
    yn->parent = xn;
    compressed_parent(yn) = xn;
    xn->children.push_back(yn);
}


/* Compute merge tree */
template<class MergeTree, class Topology, class Function, class Collapsible>
void
reeber::compute_merge_tree(MergeTree& mt, const Topology& topology, const Function& f, const Collapsible& collapsible)
{
    PROF << "compute-merge-tree";
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
            MergeTree::compressed_parent(it->second) = 0;      // reset aux
            ++it;
        }
    }
    PROF >> "compute-merge-tree";
}
