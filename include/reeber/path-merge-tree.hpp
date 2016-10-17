#include <dlog/log.h>
#include <dlog/stats.h>

template<class Vertex, class Value>
typename reeber::PathMergeTree<Vertex, Value>::Neighbor
reeber::PathMergeTree<Vertex, Value>::
add(const Vertex& x, Value v)
{
    Neighbor n = new_node();
    n->vertex = x;
    n->value = v;
    n->parent = n;
    nodes_.emplace(x,n);
    return n;
}

template<class Vertex, class Value>
void
reeber::PathMergeTree<Vertex, Value>::
merge(Neighbor u, Neighbor v)
{
    while (u != v)
    {
        if (cmp(v,u))
            std::swap(u,v);

        Neighbor up = u->parent;
        if (u != up && !cmp(v, up))     // this handles special case when v == up
            u = up;
        else if (compare_exchange(u->parent, up, v))
        {
#ifndef REEBER_USE_TBB
            auto it = std::find(up->children.begin(), up->children.end(), u);
            if (it != up->children.end())
                up->children.erase(it);
            v->children.push_back(u);
#endif
            u = up;
        }
    }
}

template<class Vertex, class Value>
typename reeber::PathMergeTree<Vertex, Value>::Neighbor
reeber::PathMergeTree<Vertex, Value>::
add_or_update(const Vertex& x, Value v)
{
    auto it = nodes().find(x);
    if (it != nodes().end())
    {
        it->second->value = v;
        return it->second;
    }
    else
        return add(x, v);
}

template<class Vertex, class Value>
typename reeber::PathMergeTree<Vertex, Value>::Neighbor
reeber::PathMergeTree<Vertex, Value>::
find_or_add(const Vertex& x, Value v)
{
    auto it = nodes().find(x);
    if (it != nodes().end())
        return it->second;
    else
        return add(x, v);
}

template<class Vertex, class Value, class Topology, class Function>
void
reeber::
compute_merge_tree2(PathMergeTree<Vertex, Value>& mt, const Topology& topology, const Function& f)
{
    dlog::prof << "compute-merge-tree2";

    typedef     typename PathMergeTree<Vertex, Value>::Neighbor     Neighbor;

    auto vertices_ = topology.vertices();
    vector<Vertex> vertices(std::begin(vertices_), std::end(vertices_));

    for_each(0, vertices.size(), [&](size_t i) { Vertex a = vertices[i]; mt.add(a, f(a)); });

    for_each(0, vertices.size(), [&](size_t i)
    {
        Vertex a = vertices[i];
        Neighbor u = mt[a];

        for (const Vertex& b : topology.link(a))
        {
            if (b < a) continue;
            Neighbor v = mt[b];
            mt.merge(u, v);
        }
#ifndef REEBER_USE_TBB
        if (u->children.size() == 1)   // degree-2 node
        {
            // collapse node down
            u->children[0]->parent = u->parent;
            u->parent->children.push_back(u->children[0]);
            mt.nodes()[a] = u->children[0];
            mt.delete_node(u);
        }
#endif
    });

    dlog::prof >> "compute-merge-tree2";
}

template<class Vertex, class Value, class Edges>
void
reeber::
merge(PathMergeTree<Vertex, Value>& mt1, PathMergeTree<Vertex, Value>& mt2, const Edges& edges)
{
    dlog::prof << "merge";

    mt1.nodes_.insert(mt2.nodes().begin(), mt2.nodes().end());
    mt2.nodes_.clear();

    for_each(0, edges.size(), [&](size_t i)
    {
        Vertex a, b;
        std::tie(a, b) = edges[i];
        mt1.merge(mt1[a], mt1[b]);
    });

    dlog::prof >> "merge";
}
