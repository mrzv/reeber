#include <dlog/log.h>

#include "format.h"

template<class Vertex, class Value>
typename reeber::MergeTree<Vertex, Value>::Neighbor
reeber::MergeTree<Vertex, Value>::
add(const Vertex& x, Value v)
{
    Neighbor n = new Node;
    nodes_[x] = n;
    n->vertex = x;
    n->value = v;
    n->cur_deepest = n;
    link(n, n, n);
    return n;
}

template<class Vertex, class Value>
void
reeber::MergeTree<Vertex, Value>::
repair(const Neighbor u)
{
    typedef     typename MergeTree::Neighbor         Neighbor;

    Neighbor s, v, s2, v2;
    std::tie(s, v) = u->parent;
    if (u == v) return;
    std::tie(s2, v2) = v->parent;
    while (!cmp(s, s2) && v != v2)
    {
        v = v2;
        std::tie(s2, v2) = v->parent;
    }
    link(u, s, v);
}

template<class Vertex, class Value>
typename reeber::MergeTree<Vertex, Value>::Neighbor
reeber::MergeTree<Vertex, Value>::
find_deepest(const Neighbor u)
{
    Neighbor u_ = u->cur_deepest;
    Neighbor v = std::get<1>(u_->parent);

    while (u_ != v)
    {
        u_ = v->cur_deepest;
        v = std::get<1>(u_->parent);
    }
    Neighbor d = u_;

    u_ = u->cur_deepest;
    v = std::get<1>(u_->parent);
    while (u_ != v)
    {
        u_ = v->cur_deepest;
        v->cur_deepest = d;
        v = std::get<1>(u_->parent);
    }

    u->cur_deepest = d;
    return d;
}


template<class MergeTree, class Topology, class Function>
void
reeber::compute_merge_tree(MergeTree& mt, const Topology& topology, const Function& f)
{
    dlog::prof << "compute-merge-tree";
    typedef     typename Topology::Vertex       Vertex;
    typedef     typename Function::Value        Value;
    typedef     std::pair<Value, Vertex>        ValueVertexPair;
    typedef     typename MergeTree::Neighbor    Neighbor;

    std::vector<ValueVertexPair>     vertices;
    vertices.reserve(topology.size());
    for (Vertex v : topology.vertices())
    {
        vertices.push_back(std::make_pair(f(v), v));
    }

    if (mt.negate())
        std::sort(vertices.begin(), vertices.end(), std::greater<ValueVertexPair>());
    else
        std::sort(vertices.begin(), vertices.end(), std::less<ValueVertexPair>());

    LOG_SEV(debug) << "Computing merge tree out of " << vertices.size() << " vertices";

    for (const ValueVertexPair& fx : vertices)
    {
        Value val; Vertex x;
        std::tie(val, x) = fx;

        Neighbor u = mt.add(x, val);

        std::vector<Neighbor> leaves;
        for (const Vertex& y : topology.link(x))
        {
            auto it = mt.nodes().find(y);
            if (it != mt.nodes().end()) leaves.push_back(mt.find_deepest(it->second));
        }
        if (!leaves.empty())
        {
            Neighbor oldest = *std::min_element(leaves.begin(), leaves.end(), [&mt](Neighbor x, Neighbor y) { return mt.cmp(x, y); });
            mt.link(u, u, oldest);
            if (leaves.size() > 1)
            {
                for (const Neighbor v : leaves) if (v != oldest) mt.link(v, u, oldest);
            }
        }
    }

    dlog::prof >> "compute-merge-tree";
}

template<class MergeTree, class Topology, class Function>
void
reeber::compute_merge_tree2(MergeTree& mt, const Topology& topology, const Function& f)
{
    typedef     typename MergeTree::Vertex          Vertex;
    typedef     typename MergeTree::Neighbor        Neighbor;

    Neighbor u, v;
    for (Vertex a : topology.vertices())
    {
        u = mt.add(a, f(a));

        for (const Vertex& b : topology.link(a))
        {
            auto it = mt.nodes().find(b);
            if (it != mt.nodes().end())
            {
                v = it->second;
                if (mt.cmp(v, u)) merge(mt, u, u, v);
                else merge(mt, v, v, u);
            }
        }
    }

    for (auto n : mt.nodes()) mt.repair(n.second);
}

template<class MergeTree, class Special>
void
reeber::sparsify(MergeTree& out, const MergeTree& in, const Special& special)
{
    dlog::prof << "sparsify";

    typedef     typename MergeTree::Vertex          Vertex;
    typedef     typename MergeTree::Neighbor        Neighbor;


    Vertex s, v;

    for (auto n : in.nodes())
    {
        if (special(n.first))
        {
            Neighbor u = n.second;
            while (!in.contains(u->vertex))
            {
                std::tie(s, v) = u->parent;
                out.add(u->vertex, u->value);
                out.add(s->vertex, s->value);
                out.add(v->vertex, v->value);
                out.link(u->vertex, s->vertex, v->vertex);
                u = v;
            }
        }
    }

    dlog::prof >> "sparsify";
}


template<class MergeTree>
void
reeber::merge(MergeTree& mt, typename MergeTree::Neighbor u, typename MergeTree::Neighbor s, typename MergeTree::Neighbor v)
{
    typedef     typename MergeTree::Neighbor          Neighbor;

    Neighbor s_u, u_, s_v, v_;
    // mt.repair(u);
    std::tie(s_u, u_) = u->parent;
    // mt.repair(v);
    std::tie(s_v, v_) = v->parent;

    while (!mt.cmp(s, s_u) && s_u != u_)
    {
        u = u_;
        // mt.repair(u);
        std::tie(s_u, u_) = u->parent;
    }
    while (!mt.cmp(s, s_v) && s_v != v_)
    {
        v = v_;
        // mt.repair(v);
        std::tie(s_v, v_) = v->parent;
    }

    if (u == v) return;

    if (mt.cmp(v, u))
    {
        mt.link(u, s, v);
        if (u != u_) merge(mt, v, s_u, u_);
    }
    else
    {
        mt.link(v, s, u);
        if (v != v_) merge(mt, u, s_v, v_);
    }
}


template<class MergeTree, class Edges>
void
reeber::merge(MergeTree& mt1, const MergeTree& mt2, const Edges& edges)
{
    typedef     typename MergeTree::Vertex          Vertex;
    typedef     typename MergeTree::Neighbor        Neighbor;

    mt1.nodes_.insert(mt2.nodes().begin(), mt2.nodes().end());

    Vertex a, b;
    Neighbor u, v;
    for (std::tuple<Vertex, Vertex> edge : edges)
    {
        std::tie(a, b) = edge;
        u = mt1.node(a);
        v = mt2.node(b);
        if (mt1.cmp(v, u)) merge(mt1, u, u, v);
        else merge(mt1, v, v, u);
    }

    for (auto n : mt1.nodes()) mt1.repair(n.second);
}


template<class MergeTree, class Functor>
void
reeber::traverse_persistence(const MergeTree& mt, const Functor& f)
{
    typedef     typename MergeTree::Neighbor        Neighbor;

    Neighbor u, s, v;
    for (auto x : mt.nodes())
    {
        u = x.second;
        std::tie(s, v) = u->parent;
        if (u != s || u == v) f(u->vertex, s->vertex, v->vertex);
    }
}
