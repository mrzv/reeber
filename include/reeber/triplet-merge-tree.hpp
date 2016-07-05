#include <dlog/log.h>

#include "format.h"

template<class Vertex, class Value>
typename reeber::TripletMergeTree<Vertex, Value>::Neighbor
reeber::TripletMergeTree<Vertex, Value>::
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
typename reeber::TripletMergeTree<Vertex, Value>::Neighbor
reeber::TripletMergeTree<Vertex, Value>::
find_or_add(const Vertex& x, Value v)
{
    Neighbor n;
    auto it = nodes().find(x);
    if (it != nodes().end()) n = it->second;
    else n = add(x, v);
    return n;
}

template<class Vertex, class Value>
void
reeber::TripletMergeTree<Vertex, Value>::
repair(const Neighbor u)
{
    typedef     typename TripletMergeTree::Neighbor         Neighbor;

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
typename reeber::TripletMergeTree<Vertex, Value>::Neighbor
reeber::TripletMergeTree<Vertex, Value>::
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


template<class Vertex, class Value, class Topology, class Function>
void
reeber::compute_merge_tree(TripletMergeTree<Vertex, Value>& mt, const Topology& topology, const Function& f)
{
    dlog::prof << "compute-merge-tree";
    typedef     std::pair<Value, Vertex>                           ValueVertexPair;
    typedef     typename TripletMergeTree<Vertex, Value>::Neighbor Neighbor;

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

template<class Vertex, class Value, class Topology, class Function>
void
reeber::compute_merge_tree2(TripletMergeTree<Vertex, Value>& mt, const Topology& topology, const Function& f)
{
    typedef     typename TripletMergeTree<Vertex, Value>::Neighbor        Neighbor;

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

template<class Vertex, class Value, class Special>
void
reeber::sparsify(TripletMergeTree<Vertex, Value>& out, TripletMergeTree<Vertex, Value>& in, const Special& special)
{
    dlog::prof << "sparsify";

    typedef     typename TripletMergeTree<Vertex, Value>::Neighbor        Neighbor;

    Neighbor u, s, v, u_, s_, v_;

    for (auto n : in.nodes())
    {
        if (special(n.first))
        {
            u = n.second;
            while (1)
            {
                std::tie(s, v) = u->parent;
                u_ = out.find_or_add(u->vertex, u->value);
                s_ = out.find_or_add(s->vertex, s->value);
                auto it = out.nodes().find(v->vertex);
                if (it != out.nodes().end())
                {
                    v_ = it->second;
                    out.link(u_, s_, v_);
                    out.link(s_, s_, v_);
                    break;
                }
                else
                {
                    v_ = out.add(v->vertex, v->value);
                    out.link(u_, s_, v_);
                    out.link(s_, s_, v_);
                    u = v;
                }
            }
        }
    }

    dlog::prof >> "sparsify";
}

template<class Vertex, class Value, class Special>
void
reeber::sparsify(TripletMergeTree<Vertex, Value>& mt, const Special& special)
{
    typedef     typename TripletMergeTree<Vertex, Value>::Neighbor        Neighbor;

    std::unordered_set<Vertex> discard;
    Neighbor u, s, v;

    for (auto n: mt.nodes()) discard.insert(n.first);

    for (auto n : mt.nodes())
    {
        if (special(n.first))
        {
            u = n.second;
            while (1)
            {
                std::tie(s, v) = u->parent;
                discard.erase(u->vertex);
                discard.erase(s->vertex);
                if (discard.find(v->vertex) == discard.end()) break;
                discard.erase(v->vertex);
                u = v;
            }
        }
    }

    for (Vertex x : discard)
    {
        delete mt.node(x);
        mt.nodes().erase(x);
    }
}

template<class Vertex, class Value>
void
reeber::merge(TripletMergeTree<Vertex, Value>& mt, typename TripletMergeTree<Vertex, Value>::Neighbor u, typename TripletMergeTree<Vertex, Value>::Neighbor s, typename TripletMergeTree<Vertex, Value>::Neighbor v)
{
    typedef     typename TripletMergeTree<Vertex, Value>::Neighbor          Neighbor;

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


template<class Vertex, class Value, class Edges>
void
reeber::merge(TripletMergeTree<Vertex, Value>& mt1, TripletMergeTree<Vertex, Value>& mt2, const Edges& edges)
{
    typedef     typename TripletMergeTree<Vertex, Value>::Neighbor        Neighbor;

    mt1.nodes_.insert(mt2.nodes().begin(), mt2.nodes().end());
    mt2.nodes_.clear();

    Vertex a, b;
    Neighbor u, v;
    for (std::tuple<Vertex, Vertex> edge : edges)
    {
        std::tie(a, b) = edge;
        u = mt1.node(a);
        v = mt1.node(b);
        if (mt1.cmp(v, u)) merge(mt1, u, u, v);
        else merge(mt1, v, v, u);
    }

    for (auto n : mt1.nodes()) mt1.repair(n.second);
}


template<class Vertex, class Value, class Functor>
void
reeber::traverse_persistence(const TripletMergeTree<Vertex, Value>& mt, const Functor& f)
{
    typedef     typename TripletMergeTree<Vertex, Value>::Neighbor        Neighbor;

    Neighbor u, s, v;
    for (auto x : mt.nodes())
    {
        u = x.second;
        std::tie(s, v) = u->parent;
        if (u != s || u == v) f(u, s, v);
    }
}
