#include <dlog/log.h>
#include <dlog/stats.h>

#include "format.h"

template<class Vertex, class Value>
typename reeber::TripletMergeTree<Vertex, Value>::Neighbor
reeber::TripletMergeTree<Vertex, Value>::
add(const Vertex& x, Value v)
{
    Neighbor n = new_node();
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
    auto it = nodes().find(x);
    if (it != nodes().end())
        return it->second;
    else
        return add(x, v);
}

template<class Vertex, class Value>
typename reeber::TripletMergeTree<Vertex, Value>::Neighbor
reeber::TripletMergeTree<Vertex, Value>::
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
std::tuple<typename reeber::TripletMergeTree<Vertex, Value>::Neighbor, typename reeber::TripletMergeTree<Vertex, Value>::Neighbor>
reeber::TripletMergeTree<Vertex, Value>::
repair(const Neighbor u)
{
    typedef     typename TripletMergeTree::Neighbor         Neighbor;

    Neighbor s, v, s2, v2;
    Neighbor ov;
    do
    {
        std::tie(s, ov) = u->parent();
        v = ov;
        if (u == v) return std::make_tuple(s,v);
        std::tie(s2, v2) = v->parent();
        while (!cmp(s, s2) && v != v2)
        {
            v = v2;
            std::tie(s2, v2) = v->parent();
        }
    } while (!cas_link(u,s,ov,s,v));

    return std::make_tuple(s,v);
}

template<class Vertex, class Value>
typename reeber::TripletMergeTree<Vertex, Value>::Neighbor
reeber::TripletMergeTree<Vertex, Value>::
find_deepest(const Neighbor u)
{
    Neighbor u_ = u->cur_deepest;
    Neighbor v = std::get<1>(u_->parent());

    while (u_ != v)
    {
        u_ = v->cur_deepest;
        v = std::get<1>(u_->parent());
    }
    Neighbor d = u_;

    u_ = u->cur_deepest;
    v = std::get<1>(u_->parent());
    while (u_ != v)
    {
        u_ = v->cur_deepest;
        v->cur_deepest = d;
        v = std::get<1>(u_->parent());
    }

    u->cur_deepest = d;
    return d;
}


template<class Vertex, class Value, class Topology, class Function>
void
reeber::compute_merge_tree(TripletMergeTree<Vertex, Value>& mt, const Topology& topology, const Function& f)
{
    dlog::prof << "compute-merge-tree";

    typedef     typename TripletMergeTree<Vertex, Value>::Neighbor Neighbor;
    typedef     typename TripletMergeTree<Vertex, Value>::Node::ValueVertex ValueVertex;

    std::vector<ValueVertex>     vertices;
    vertices.reserve(topology.size());
    for (Vertex v : topology.vertices())
    {
        vertices.push_back(std::make_pair(f(v), v));
    }

    if (mt.negate())
        std::sort(vertices.begin(), vertices.end(), std::greater<ValueVertex>());
    else
        std::sort(vertices.begin(), vertices.end(), std::less<ValueVertex>());

    LOG_SEV(debug) << "Computing merge tree out of " << vertices.size() << " vertices";

    for (const ValueVertex& fx : vertices)
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

template<class Vertex, class Value, class Special>
void
reeber::remove_degree_two(TripletMergeTree<Vertex, Value>& mt, const Special& special)
{
    dlog::prof << "remove-degree-two";

    typedef     typename TripletMergeTree<Vertex, Value>::Neighbor          Neighbor;
    typedef     typename TripletMergeTree<Vertex, Value>::Node::ValueVertex ValueVertex;

    set<Vertex> keep;

    for_each_range(mt.nodes(), [&](const std::pair<Vertex,Neighbor>& n)
    {
        Neighbor u, s, v;
        u = n.second;
        std::tie(s, v) = u->parent();
        if (u != s || special(u->vertex))
        {
            while (1)
            {
                keep.insert(u->vertex);
                keep.insert(s->vertex);
                if (keep.find(v->vertex) != keep.end()) break;
                keep.insert(v->vertex);
                u = v;
                std::tie(s, v) = u->parent();
            }
        }
    });

    // Although the standard guarantees that this works only starting with
    // C++14, according to this issue, all compilers support it with C++11:
    // http://wg21.cmeerw.net/lwg/issue2356
    auto it = mt.nodes().begin();
    while (it != mt.nodes().end())
    {
        if (keep.find(it->first) == keep.end())
        {
            Neighbor u = it->second;
            Neighbor v = std::get<1>(u->parent());
            v->vertices.push_back(ValueVertex(u->value, u->vertex));
            mt.delete_node(it->second);
            it = map_erase(mt.nodes(), it);
        } else
            ++it;
    }

    dlog::prof >> "remove-degree-two";
}

template<class Vertex, class Value, class Topology, class Function>
void
reeber::compute_merge_tree2(TripletMergeTree<Vertex, Value>& mt, const Topology& topology, const Function& f)
{
    dlog::prof << "compute-merge-tree2";

    typedef     typename TripletMergeTree<Vertex, Value>::Neighbor        Neighbor;

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
            auto it = mt.nodes().find(b);
            Neighbor v = it->second;
            if (mt.cmp(v, u)) merge(mt, u, u, v);
            else merge(mt, v, v, u);
        }
    });

    for_each_range(mt.nodes(), [&](const std::pair<Vertex,Neighbor>& n) { mt.repair(n.second); });

    dlog::prof >> "compute-merge-tree2";
}

template<class Vertex, class Value, class Special>
reeber::set<Vertex>
reeber::sparsify_keep(TripletMergeTree<Vertex, Value>& mt, const Special& special)
{
    typedef     typename TripletMergeTree<Vertex, Value>::Neighbor        Neighbor;

    set<Vertex> keep;
    for_each_range(mt.nodes(), [&](const std::pair<Vertex,Neighbor>& n)
    {
        Neighbor u, s, v;
        if (special(n.first))
        {
            u = n.second;
            while (1)
            {
                std::tie(s, v) = u->parent();
                keep.insert(u->vertex);
                keep.insert(s->vertex);
                if (keep.find(v->vertex) != keep.end()) break;
                u = v;
            }
        }
    });
    return keep;
}

template<class Vertex, class Value, class Special>
void
reeber::sparsify(TripletMergeTree<Vertex, Value>& out, TripletMergeTree<Vertex, Value>& in, const Special& special)
{
    dlog::prof << "sparsify";

    typedef     typename TripletMergeTree<Vertex, Value>::Neighbor        Neighbor;

    set<Vertex> keep = sparsify_keep(in, special);

    for_each_range(keep, [&](Vertex x) { out.add(x, in[x]->value); });

    for_each_range(out.nodes(), [&](const std::pair<Vertex,Neighbor>& n)
    {
        Neighbor s,v;
        std::tie(s,v) = in[n.first]->parent();
        Neighbor ou = n.second;
        Neighbor os = out[s->vertex];
        Neighbor ov = out[v->vertex];
        out.link(ou, os, ov);
    });

    dlog::prof >> "sparsify";
}

template<class Vertex, class Value, class Special>
void
reeber::sparsify(TripletMergeTree<Vertex, Value>& mt, const Special& special)
{
    dlog::prof << "sparsify";

    set<Vertex> keep = sparsify_keep(mt, special);

    // Although the standard guarantees that this works only starting with
    // C++14, according to this issue, all compilers support it with C++11:
    // http://wg21.cmeerw.net/lwg/issue2356
    auto it = mt.nodes().begin();
    while (it != mt.nodes().end())
    {
        if (keep.find(it->first) == keep.end())
        {
            mt.delete_node(it->second);
            it = map_erase(mt.nodes(), it);
        } else
            ++it;
    }

    dlog::prof >> "sparsify";
}

template<class Vertex, class Value>
void
reeber::merge(TripletMergeTree<Vertex, Value>& mt, typename TripletMergeTree<Vertex, Value>::Neighbor u, typename TripletMergeTree<Vertex, Value>::Neighbor s, typename TripletMergeTree<Vertex, Value>::Neighbor v)
{
    typedef     typename TripletMergeTree<Vertex, Value>::Neighbor          Neighbor;

    Neighbor s_u, u_, s_v, v_;
    std::tie(s_u, u_) = mt.repair(u);
    std::tie(s_v, v_) = mt.repair(v);

    while (!mt.cmp(s, s_u) && s_u != u_)
    {
        u = u_;
        std::tie(s_u, u_) = mt.repair(u);
    }
    while (!mt.cmp(s, s_v) && s_v != v_)
    {
        v = v_;
        std::tie(s_v, v_) = mt.repair(v);
    }

    if (u == v) return;

    if (mt.cmp(v, u))
    {
        bool success = mt.cas_link(u, s_u, u_, s, v);
        if (!success)
            merge(mt, u, s, v);     // rinse and repeat
        else if (u != u_)
            merge(mt, v, s_u, u_);
    }
    else
    {
        bool success = mt.cas_link(v, s_v, v_, s, u);
        if (!success)
            merge(mt, u, s, v);     // rinse and repeat
        else if (v != v_)
            merge(mt, u, s_v, v_);
    }
}


template<class Vertex, class Value, class Edges>
void
reeber::merge(TripletMergeTree<Vertex, Value>& mt1, TripletMergeTree<Vertex, Value>& mt2, const Edges& edges)
{
    dlog::prof << "merge";

    typedef     typename TripletMergeTree<Vertex, Value>::Neighbor        Neighbor;

    mt1.nodes_.insert(mt2.nodes().begin(), mt2.nodes().end());
    mt2.nodes_.clear();

    for_each(0, edges.size(), [&](size_t i)
    {
        Vertex a, b, c;
        std::tie(a, b, c) = edges[i];
        Neighbor u = mt1.node(a);
        Neighbor s = mt1.node(b);
        Neighbor v = mt1.node(c);
        if (mt1.cmp(v, u)) merge(mt1, u, s, v);
        else merge(mt1, v, s, u);
    });

    for_each_range(mt1.nodes(), [&](const std::pair<Vertex,Neighbor>& n) { mt1.repair(n.second); });

    dlog::prof >> "merge";
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
        std::tie(s, v) = u->parent();
        if (u != s || u == v) f(u, s, v);
    }
}
