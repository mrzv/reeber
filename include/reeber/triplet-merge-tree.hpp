#include <dlog/log.h>
#include <dlog/stats.h>

#include "format.h"

template<class Vertex, class Value>
typename reeber::TripletMergeTree<Vertex, Value>::Neighbor
reeber::TripletMergeTree<Vertex, Value>::
add(const Vertex& x, Value v)
{
    Neighbor n = new_node();
    n->vertex = x;
    n->value = v;
    n->cur_deepest = n;
    link(n, n, n);
    nodes_.emplace(x,n);
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
typename reeber::TripletMergeTree<Vertex, Value>::Neighbor
reeber::TripletMergeTree<Vertex, Value>::
representative(Neighbor u, Neighbor a) const
{
    Neighbor s, v;
    std::tie(s, v) = u->parent();
    while (!cmp(a, s) && s != v)
    {
        u = v;
        std::tie(s, v) = u->parent();
    }
    return u;
}

template<class Vertex, class Value>
std::tuple<typename reeber::TripletMergeTree<Vertex, Value>::Neighbor, typename reeber::TripletMergeTree<Vertex, Value>::Neighbor>
reeber::TripletMergeTree<Vertex, Value>::
repair(const Neighbor u)
{
    Neighbor s, v, ov;
    do
    {
        std::tie(s, ov) = u->parent();
        v = representative(u, s);
        if (u == v) return std::make_tuple(s,v);
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
        if (u == v) keep.insert(u->vertex);
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

template<class Vertex, class Value>
void
reeber::repair(TripletMergeTree<Vertex, Value>& mt)
{
    using Neighbor = typename TripletMergeTree<Vertex, Value>::Neighbor;
    for_each_range(mt.nodes(), [&](const std::pair<Vertex,Neighbor>& n) { mt.repair(n.second); });
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
            Neighbor v = mt[b];
            mt.merge(u, v);
        }
    });

    repair(mt);

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
        Neighbor s, v;
        if (special(n.first))
        {
            Neighbor u = n.second;
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
reeber::TripletMergeTree<Vertex, Value>::
merge(Neighbor u, Neighbor s, Neighbor v)
{
    while(true)
    {
        u = representative(u, s);
        v = representative(v, s);
        if (u == v)
            break;

        Neighbor s_u, u_;
        Neighbor s_v, v_;
        std::tie(s_u, u_) = u->parent();
        std::tie(s_v, v_) = v->parent();

        // check that s_u and s_v haven't changed since running representative
        if (s_u != u_ && !cmp(s, s_u))
            continue;
        if (s_v != v_ && !cmp(s, s_v))
            continue;

        if (cmp(v, u))
        {
            std::swap(u, v);
            std::swap(s_u, s_v);
            std::swap(u_, v_);
        }

        bool success = cas_link(v, s_v, v_, s, u);
        if (success)
        {
            if (v == v_)
                break;

            s = s_v;
            v = v_;
        } // else: rinse and repeat
    }
}

template<class Vertex, class Value>
void
reeber::TripletMergeTree<Vertex, Value>::
merge(Neighbor u, Neighbor v)
{
    if (cmp(u, v))
        merge(v, v, u);
    else
        merge(u, u, v);
}

template<class Vertex, class Value, class Edges>
void
reeber::merge(TripletMergeTree<Vertex, Value>& mt1, TripletMergeTree<Vertex, Value>& mt2, const Edges& edges)
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

    repair(mt1);

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
