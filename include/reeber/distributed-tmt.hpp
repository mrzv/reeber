#include <reeber/range/filtered.h>

// local both in vertices and in link
template<class Topology, class LocalTest>
struct reeber::detail::LocalTopology
{
    using Vertices          = decltype(std::declval<Topology>().vertices());
    using FilteredVertices  = reeber::range::filtered_range<Vertices, LocalTest>;
    using Link              = decltype(std::declval<Topology>().link(*std::begin(std::declval<Vertices>())));
    using FilteredLink      = reeber::range::filtered_range<Link, LocalTest>;
    using Vertex            = typename std::remove_reference<decltype(*std::begin(std::declval<Vertices>()))>::type;

                LocalTopology(const Topology& topology_, const LocalTest& local_test_):
                    topology(topology_), local_test(local_test_)        {}

    FilteredVertices    vertices() const     { return topology.vertices() | reeber::range::filtered(local_test); }
    FilteredLink        link(Vertex v) const { return topology.link(v)    | reeber::range::filtered(local_test); }

    const Topology&     topology;
    const LocalTest&    local_test;
};

// TODO: this needs to use gids as a mechanism to decide what to prune, not boxes
// NB: this assumes binary reduction
template<class Block, class Vertex, class Value, class Partners, class GidGenerator>
struct reeber::detail::MergeSparsify
{
    using TripletMergeTree  = reeber::TripletMergeTree<Vertex,Value>;
    using EdgeMap           = reeber::EdgeMap<Vertex,Value>;
    using EdgeMaps          = reeber::EdgeMaps<Vertex,Value>;

    TripletMergeTree Block::*       tmt;
    EdgeMaps Block::*               edge_maps;
    const GidGenerator&             gid_generator;

            MergeSparsify(TripletMergeTree Block::* tmt_,
                          EdgeMaps Block::*         edge_maps_,
                          const GidGenerator&       gid_generator_):
                tmt(tmt_), edge_maps(edge_maps_),
                gid_generator(gid_generator_)       {}

    void    operator()(Block* b, const diy::ReduceProxy& srp, const Partners& partners) const
    {
        using Edge              = Edge<Vertex>;

        LOG_SEV(debug) << "Entered merge_sparsify()";

        unsigned round    = srp.round();
        LOG_SEV(debug) << "Round: " << round;
        LOG_SEV_IF(srp.master()->communicator().rank() == 0, info) << "round = " << srp.round();

        auto& edges = (b->*edge_maps)[srp.gid()];

        // receive trees, merge, and sparsify
        int in_size = srp.in_link().size();
        LOG_SEV(debug) << "  incoming link size: " << in_size;
        if (in_size)
        {
            dlog::prof << "dequeue";
            std::vector<TripletMergeTree>  trees;
            for (int i = 0; i < in_size; ++i)
                trees.emplace_back((b->*tmt).negate());
            EdgeMap out_edges;
            int local_pos = -1;
            for (int i = 0; i < in_size; ++i)
            {
              int nbr_gid = srp.in_link().target(i).gid;
              if (nbr_gid == srp.gid())
              {
                  local_pos = i;
                  trees[i].swap(b->*tmt);
                  LOG_SEV(debug) << "  swapped in tree of size: " << trees[i].size();
              } else
              {
                  srp.dequeue(nbr_gid, trees[i]);
                  srp.dequeue(nbr_gid, out_edges);
                  LOG_SEV(debug) << "  received tree of size: " << trees[i].size();
              }
            }
            LOG_SEV(debug) << "  trees and bounds received";
            dlog::prof >> "dequeue";

            dlog::prof << "compute edges";
            std::vector<Edge> merge_edges, discard;
            auto& mt = trees[local_pos];
            for (auto& kv : edges)
            {
                // Determine if (u,v) is an edge to this neighbor; if so, find
                // u - s - v, and determine whether s is local or remote (i.e.,
                // if we are mering (u,s) or (s,v))

                auto& uv = kv.first;
                Vertex u, v;
                std::tie(u,v) = uv;
                auto vu = std::make_tuple(v,u);
                if (out_edges.find(vu) != out_edges.end())
                {
                    map_erase(out_edges, vu);
                    Vertex s = std::get<1>(kv.second);

                    if (mt.contains(s))
                        merge_edges.emplace_back(s,v);
                    else
                        merge_edges.emplace_back(u,s);
                    discard.push_back(uv);
                }
            }
            for (Edge e : discard)
                map_erase(edges, e);
            edges.insert(out_edges.begin(), out_edges.end());
            dlog::prof >> "compute edges";

            trees[0].swap(b->*tmt);
            reeber::merge(b->*tmt, trees[1], merge_edges);

            trees.clear();
            LOG_SEV(debug) << "  trees merged: " << (b->*tmt).size();
        }

        dlog::prof << "compute edge_vertices";
        std::unordered_set<Vertex> edge_vertices;
        for (auto& kv : edges)
        {
            Vertex u_ = std::get<0>(kv.first);
            edge_vertices.insert(u_);
            Vertex s = std::get<1>(kv.second);
            edge_vertices.insert(s);
            // s might not be in our "global" domain, but that's Ok; we ignore that optimization
        }
        dlog::prof >> "compute edge_vertices";

        auto gid_gen    = gid_generator(b);
        int  gid        = srp.gid();
        auto local_test = [&gid_gen,gid](Vertex u)  { return gid_gen(u) == gid; };

        if (in_size)
            reeber::sparsify(b->*tmt, [&edge_vertices, &local_test](Vertex u)
                               { return local_test(u) || edge_vertices.find(u) != edge_vertices.end(); });

        // send (without the vertices) to the neighbors
        int out_size = srp.out_link().size();
        if (out_size == 0)        // final round: create the final local-global tree, nothing needs to be sent
        {
            LOG_SEV(debug) << "Sparsifying final tree of size: " << (b->*tmt).size();
            reeber::sparsify(b->*tmt, local_test);
            LOG_SEV(debug) << "[" << b->gid << "] " << "Final tree size: " << (b->*tmt).size();
            return;
        }

        TripletMergeTree mt_out((b->*tmt).negate());
        reeber::sparsify(mt_out, b->*tmt, [&edge_vertices](Vertex u) { return edge_vertices.find(u) != edge_vertices.end(); });

        dlog::prof << "enqueue";
        for (int i = 0; i < out_size; ++i)
        {
          diy::BlockID nbr_bid = srp.out_link().target(i);
          if (nbr_bid.gid != srp.gid())
          {
            srp.enqueue(nbr_bid, mt_out, &save_no_vertices);
            srp.enqueue(nbr_bid, edges);
          }
        }
        dlog::prof >> "enqueue";
    }

    static void save_no_vertices(diy::BinaryBuffer& bb, const TripletMergeTree& mt)
    {
        reeber::Serialization<TripletMergeTree>::save(bb, mt, false);
    }
};


// "Resolve" outgoing edges stored as keys in edge_maps[target_gid]. For each
// one, find its canonical representation (representatives in both blocks) as
// well as the value and vertex through which the connection is made
template<class Block, class Vertex, class Value>
void
reeber::
resolve_edges(diy::Master&                            master,
              TripletMergeTree<Vertex,Value> Block::* tmt_,
              EdgeMaps<Vertex,Value> Block::*         edge_maps_)
{
    using ValueVertex   = std::tuple<Value, Vertex>;
    using RelabelMap    = std::unordered_map<Vertex, ValueVertex>;
    using RelabelMaps   = std::unordered_map<int, RelabelMap>;
    using EdgeMap       = reeber::EdgeMap<Vertex,Value>;

    master.foreach([&](Block* b, const diy::Master::ProxyWithLink& cp)
    {
        auto&   tmt         = b->*tmt_;
        auto&   edge_maps   = b->*edge_maps_;

        size_t edge_count = 0;

        RelabelMaps relabel;
        for (auto& em_kv : edge_maps)
        {
            int   v_gid    = em_kv.first;
            auto& edge_map = em_kv.second;

            EdgeMap new_edge_map;
            for (auto& kv : edge_map)
            {
                Vertex u,v;
                std::tie(u,v) = kv.first;

                auto  u_node = tmt[u];
                auto  u_     = tmt.representative(u_node, u_node)->vertex;
                Value f_u    = u_node->value;

                ValueVertex uval { f_u, u };

                ++edge_count;

                relabel[v_gid][u] = { f_u, u_ };
                auto it = new_edge_map.find({ u_, v });
                if (it == new_edge_map.end())
                    new_edge_map.emplace(std::make_tuple(u_, v), uval);
                else if (tmt.cmp(uval, it->second))
                    it->second = uval;
            }
            edge_map.swap(new_edge_map);
        }

        auto* l = cp.link();
        for (int i = 0; i < l->size(); i++)
        {
            int target_gid = l->target(i).gid;
            cp.enqueue(l->target(i), relabel[target_gid]);
        }

        cp.collectives()->clear();
        cp.all_reduce(edge_count, std::plus<size_t>());
    });

    master.exchange();

    master.foreach([&](Block* b, const diy::Master::ProxyWithLink& cp)
    {
        auto&   tmt         = b->*tmt_;
        auto&   edge_maps   = b->*edge_maps_;
        auto&   edges       = edge_maps[cp.gid()];

        size_t edge_count = cp.get<size_t>();
        cp.scratch(edge_count);

        size_t  new_edge_count = 0;
        auto* l = cp.link();
        for (int i = 0; i < l->size(); i++)
        {
            int target_gid = l->target(i).gid;

            RelabelMap relabel;
            cp.dequeue(target_gid, relabel);

            EdgeMap new_edges;
            for (auto& kv : edge_maps[target_gid])
            {
                Vertex u, u_, v, v_;
                Value  f_u, f_v;

                std::tie(u_, v)   = kv.first;
                std::tie(f_u, u)  = kv.second;
                std::tie(f_v, v_) = relabel[v];

                ValueVertex uval { f_u, u },
                            vval { f_v, v };

                // ensure u < v
                if (tmt.cmp(vval, uval))
                    std::swap(uval, vval);

                auto u_v_ = std::make_tuple(u_,v_);
                auto it = new_edges.find(u_v_);
                if (it != new_edges.end())
                {
                    if (tmt.cmp(vval, it->second))
                        it->second = vval;
                } else
                    new_edges[u_v_] = vval;
            }
            edge_maps[target_gid].clear();
            edges.insert(new_edges.begin(), new_edges.end());
            new_edge_count += new_edges.size();
        }

        cp.all_reduce(new_edge_count, std::plus<size_t>());
    });
    master.exchange();      // process collectives: add up the edges
    for (unsigned i = 0; i < master.size(); ++i)
    {
        if (master.gid(i) == 0)
        {
            auto cp = master.proxy(i);
            auto edges_before = cp.get<size_t>();
            auto edges_after  = cp.get<size_t>();
            LOG_SEV(info) << "total edges before = " << edges_before << "; total edges after = " << edges_after;
        }
        master.proxy(i).collectives()->clear();
    }
}
