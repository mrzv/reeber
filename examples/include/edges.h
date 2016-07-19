#ifndef REEBER_EDGES_H
#define REEBER_EDGES_H

#include <diy/link.hpp>

#include <reeber/box.h>
#include <reeber/grid.h>

template<class Block_>
struct EnqueueEdges
{
    typedef     Block_                                Block;
    typedef     typename Block::OffsetGrid            Grid;
    typedef     typename Block::Box                   Box;
    typedef     typename Grid::Index                  Index;
    typedef     typename Grid::Vertex                 Vertex;
    typedef     typename Grid::Value                  Value;
    typedef     typename Block::TripletMergeTree      TripletMergeTree;
    typedef     typename Block::EdgeMap               EdgeMap;
    typedef     typename TripletMergeTree::Neighbor   Neighbor;
    typedef     typename TripletMergeTree::Node       Node;
    typedef     Grid                                  Block::*GridPtr;
    typedef     Box                                   Block::*BoxPtr;
    typedef     TripletMergeTree                      Block::*MtPtr;
    typedef     std::vector<EdgeMap>                  Block::*EdgeMapsPtr;

                EnqueueEdges(GridPtr grid_, BoxPtr local_, MtPtr mt_, EdgeMapsPtr edge_maps_, bool wrap_):
                    grid(grid_), local(local_), mt(mt_), edge_maps(edge_maps_), wrap(wrap_)          {}

    void        operator()(void* b_, const diy::Master::ProxyWithLink& cp, void*) const
    {
        typedef     diy::RegularGridLink                        RGLink;

        Block*      b = static_cast<Block*>(b_);
        RGLink*     l = static_cast<RGLink*>(cp.link());

        {
            Box expanded(b->*local);
            if (!wrap)
            {
                for (unsigned i = 0; i < Box::dimension(); ++i)
                {
                    if (expanded.from()[i] > 0) expanded.from()[i]--;
                    if (expanded.to()[i] < (b->*local).grid_shape()[i] - 1) expanded.to()[i]++;
                }
            }
            else
            {
                expanded.from() -= Vertex::one();
                expanded.to() += Vertex::one();
            }

            std::vector<std::unordered_map<Index, std::tuple<Value, Index>>> relabel(l->size());
            (b->*edge_maps).resize(l->size());

            for (unsigned axis = 0; axis < Box::dimension(); ++axis)
            {
                for (bool upper : {false, true})
                {
                    Box side = (b->*local).side(axis, upper);
                    r::VerticesIterator<Vertex> it = r::VerticesIterator<Vertex>::begin(side.from(), side.to()),
                                                end = r::VerticesIterator<Vertex>::end(side.from(), side.to());
                    while (it != end)
                    {
                        Index u = (b->*local).position_to_vertex()(*it);
                        Neighbor u_node = (b->*mt).node(u);
                        Index s = std::get<0>(u_node->parent())->vertex;
                        Index u_ = std::get<1>(u_node->parent())->vertex;
                        if (u != s) u_ = u;
                        for (Vertex vp : expanded.position_link(u))
                        {
                            // ensure that vertex is inside of domain
                            vp = (b->*local).positive_position(vp);
                            Index v = (b->*local).position_to_vertex()(vp);
                            if ((b->*local).contains(vp)) continue;
                            int neighbor = -1;
                            for (int i = 0; i < l->size(); i++)
                            {
                                auto bounds = l->bounds(i);
                                auto shape = (b->*local).grid_shape();
                                Box box(shape, bounds.min, bounds.max);
                                auto vp = box.position(v);
                                if (box.contains(vp))
                                {
                                    relabel[i][u] =  std::make_tuple(u_node->value, u_);
                                    neighbor = i;
                                    break;
                                }
                            }
                            if (neighbor != -1)
                            {
                                if ((b->*edge_maps)[neighbor].find(std::make_tuple(u_, v)) != (b->*edge_maps)[neighbor].end())
                                {
                                    Neighbor u_node_old = (b->*mt).node(std::get<1>((b->*edge_maps)[neighbor][std::make_tuple(u_, v)]));
                                    if ((b->*mt).cmp(u_node, u_node_old)) (b->*edge_maps)[neighbor][std::make_tuple(u_, v)] = std::make_tuple(u_node->value, u);
                                }
                                else (b->*edge_maps)[neighbor][std::make_tuple(u_, v)] = std::make_tuple(u_node->value, u);
                            }
                            else LOG_SEV(warning) << "Didn't find neighbor";
                        }
                        ++it;
                    }
                }
            }

            for (int i = 0; i < l->size(); i++) cp.enqueue(l->target(i), relabel[i]);
        }
    }

    GridPtr grid;
    BoxPtr  local;
    MtPtr   mt;
    EdgeMapsPtr edge_maps;
    bool wrap;
};

template<class Block_>
struct DequeueEdges
{
    typedef     Block_                                Block;
    typedef     typename Block::OffsetGrid            Grid;
    typedef     typename Block::Box                   Box;
    typedef     typename Grid::Index                  Index;
    typedef     typename Grid::Vertex                 Vertex;
    typedef     typename Grid::Value                  Value;
    typedef     typename Block::TripletMergeTree      TripletMergeTree;
    typedef     typename Block::EdgeMap               EdgeMap;
    typedef     typename TripletMergeTree::Neighbor   Neighbor;
    typedef     typename TripletMergeTree::Node       Node;
    typedef     Grid                                  Block::*GridPtr;
    typedef     Box                                   Block::*BoxPtr;
    typedef     TripletMergeTree                      Block::*MtPtr;
    typedef     std::vector<EdgeMap>                  Block::*EdgeMapsPtr;
    typedef     EdgeMap                               Block::*EdgePtr;


                DequeueEdges(GridPtr grid_, BoxPtr local_, MtPtr mt_, EdgeMapsPtr edge_maps_, EdgePtr edges_):
                    grid(grid_), local(local_), mt(mt_), edge_maps(edge_maps_), edges(edges_)          {}

    void        operator()(void* b_, const diy::Master::ProxyWithLink& cp, void*) const
    {
        typedef     diy::RegularGridLink                        RGLink;

        Block*          b = static_cast<Block*>(b_);
        RGLink*         l = static_cast<RGLink*>(cp.link());

        std::unordered_map<Index, std::tuple<Value, Index>> relabel;
        for (int i = 0; i < l->size(); i++)
        {
            EdgeMap new_edges;
            cp.dequeue(l->target(i).gid, relabel);
            for (auto& kv : (b->*edge_maps)[i])
            {
                Index u, u_, v, v_;
                Value u_val, v_val;
                std::tie(u_, v) = kv.first;
                std::tie(u_val, u) = kv.second;
                std::tie(v_val, v_) = relabel[v];

                Node u_node, v_node;
                u_node.vertex = u;
                u_node.value = u_val;
                v_node.vertex = v;
                v_node.value = v_val;

                if (new_edges.find(std::make_tuple(u_, v_)) != new_edges.end())
                {
                    Node cur_node;
                    std::tie(cur_node.value, cur_node.vertex) = new_edges[std::make_tuple(u_, v_)];
                    if ((b->*mt).cmp(u_node, v_node) && (b->*mt).cmp(v_node, cur_node))
                        new_edges[std::make_tuple(u_, v_)] = std::make_tuple(v_node.value, v_node.vertex);
                    else if ((b->*mt).cmp(v_node, u_node) && (b->*mt).cmp(u_node, cur_node))
                        new_edges[std::make_tuple(u_, v_)] = std::make_tuple(u_node.value, u_node.vertex);
                }
                else
                {
                    if ((b->*mt).cmp(u_node, v_node)) new_edges[std::make_tuple(u_, v_)] = std::make_tuple(v_node.value, v_node.vertex);
                    else new_edges[std::make_tuple(u_, v_)] = std::make_tuple(u_node.value, u_node.vertex);
                }
            }
            (b->*edge_maps)[i].clear();
            (b->*edges).insert(new_edges.begin(), new_edges.end());
        }
    }

    GridPtr grid;
    BoxPtr  local;
    MtPtr   mt;
    EdgeMapsPtr edge_maps;
    EdgePtr edges;
};

#endif
