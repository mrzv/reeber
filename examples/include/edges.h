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
    typedef     EdgeMap                               Block::*EdgePtr;

                EnqueueEdges(GridPtr grid_, BoxPtr local_, MtPtr mt_, EdgePtr edges_, bool wrap_):
                    grid(grid_), local(local_), mt(mt_), edges(edges_), wrap(wrap_)          {}

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

            std::unordered_map<Index, std::tuple<Value, Index>> relabel;

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
                        Index s = std::get<0>(u_node->parent)->vertex;
                        Index u_ = std::get<1>(u_node->parent)->vertex;
                        if (u != s) u_ = u;
                        relabel[u] = std::make_tuple(u_node->value, u_);
                        for (Vertex vp : expanded.position_link(u))
                        {
                            if (!(b->*local).contains(vp))
                            {
                                // ensure that vertex is inside of domain
                                vp = (b->*local).positive_position(vp);
                                Index v = (b->*local).position_to_vertex()(vp);
                                if ((b->*edges).find(std::make_tuple(u_, v)) != (b->*edges).end())
                                {
                                    Neighbor u_node_old = (b->*mt).node(std::get<1>((b->*edges)[std::make_tuple(u_, v)]));
                                    if ((b->*mt).cmp(u_node, u_node_old)) (b->*edges)[std::make_tuple(u_, v)] = std::make_tuple(u_node->value, u);
                                }
                                else (b->*edges)[std::make_tuple(u_, v)] = std::make_tuple(u_node->value, u);
                            }
                        }
                        ++it;
                    }
                }
            }

            for (int i = 0; i < l->size(); i++) cp.enqueue(l->target(i), relabel);
        }
    }

    GridPtr grid;
    BoxPtr  local;
    MtPtr   mt;
    EdgePtr edges;
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
    typedef     EdgeMap                               Block::*EdgePtr;


                DequeueEdges(GridPtr grid_, BoxPtr local_, MtPtr mt_, EdgePtr edges_):
                    grid(grid_), local(local_), mt(mt_), edges(edges_)          {}

    void        operator()(void* b_, const diy::Master::ProxyWithLink& cp, void*) const
    {
        typedef     diy::RegularGridLink                        RGLink;

        Block*          b = static_cast<Block*>(b_);
        RGLink*         l = static_cast<RGLink*>(cp.link());

        std::unordered_map<Index, std::tuple<Value, Index>> relabel;
        EdgeMap new_edges;
        for (int i = 0; i < l->size(); i++)
        {
            cp.dequeue(l->target(i).gid, relabel);
            for (auto& kv : b->*edges)
            {
                Index u_, v, u;
                Value u_val;
                std::tie(u_, v) = kv.first;
                std::tie(u_val, u) = kv.second;
                if (relabel.find(v) != relabel.end())
                {
                    Index v_;
                    Value v_val;
                    std::tie(v_val, v_) = relabel[v];

                    Node u_node, v_node;
                    u_node.vertex = u;
                    u_node.value = u_val;
                    v_node.vertex = v;
                    v_node.value = v_val;

                    if (new_edges.find(std::make_tuple(u_, v_)) != (b->*edges).end())
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
            }
        }
        (b->*edges).swap(new_edges);
    }

    GridPtr grid;
    BoxPtr  local;
    MtPtr   mt;
    EdgePtr edges;
};

#endif
