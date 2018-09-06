#pragma once

#include <string>
#include <map>

#include "fab-tmt-block.h"

template<class Vertex, class Real>
using LocalIntegral = std::map<Vertex, Real>;

template<class Real, unsigned D>
LocalIntegral<reeber::AmrVertexId, Real> get_local_integral(FabTmtBlock<Real, D>* b, Real rho_min, Real rho_max, const std::string& fname)
{
    using AmrVertexId = reeber::AmrVertexId;
    using Block = FabTmtBlock<Real, D>;
    using Node = typename Block::Node;

    Real scaling_factor = b->scaling_factor();

    auto deepest_vertices = b->get_final_deepest_vertices();
    std::map<AmrVertexId, Real> local_integral;
    const bool negate = b->negate_;
    const auto& nodes = b->get_merge_tree().nodes();
    for(const auto& vertex_node_pair : nodes)
    {
        AmrVertexId current_vertex = vertex_node_pair.first;

        assert(current_vertex == vertex_node_pair.second->vertex);

        if (current_vertex.gid != b->gid)
            continue;

        Node* current_node = vertex_node_pair.second;
        AmrVertexId root = b->final_vertex_to_deepest_[current_vertex];
        // save only information about local vertices

        Real root_value = nodes.at(root)->value;
        if ((negate and root_value < rho_max) or (not negate and root_value > rho_min))
            continue;
        local_integral[root] += scaling_factor * current_node->value;
        for(const auto& value_vertex_pair : current_node->vertices) {
            local_integral[root] += scaling_factor * value_vertex_pair.first;
        }
    }
    return local_integral;
}