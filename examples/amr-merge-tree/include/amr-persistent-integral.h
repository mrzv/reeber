#pragma once

#include <string>
#include <map>

#include "fab-tmt-block.h"

template<class Real>
Real contribution_to_integral(Real v, Real rho_min, Real rho_max)
{
    if (v >= rho_min and v <= rho_max)
        return v;
    else
        return 0;
}

template<class Real>
using LocalIntegral = std::map<reeber::AmrVertexId, Real>;

template<class Real, unsigned D>
LocalIntegral<Real> get_local_integral(FabTmtBlock<Real, D>* b, Real rho_min, Real rho_max, const std::string& fname)
{
    using AmrVertexId = reeber::AmrVertexId;
    using Block = FabTmtBlock<Real, D>;
    using Node = typename Block::Node;
    using Neighbor = typename Block::Neighbor;

    auto deepest_vertices = b->get_current_deepest_vertices();
    std::map<AmrVertexId, Real> local_integral;
    const bool negate = b->negate_;
    const auto& nodes = b->get_merge_tree().nodes();
    for(const auto& vertex_node_pair : nodes)
    {
        AmrVertexId current_vertex = vertex_node_pair.first;
        Node* current_node = vertex_node_pair.second;
        AmrVertexId root = b->vertex_to_deepest_[current_vertex];
        // save only information about local components
        if (root.gid != b->gid)
            continue;
        Real root_value = nodes.at(root)->value;
        if ((negate and root_value < rho_max) or (not negate and root_value > rho_min))
            continue;
        local_integral[root] += contribution_to_integral(current_node->value, rho_min, rho_max);
        for(const auto& value_vertex_pair : current_node->vertices) {
            local_integral[root] += contribution_to_integral(value_vertex_pair.first, rho_min, rho_max);
        }
    }
    return local_integral;
}