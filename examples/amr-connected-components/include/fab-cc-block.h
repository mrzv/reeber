#pragma once

#include <utility>
#include <numeric>
#include <stack>
#include <boost/functional/hash.hpp>

#include "diy/serialization.hpp"
#include "diy/grid.hpp"
#include "diy/link.hpp"
#include "diy/fmt/format.h"
#include "diy/fmt/ostream.h"
#include "diy/point.hpp"
#include <diy/master.hpp>

#include "disjoint-sets.h"

#include "reeber/amr-vertex.h"
#include "reeber/triplet-merge-tree.h"
#include "reeber/triplet-merge-tree-serialization.h"
#include "reeber/grid.h"
#include "reeber/grid-serialization.h"
#include "reeber/masked-box.h"
#include "reeber/edges.h"

#include "../../amr-merge-tree/include/fab-block.h"
#include "reeber/amr_helper.h"

#include "fab-connected-component.h"

namespace r = reeber;

template<class Real, unsigned D>
struct FabComponentBlock {
    using Shape = diy::Point<int, D>;

    using Grid = r::Grid<Real, D>;
    using GridRef = r::GridRef<Real, D>;

    using Component = FabConnectedComponent<Real>;

    using TripletMergeTree = typename Component::TripletMergeTree;
    using AmrVertexId = typename Component::AmrVertexId;
    using VertexValue = typename Component::VertexValue;
    using Value = typename Grid::Value;
    using MaskedBox = r::MaskedBox<D>;
    using Vertex = typename MaskedBox::Position;

    using AmrVertexContainer = std::vector<AmrVertexId>;
    using AmrVertexSet = typename Component::AmrVertexSet;

    using VertexValueMap = typename Component::VertexValueMap;

    using AmrEdge = r::AmrEdge;
    using AmrEdgeContainer = r::AmrEdgeContainer;
    using AmrEdgeSet = std::set<AmrEdge>;
    using VertexEdgesMap = std::map<AmrVertexId, AmrEdgeContainer>;

    using GidSet = typename Component::GidSet;
    using GidVector = std::vector<int>;

    using RealType = Real;

//    using UnionFind = typename Component::UnionFind;
    using VertexVertexMap = std::unordered_map<AmrVertexId, AmrVertexId>;
//    using VertexSizeMap = typename UnionFind::VertexSizeMap;

    using Neighbor = typename TripletMergeTree::Neighbor;
    using Node = typename TripletMergeTree::Node;
    using VertexNeighborMap =typename TripletMergeTree::VertexNeighborMap;

    using DiagramPoint = std::pair<Real, Real>;
    using Diagram = std::vector<DiagramPoint>;

    using LocalIntegral = std::unordered_map<AmrVertexId, Real>;
    // data

    int gid;
    MaskedBox local_;
    GridRef fab_;

    // if relative threshold is given, we cannot determine
    // LOW values in constructor. Instead, we mark all unmasked vertices
    // ACTIVE and save the average of their values in local_sum_ and the number in local_n_unmasked_
    // Pointer to grid data is saved in fab_ and after all blocks exchange their local averages
    // we resume initialization
    Real sum_{0};
    size_t n_unmasked_{0};
    std::unordered_map<AmrVertexId, Real> vertex_values_;

//    UnionFind disjoint_sets_;   // keep topology of graph of connected components
    std::vector<Component> components_;

    VertexVertexMap vertex_to_deepest_;

    diy::DiscreteBounds domain_;

    int done_{0};

    std::unordered_map<int, AmrEdgeContainer> gid_to_outgoing_edges_;

    VertexValueMap original_integral_values_;
    VertexValueMap global_integral_;

    std::unordered_map<AmrVertexId, Diagram> local_diagrams_;

    bool negate_;

    TripletMergeTree merge_tree_;
//    TripletMergeTree original_tree_;
    // tracking how connected components merge - disjoint sets data structure

    int round_{0};

    LocalIntegral local_integral_;

    // methods

    // simple getters/setters
    const diy::DiscreteBounds& domain() const
    { return domain_; }

    int refinement() const
    { return local_.refinement(); }

    int level() const
    { return local_.level(); }

    FabComponentBlock(diy::GridRef<Real, D>& fab_grid,
                      int _ref,
                      int _level,
                      const diy::DiscreteBounds& _domain,
                      const diy::DiscreteBounds& bounds,
                      const diy::DiscreteBounds& core,
                      int _gid,
                      diy::AMRLink *amr_link,
                      Real rho,                                           // threshold for LOW value
                      bool _negate,
                      bool is_absolute_threshold);

    FabComponentBlock() :
            fab_(nullptr, diy::Point<int, D>::zero())
    {}

    void init(Real absolute_rho, diy::AMRLink *amr_link);

    // compare w.r.t negate_ flag
    bool cmp(Real a, Real b) const;

    void set_low(const diy::Point<int, D>& v_bounds,
                 const Real& absolute_rho);

    void set_mask(const diy::Point<int, D>& v_bounds,
                  diy::AMRLink *l,
                  const Real& rho,
                  bool is_absolute_threshold);

    void compute_outgoing_edges(diy::AMRLink *l, VertexEdgesMap& vertex_to_outgoing_edges);

    void compute_original_connected_components(const VertexEdgesMap& vertex_to_outgoing_edges);

    void delete_low_edges(int sender_gid, AmrEdgeContainer& edges_from_sender,
                          const VertexVertexMap& received_vertex_to_deepest);

//
    void adjust_outgoing_edges();

    void sparsify_prune_original_tree() {}

    int get_n_components_for_gid(int gid) const;

    int are_all_components_done() const;

    bool is_deepest_computed(const AmrVertexId& v) const;

    std::vector<AmrVertexId> get_current_deepest_vertices() const;

    void compute_final_connected_components();

//    std::unordered_map<AmrVertexId, AmrVertexId> compute_connectivity(const AmrVertexSet& deepest);
    void update_connectivity(const AmrVertexContainer& deepest);

    AmrVertexContainer component_of(AmrVertexId deepest);

    void compute_local_integral(Real theta);

    bool check_symmetry(int gid, const std::vector<Component>& received_components);

    TripletMergeTree& get_merge_tree()
    { return merge_tree_; }

    const TripletMergeTree& get_merge_tree() const
    { return merge_tree_; }

    Real scaling_factor() const;

    Component& get_component_by_deepest(const AmrVertexId& deepest)
    {
        auto res_iter = std::find_if(components_.begin(), components_.end(),
                                     [deepest](const Component& c) { return c.original_deepest() == deepest; });
        if (res_iter == components_.end())
        {
            throw std::runtime_error("error in find_componenent, bad deepest");
        }
        return *res_iter;
    }

    void sanity_check_fin() const;

    static void *create()
    {
        return new FabComponentBlock;
    }

    static void destroy(void *b)
    {
        delete static_cast<FabComponentBlock *>(b);
    }

    static void save(const void *b, diy::BinaryBuffer& bb);

    static void load(void *b, diy::BinaryBuffer& bb);
};


#include "fab-cc-block.hpp"
