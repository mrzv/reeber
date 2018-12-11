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
struct FabComponentBlock
{
    using Shape = diy::Point<int, D>;

    using Grid = r::Grid<Real, D>;
    using GridRef = r::GridRef<Real, D>;
    // index of point: first = index inside box, second = index of a box
    using AmrVertexId = r::AmrVertexId;

    using Component = FabConnectedComponent<Real>;

    using VertexValue = typename Component::VertexValue;

    using Value = typename Grid::Value;
    using MaskedBox = r::MaskedBox<D>;
    using Vertex = typename MaskedBox::Position;
    using AmrVertexContainer = std::vector<AmrVertexId>;
    using AmrVertexSet = typename Component::AmrVertexSet;

    using AmrEdge = r::AmrEdge;
    using AmrEdgeContainer = r::AmrEdgeContainer;
    using AmrEdgeSet = std::set<AmrEdge>;
    using VertexEdgesMap = std::map<AmrVertexId, AmrEdgeContainer>;

    using GidContainer = std::set<int>;
    using GidVector = std::vector<int>;

    using RealType = Real;

    using UnionFind = DisjointSets<AmrVertexId>;
    using VertexVertexMap = UnionFind::VertexVertexMap;
    using VertexSizeMap = UnionFind::VertexSizeMap;


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

    UnionFind disjoint_sets_;   // keep topology of graph of connected components
    std::vector<Component> components_;

    VertexVertexMap vertex_to_deepest_;

    diy::DiscreteBounds domain_;

    int done_{0};

    //    // will be changed in each communication round
    //    // only for baseiline algorithm
    //    AmrEdgeSet outgoing_edges_;

    // is pre-computed once. does not change
    // only for baseiline algorithm
//    AmrEdgeContainer initial_edges_;

    std::unordered_map<int, AmrEdgeContainer> gid_to_outgoing_edges_;

    bool negate_;

    // tracking how connected components merge - disjoint sets data structure

    int round_{0};

    // methods

    // simple getters/setters
    const diy::DiscreteBounds& domain() const { return domain_; }

    int refinement() const { return local_.refinement(); }

    int level() const { return local_.level(); }

    FabComponentBlock(diy::GridRef<Real, D>& fab_grid,
                      int _ref,
                      int _level,
                      const diy::DiscreteBounds& _domain,
                      const diy::DiscreteBounds& bounds,
                      const diy::DiscreteBounds& core,
                      int _gid,
                      diy::AMRLink* amr_link,
                      Real rho,                                           // threshold for LOW value
                      bool _negate,
                      bool is_absolute_threshold) :
            gid(_gid),
            local_(project_point<D>(core.min), project_point<D>(core.max), project_point<D>(bounds.min),
                   project_point<D>(bounds.max), _ref, _level, gid, fab_grid.c_order()),
            fab_(fab_grid.data(), fab_grid.shape(), fab_grid.c_order()),
            domain_(_domain),
            negate_(_negate)
    {
        bool debug = false;

        std::string debug_prefix = "FabComponentBlock ctor, gid = " + std::to_string(gid);

        if(debug) fmt::print("{} setting mask\n", debug_prefix);

        diy::for_each(local_.mask_shape(), [this, amr_link, rho, is_absolute_threshold](const Vertex& v) {
            this->set_mask(v, amr_link, rho, is_absolute_threshold);
        });

        //        if (debug) fmt::print("gid = {}, checking mask\n", gid);
        int max_gid = 0;
        for(int i = 0; i < amr_link->size(); ++i)
        {
            max_gid = std::max(max_gid, amr_link->target(i).gid);
        }

        //local_.check_mask_validity(max_gid);

        if(is_absolute_threshold)
        {
            init(rho, amr_link);
        }
    }

    FabComponentBlock() :
            fab_(nullptr, diy::Point<int, D>::zero())
    {}

    void init(Real absolute_rho, diy::AMRLink* amr_link);

    // compare w.r.t negate_ flag
    bool cmp(Real a, Real b) const;

    void set_low(const diy::Point<int, D>& v_bounds,
                 const Real& absolute_rho);

    void set_mask(const diy::Point<int, D>& v_bounds,
                  diy::AMRLink* l,
                  const Real& rho,
                  bool is_absolute_threshold);

    // return true, if both edge vertices are in the current neighbourhood
    // no checking of mask is performed, if a vertex is LOW, function will return true.
    // Such an edge must be silently ignored in the merge procedure.
    bool edge_exists(const AmrEdge& e) const;

    // return true, if one of the edge's vertices is inside current neighbourhood
    // and the other is outside
    bool edge_goes_out(const AmrEdge& e) const;

    void compute_outgoing_edges(diy::AMRLink* l, VertexEdgesMap& vertex_to_outgoing_edges);

    void compute_original_connected_components(const VertexEdgesMap& vertex_to_outgoing_edges);

//    void compute_final_connected_components();
//
    void delete_low_edges(int sender_gid, AmrEdgeContainer& edges_from_sender, const VertexVertexMap& received_vertex_to_deepest);
//
    void adjust_outgoing_edges();

    bool is_component_connected_to_any_internal(const AmrVertexId& deepest);

    void sparsify_prune_original_tree() {}

//    void add_received_original_vertices(const VertexVertexMap& received_vertex_to_deepest);

    int get_n_components_for_gid(int gid) const;

    int are_all_components_done() const;

    std::vector<AmrVertexId> get_current_deepest_vertices() const;

//    int n_undone_components() const;

    int is_done_simple(const std::vector<FabComponentBlock::AmrVertexId>& vertices_to_check);

    void compute_local_integral(Real rho, Real theta);

    bool check_symmetry(int gid, const std::vector<Component>& received_components);

    Real scaling_factor() const;

    Component& get_component_by_deepest(const AmrVertexId& deepest)
    {
        auto res_iter = std::find_if(components_.begin(), components_.end(), [deepest](const Component& c) { return c.original_deepest_ == deepest; });
        if (res_iter == components_.end())
            throw std::runtime_error("error in find_componenent, bad deepest");
        return  *res_iter;
    }

    std::vector<AmrVertexId> get_original_deepest_vertices() const;

//    const AmrEdgeContainer& get_all_outgoing_edges() { return initial_edges_; }

//    // v must be the deepest vertex in a local connected component
//    // cannot be const - path compression!
//    AmrVertexId find_component_in_disjoint_sets(AmrVertexId v);

    static void* create()
    {
        return new FabComponentBlock;
    }

    static void destroy(void* b)
    {
        delete static_cast<FabComponentBlock*>(b);
    }

    static void save(const void* b, diy::BinaryBuffer& bb);

    static void load(void* b, diy::BinaryBuffer& bb);
};


#include "fab-cc-block.hpp"
