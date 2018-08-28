#pragma once

#include <utility>
#include <numeric>
#include <boost/functional/hash.hpp>

#include "diy/serialization.hpp"
#include "diy/grid.hpp"
#include "diy/link.hpp"
#include "diy/fmt/format.h"
#include "diy/fmt/ostream.h"
#include "diy/point.hpp"
#include <diy/master.hpp>

#include "reeber/amr-vertex.h"
#include "reeber/triplet-merge-tree.h"
#include "reeber/triplet-merge-tree-serialization.h"
#include "reeber/grid.h"
#include "reeber/grid-serialization.h"
#include "reeber/masked-box.h"
#include "reeber/edges.h"

#include "fab-block.h"
#include "reeber/amr_helper.h"

namespace r = reeber;

template<class T, unsigned D>
struct FabTmtBlock;

template<class T, unsigned D>
void set_mask(typename FabTmtBlock<T, D>::MaskedBox& local,
              const diy::GridRef<T, D>& fab,
              const diy::Point<int, D>& v_bounds,
              diy::AMRLink* l,
              const diy::DiscreteBounds& domain,
              const T& rho);


template<class Real, unsigned D>
struct FabTmtBlock
{
    using Shape = diy::Point<int, D>;

    using Grid = r::Grid<Real, D>;
    // index of point: first = index inside box, second = index of a box
    using AmrVertexId = r::AmrVertexId;
    using Value = typename Grid::Value;
    using MaskedBox = r::MaskedBox<D>;
    using Vertex = typename MaskedBox::Position;
    using TripletMergeTree = r::TripletMergeTree<r::AmrVertexId, Value>;
    using Size = std::vector<Real>;
    using AmrVertexContainer = std::vector<AmrVertexId>;

    using Neighbor = typename TripletMergeTree::Neighbor;
    using Node = typename TripletMergeTree::Node;
    using VertexNeighborMap =typename TripletMergeTree::VertexNeighborMap;

    using AmrEdge = r::AmrEdge;
    using AmrEdgeContainer = r::AmrEdgeContainer;
    using AmrEdgeSet = std::set<AmrEdge>;
    using VertexEdgesMap = std::map<AmrVertexId, AmrEdgeContainer>;
    using VertexVertexMap = std::map<AmrVertexId, AmrVertexId>;
    using VertexSizeMap = std::map<AmrVertexId, int>;

    using GidContainer = std::set<int>;
    using GidVector = std::vector<int>;

    template<class Vertex_, class Node>
    struct TmtConnectedComponent
    {
        // types
        using Vertex = Vertex_;
        using Neighbor = Node*;

        // fields
        AmrVertexId root_;
        GidContainer current_neighbors_;
        GidContainer processed_neighbors_;
#ifdef SEND_COMPONENTS
        AmrEdgeContainer outgoing_edges_;
        TripletMergeTree merge_tree_;
#endif

        // methods

        TmtConnectedComponent()
        {
        }

        template<class EC>
        TmtConnectedComponent(const AmrVertexId& root, const EC& _edges) :
                root_(root),
                current_neighbors_({ root.gid })
#ifdef SEND_COMPONENTS
                , outgoing_edges_(_edges.cbegin(), _edges.cend())
#endif
        {
#ifdef SEND_COMPONENTS
            bool debug = false;

            std::transform(outgoing_edges_.begin(), outgoing_edges_.end(),
                           std::inserter(current_neighbors_, current_neighbors_.begin()),
                           [this](const AmrEdge& e) {
                               assert(std::get<0>(e).gid == this->root_.gid);
                               return std::get<1>(e).gid;
                           });

            if (debug)
                fmt::print("entered TmtConnectedComponentConstructor, root = {}, #edges = {}\n", root,
                           outgoing_edges_.size());
#endif
        }

#ifdef SEND_COMPONENTS
        template<class EC>
        void add_edges(const EC& more_edges)
        {
            outgoing_edges_.insert(outgoing_edges_.end(), more_edges.cbegin(), more_edges.cend());

            std::transform(more_edges.cbegin(), more_edges.cend(),
                           std::inserter(current_neighbors_, current_neighbors_.begin()),
                           [this](const AmrEdge& e) {
                               assert(std::get<0>(e).gid == this->root_.gid);
                               return std::get<1>(e).gid;
                           });
        }
#endif
        int is_done() const
        {
            return current_neighbors_ == processed_neighbors_;
        }

    };

    using Component = TmtConnectedComponent<reeber::AmrVertexId, Node>;


    // data

    int gid;
    MaskedBox local_;
    TripletMergeTree mt_;
    TripletMergeTree original_tree_;

    // this vector is not serialized, because we send trees component-wise
    std::vector<Component> components_;

    diy::DiscreteBounds domain_;

    int done_ { 0 };

    //    // will be changed in each communication round
    //    // only for baseiline algorithm
    //    AmrEdgeSet outgoing_edges_;

    // is pre-computed once. does not change
    // only for baseiline algorithm
    AmrEdgeContainer initial_edges_;

    std::map<int, AmrEdgeContainer> gid_to_outgoing_edges_;

    std::set<int> new_receivers_;
    std::set<int> processed_receiveres_;

    GidVector original_link_gids_;

    bool negate_;

    // to store information about local connected component in a serializable way
    VertexVertexMap vertex_to_deepest_;

    // tracking how connected components merge - disjoint sets data structure
    VertexVertexMap components_disjoint_set_parent_;
    VertexSizeMap components_disjoint_set_size_;

    int round_ { 0 };

    // methods

    // simple getters/setters
    const diy::DiscreteBounds& domain() const
    { return domain_; }

    int refinement() const
    { return local_.refinement(); }

    const GidVector& get_original_link_gids() const
    { return original_link_gids_; }


    FabTmtBlock(const diy::GridRef<Real, D>& fab_grid,
                int _ref,
                int _level,
                const diy::DiscreteBounds& _domain,
                const diy::DiscreteBounds& bounds,
                const diy::DiscreteBounds& core,
                int _gid,
                diy::AMRLink* amr_link,
                Real rho,                                           // threshold for LOW value
                bool _negate) :
            gid(_gid),
            local_(project_point<D>(core.min), project_point<D>(core.max), project_point<D>(bounds.min),
                   project_point<D>(bounds.max), _ref, _level, gid, fab_grid.c_order()),
            mt_(_negate),
            original_tree_(_negate),
            domain_(_domain),
            processed_receiveres_({ gid }),
            negate_(_negate)
    {
        bool debug = true;

        std::string debug_prefix = "FabTmtBlock ctor, gid = " + std::to_string(gid);

        if (debug) fmt::print("{} setting mask\n", debug_prefix);

        diy::for_each(local_.mask_shape(), [this, amr_link, &fab_grid, rho](const Vertex& v) {
            set_mask<Real, D>(this->local_, fab_grid, v, amr_link, this->domain(), rho);
        });

        //        if (debug) fmt::print("gid = {}, checking mask\n", gid);
        int max_gid = 0;
        for (int i = 0; i < amr_link->size(); ++i) {
            max_gid = std::max(max_gid, amr_link->target(i).gid);
        }

        local_.check_mask_validity(max_gid);

        reeber::compute_merge_tree2(mt_, local_, fab_grid);

        if (debug) fmt::print("{} local tree computed\n", debug_prefix);

        mt_.make_deep_copy(original_tree_);

        if (debug) fmt::print("{} local tree copied\n", debug_prefix);

        VertexEdgesMap vertex_to_outgoing_edges;

        compute_outgoing_edges(amr_link, vertex_to_outgoing_edges);

        if (debug) fmt::print("{} outgoing edges computed\n", debug_prefix);

        compute_connected_components(vertex_to_outgoing_edges);

        if (debug) fmt::print("{} connected components computed\n", debug_prefix);

        for (int i = 0; i < amr_link->size(); ++i) {
            if (amr_link->target(i).gid != gid) {
                new_receivers_.insert(amr_link->target(i).gid);
                original_link_gids_.push_back(amr_link->target(i).gid);
            }
        }

        if (debug)
            fmt::print( "{}, constructed, refinement = {}, level = {}, local = {}, domain.max = {}, #components = {}\n", debug_prefix, _ref, _level, local_, domain().max, components_.size());
        if (debug)
            fmt::print("{},  constructed, tree.size = {}, new_receivers.size = {}, new_receivers = {}\n",
                       debug_prefix, mt_.size(), new_receivers_.size(), container_to_string(new_receivers_));

        assert(mt_.size() == original_tree_.size());
    }

    FabTmtBlock() {}

    const TripletMergeTree& get_merge_tree() const { return mt_; }

    // return true, if both edge vertices are in the current neighbourhood
    // no checking of mask is performed, if a vertex is LOW, function will return true.
    // Such an edge must be silently ignored in the merge procedure.
    bool edge_exists(const AmrEdge& e) const;

    // return true, if one of the edge's vertices is inside current neighbourhood
    // and the other is outside
    bool edge_goes_out(const AmrEdge& e) const;

//    void adjust_outgoing_edges();
    // diy stuff

    bool deepest_computed(Neighbor n) const
    { return deepest_computed(n->vertex); }

    bool deepest_computed(const AmrVertexId& v) const
    { return vertex_to_deepest_.find(v) != vertex_to_deepest_.cend(); }

    AmrVertexId deepest(Neighbor n) const
    { return deepest(n->vertex); }

    AmrVertexId deepest(const AmrVertexId& v) const;


    void set_deepest(const AmrVertexId& v, const AmrVertexId& deepest)
    { vertex_to_deepest_[v] = deepest; }

    void create_component(const AmrVertexId& deepest_vertex, const AmrEdgeContainer& edges);

    Component& find_component(const AmrVertexId& deepest_vertex);

    void compute_outgoing_edges(diy::AMRLink* l, VertexEdgesMap& vertex_to_outgoing_edges);

    void compute_connected_components(const VertexEdgesMap& vertex_to_outgoing_edges);

//    void delete_low_edges(int sender_gid, AmrEdgeContainer& edges_from_sender);

    void adjust_original_gids(int sender_gid, FabTmtBlock::GidVector& edges_from_sender);

    // disjoint-sets related methods
    bool are_components_connected(const AmrVertexId& deepest_a, const AmrVertexId& deepest_b);

    bool is_component_connected_to_any_internal(const AmrVertexId& deepest);

    void connect_components(const AmrVertexId& deepest_a, const AmrVertexId& deepest_b);

    void add_component_to_disjoint_sets(const AmrVertexId& deepest_vertex);

    int is_done_simple(const std::vector<FabTmtBlock::AmrVertexId>& vertices_to_check);

#ifdef SEND_COMPONENTS
    int are_all_components_done() const;
#endif
    std::vector<AmrVertexId> get_deepest_vertices() const;

    const AmrEdgeContainer& get_all_outgoing_edges()
    { return initial_edges_; }

    // v must be the deepest vertex in a local connected component
    // cannot be const - path compression!
    AmrVertexId find_component_in_disjoint_sets(AmrVertexId v);

    bool gid_must_be_in_link(int gid) const;

    static void* create()
    {
        return new FabTmtBlock;
    }

    static void destroy(void* b)
    {
        delete static_cast<FabTmtBlock*>(b);
    }

    static void save(const void* b, diy::BinaryBuffer& bb);

    static void load(void* b, diy::BinaryBuffer& bb);
};


namespace diy {

    //    template<class R, unsigned D>
    //    struct Serialization<typename FabTmtBlock<R, D>::Component>
    //    {
    //        using Component = FabTmtBlock<R, D>::Component;
    //
    //        static void save(diy::BinaryBuffer& bb, const Component& c)
    //        {
    //        }
    //
    //        static void load(diy::BinaryBuffer& bb, Component& c)
    //        {
    //        }
    //    };
//
    template<class R, unsigned D>
    struct Serialization<FabTmtBlock<R, D>>
    {
        static void save(diy::BinaryBuffer& bb, const FabTmtBlock<R, D>& b)
        {
            FabTmtBlock<R, D>::save(b, bb);
        }

        static void load(diy::BinaryBuffer& bb, FabTmtBlock<R, D>& b)
        {
            FabTmtBlock<R,D>::load(b, bb);
        }
    };
}

#include "fab-tmt-block.hpp"
