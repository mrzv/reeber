#pragma once

#include <utility>
#include <vector>
#include <unordered_set>
#include <numeric>

#include "diy/serialization.hpp"

#include "disjoint-sets.h"

#include "reeber/amr-vertex.h"
#include "reeber/triplet-merge-tree.h"
#include "reeber/edges.h"
#include "reeber/amr_helper.h"



template<class Real_>
class FabConnectedComponent
{
public:
    // types
    using Real = Real_;

    struct VertexValue
    {
        reeber::AmrVertexId vertex;
        Real value;
    };

    using AmrVertexId = reeber::AmrVertexId;
    using AmrEdge = reeber::AmrEdge;
    using AmrEdgeContainer = reeber::AmrEdgeContainer;
    using AmrVertexSet = std::unordered_set<AmrVertexId>;
    using GidSet = std::unordered_set<int>;
    using TripletMergeTree = reeber::TripletMergeTree<AmrVertexId, Real>;

private:
    // fields
    bool negate_;
    AmrVertexId global_deepest_;    // will be updated in each communication round
    AmrVertexId original_deepest_;
    Real global_deepest_value_;

    AmrVertexSet current_neighbors_;
    AmrVertexSet processed_neighbors_;

    GidSet current_gids_;
    GidSet processed_gids_;

    AmrEdgeContainer edges_;

public:
    // methods
    FabConnectedComponent();
    FabConnectedComponent(bool negate, const AmrVertexId& deepest, Real deepest_value);

    // getters
    AmrVertexId global_deepest() const { return global_deepest_; }
    AmrVertexId original_deepest() const { return original_deepest_; }
    Real global_deepest_value() const { return global_deepest_value_; }
    const AmrVertexSet& current_neighbors() const { return current_neighbors_; }
    const AmrVertexSet& processed_neighbors() const { return processed_neighbors_; }
    const GidSet& current_gids() const { return current_gids_; }
    const GidSet& processed_gids() const { return processed_gids_;}
    const AmrEdgeContainer edges() const { return edges_; }

    bool cmp(Real x, Real y) const;

    int is_done() const;

    void set_global_deepest(const VertexValue& vv);

    void add_current_neighbor(const AmrVertexId& new_current_neighbor);
    void set_current_neighbors(const AmrVertexSet& new_current_neighbhors);

    int must_send_to_gid(int gid) const;

    void mark_gid_processed(int _gid);
    void mark_neighbor_processed(AmrVertexId v);

    void add_edge(const AmrEdge& e);

    // access to tree
    TripletMergeTree tree_;

    std::string to_string() const;

    friend diy::Serialization<FabConnectedComponent<Real>>;
};


template<class Real>
std::ostream& operator<<(std::ostream& os, const FabConnectedComponent<Real>& c)
{
    os << c.to_string();
    return os;
}


namespace diy
{
    template<class R>
    struct Serialization<FabConnectedComponent<R>>
    {
        using Component = FabConnectedComponent<R>;

        static void save(BinaryBuffer& bb, const Component& c)
        {
            diy::save(bb, c.global_deepest_);
            diy::save(bb, c.original_deepest_);
            diy::save(bb, c.negate_);
//            diy::save(bb, c.global_integral_value_);
//            diy::save(bb, c.original_integral_value_);
            diy::save(bb, c.global_deepest_value_);
//            diy::save(bb, c.original_deepest_value_);
            diy::save(bb, c.current_neighbors_);
            diy::save(bb, c.processed_neighbors_);
        }

        static void load(BinaryBuffer& bb, Component& c)
        {
            diy::load(bb, c.global_deepest_);
            diy::load(bb, c.original_deepest_);
            diy::load(bb, c.negate_);
//            diy::load(bb, c.global_integral_value_);
//            diy::load(bb, c.original_integral_value_);
            diy::load(bb, c.global_deepest_value_);
//            diy::load(bb, c.original_deepest_value_);
            diy::load(bb, c.current_neighbors_);
            diy::load(bb, c.processed_neighbors_);
        }
    };


};


#include "fab-connected-component.hpp"
