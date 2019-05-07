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

#include "small_set.h"


template<class Real_>
class FabConnectedComponent {
public:
    // types
    using Real = Real_;

    struct VertexValue {
        reeber::AmrVertexId vertex;
        Real value;
    };

    using AmrVertexId = reeber::AmrVertexId;
    using AmrEdge = reeber::AmrEdge;
    using AmrEdgeContainer = reeber::AmrEdgeContainer;
    using AmrVertexSet = std::unordered_set<AmrVertexId>;
//    using AmrVertexSet = SmallSet<AmrVertexId>;
    using GidSet = std::unordered_set<int>;
    using TripletMergeTree = reeber::TripletMergeTree<AmrVertexId, Real>;
    using VertexValueMap = std::unordered_map<AmrVertexId, Real>;
    using UnionFind = DisjointSets<AmrVertexId>;

    using ExtraValues = std::map<std::string, Real>;

private:
    // fields
    bool negate_;
    AmrVertexId original_deepest_;
    Real global_deepest_value_;

    AmrVertexSet current_neighbors_;

    GidSet current_gids_;
    GidSet processed_gids_;

    AmrEdgeContainer edges_;

    std::size_t n_prev_current_neighbors_ { 0 };

    ExtraValues extra_integral_values_;

public:
    // methods
    // ctors
    FabConnectedComponent();
    FabConnectedComponent(bool negate, const AmrVertexId& deepest, Real deepest_value, const ExtraValues& extra_integral_values);

    // getters
    AmrVertexId         original_deepest()  const { return original_deepest_; }
    const AmrVertexSet& current_neighbors() const { return current_neighbors_; }
    const GidSet&       current_gids()      const { return current_gids_; }
    const GidSet&       processed_gids()    const { return processed_gids_; }
    const ExtraValues&  extra_values()      const { return extra_integral_values_; }

    const AmrEdgeContainer edges()          const { return edges_; }

    bool cmp(Real x, Real y) const;

    bool is_done_sending() const;

    void add_current_neighbor(const AmrVertexId& new_current_neighbor);
    void set_current_neighbors(const AmrVertexSet& new_current_neighbors);

    bool must_send_to_gid(int gid) const;
    bool must_send_tree_to_gid(int gid) const;

    int must_send_neighbors() const;

    void mark_all_gids_processed();

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


namespace diy {
    template<class R>
    struct Serialization<FabConnectedComponent<R>> {
        using Component = FabConnectedComponent<R>;

        static void save(BinaryBuffer& bb, const Component& c)
        {
            diy::save(bb, c.original_deepest_);
            diy::save(bb, c.negate_);
            diy::save(bb, c.global_deepest_value_);
            diy::save(bb, c.current_neighbors_);
            diy::save(bb, c.current_gids_);
            diy::save(bb, c.processed_gids_);
        }

        static void load(BinaryBuffer& bb, Component& c)
        {
            diy::load(bb, c.original_deepest_);
            diy::load(bb, c.negate_);
            diy::load(bb, c.global_deepest_value_);
            diy::load(bb, c.current_neighbors_);
            diy::load(bb, c.current_gids_);
            diy::load(bb, c.processed_gids_);
        }
    };
}

#include "fab-connected-component.hpp"
