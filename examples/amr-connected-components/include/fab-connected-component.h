#pragma once

#include <utility>
#include <vector>
#include <unordered_set>
#include <numeric>

#include "diy/serialization.hpp"

#include "disjoint-sets.h"

#include "reeber/amr-vertex.h"

//#include "../../amr-merge-tree/include/fab-block.h"
#include "reeber/amr_helper.h"



template<class Real>
struct FabConnectedComponent
{
    // types
    struct VertexValue
    {
        reeber::AmrVertexId vertex;
        Real value;
    };

    using AmrVertexId = reeber::AmrVertexId;
    using AmrVertexSet = std::unordered_set<AmrVertexId>;

    // fields
    bool negate_;
    AmrVertexId global_deepest_;    // will be updated in each communication round
    AmrVertexId original_deepest_;

    Real global_integral_value_;
    Real original_integral_value_;
    Real global_deepest_value_;
    Real original_deepest_value_;

    AmrVertexSet current_neighbors_;
    AmrVertexSet processed_neighbors_;

    // methods

    FabConnectedComponent();

    FabConnectedComponent(bool negate, const AmrVertexId& deepest, Real total_value, Real deepest_value);

    bool cmp(Real x, Real y) const;

    int is_done() const;

    void set_global_deepest(const VertexValue& vv);

    void set_current_neighbors(const AmrVertexSet& new_current_neighbhors);

    int must_send_to_gid(int gid) const;

    void mark_gid_processed(int _gid);

    std::string to_string() const;

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
            diy::save(bb, c.global_integral_value_);
            diy::save(bb, c.original_integral_value_);
            diy::save(bb, c.global_deepest_value_);
            diy::save(bb, c.original_deepest_value_);
            diy::save(bb, c.current_neighbors_);
            diy::save(bb, c.processed_neighbors_);
        }

        static void load(BinaryBuffer& bb, Component& c)
        {
            diy::load(bb, c.global_deepest_);
            diy::load(bb, c.original_deepest_);
            diy::load(bb, c.negate_);
            diy::load(bb, c.global_integral_value_);
            diy::load(bb, c.original_integral_value_);
            diy::load(bb, c.global_deepest_value_);
            diy::load(bb, c.original_deepest_value_);
            diy::load(bb, c.current_neighbors_);
            diy::load(bb, c.processed_neighbors_);
        }
    };


};


#include "fab-connected-component.hpp"
