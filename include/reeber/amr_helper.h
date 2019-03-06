#pragma once

#include <assert.h>
#include <string>
#include <sstream>

#include "format.h"
#include "diy/constants.h"
#include "diy/point.hpp"
#include "diy/grid.hpp"
#include "diy/link.hpp"

#include "reeber/amr-vertex.h"

template<class Cont>
std::string container_to_string(const Cont& v)
{
    std::stringstream ss;
    ss << "[";
    for (const auto& x : v) {
        ss << x << ", ";
    }
    ss << "]";
    return ss.str();
}


inline std::ostream& operator<<(std::ostream& os, diy::BlockID bid)
{
    os << "BlockID(proc = " << bid.proc << ", gid = " << bid.gid << ")";
    return os;
}

inline std::set<diy::BlockID> link_unique(diy::AMRLink *amr_link, int gid)
{
    std::set<diy::BlockID> result;
    for (int i = 0; i < amr_link->size(); ++i)
    {
        if (amr_link->target(i).gid != gid)
        {
            result.insert(amr_link->target(i));
        }
    }
    return result;
}

template<unsigned D, size_t DD=DIY_MAX_DIM>
diy::Point<int, D> point_from_dynamic_point(const diy::DynamicPoint<int, DD>& pt)
{
    static_assert(D<=DD, "cast to point drops dimension");
    diy::Point<int, D> result;
    for(int i = 0; i < D; ++i)
        result[i] = pt[i];
    return result;
}

// take point p from level with get_refinement = point_refinement
// return coarsened point for get_refinement = target_refinement
template<class C, unsigned int D>
inline diy::Point<C, D> coarsen_point(const diy::Point<C, D>& p, int point_refinement, int target_refinement)
{
    if (point_refinement == target_refinement) {
        return p;
    }

    //assert(point_refinement >= target_refinement and point_refinement % target_refinement == 0);

    double factor = (double) target_refinement / (double) point_refinement;
    diy::Point<C, D> result;
    for (unsigned int i = 0; i < D; ++i) {
        result[i] = floor(factor * p[i]);
    }
    return result;
}

// take point p from level with get_refinement ref and domain (not refined; domains is assumed to start from origin)
// return wrapped point
template<class C, unsigned int D>
inline diy::Point<C, D> wrap_point(const diy::Point<C, D>& p, const diy::DiscreteBounds& domain, int ref)
{
    // we assume that domain starts from the origin
    //assert(domain.min == decltype(domain.min)::zero());

    using Point = diy::Point<C, D>;

    Point result = p;
    for (unsigned i = 0; i < D; ++i) {
        int d = ref * (domain.max[i] + 1);
        if (p[i] < 0) {
            result[i] += d;
        } else if (p[i] >= d) {
            result[i] -= d;
        }
    }
    return result;
}

// return the first D coordinates of p
template<unsigned int D, class C>
inline diy::Point<C, D> project_point(const typename diy::Point<C, DIY_MAX_DIM>& p)
{
    diy::Point<C, D> result;
    for(unsigned int i = 0; i < D; ++i)
        result[i] = p[i];
    return result;
}

// return the first D coordinates of p
template<size_t D, class C>
inline diy::DynamicPoint<C, D> project_point(const typename diy::DynamicPoint<C, DIY_MAX_DIM>& p)
{
    diy::DynamicPoint<C, D> result { D };
    for(unsigned int i = 0; i < D; ++i)
        result[i] = p[i];
    return result;
}

// v: point in global coordinates from some level
// v_refinement: refinement of the level of v
// link_idx: index of a block in which v must be contained
// l: AMR link
// return index of v in the box with index link_idx
// v must correspond to unique vertex in the box!
// (i.e., the box cannot be a neigbour from above!)
template<unsigned D>
size_t get_vertex_id(const diy::Point<int, D>& v, int v_refinement, int link_idx, diy::AMRLink* l, bool c_order)
{
    using Position = diy::Point<int, D>;

    Position from = point_from_dynamic_point<D>(l->bounds(link_idx).min);
    Position to = point_from_dynamic_point<D>(l->bounds(link_idx).max);
    // TODO: vector refinement
    int refinement = l->refinement(link_idx)[0];


    bool debug = false;

    diy::GridRef<void*, D> grid(nullptr, to - from + Position::one(), c_order);

    if (debug) {
        fmt::print("enter get_vertex_id v = {}, v_ref = {}\n", v, v_refinement);
    }

    //assert(v_refinement >= refinement);

    // bring point to my level
    Position vv = coarsen_point(v, v_refinement, refinement);

    // bring vv to local coordinates
    vv -= from;
    if (debug) { fmt::print("in local coords vv = {}\n", vv); }
    if (debug) { fmt::print("index = {}\n", grid.index(vv)); }

    //assert(0 <= grid.index(vv) and grid.index(vv) < grid.size());

    return grid.index(vv);
}

// return (in global coordinates on the high level)
// bounding vertices of a box that corresponds to v from lower level
template<unsigned D>
std::tuple<diy::Point<int, D>, diy::Point<int, D>>
refine_vertex(const diy::Point<int, D>& v, int v_refinement, int target_refinement)
{
    using Position = diy::Point<int, D>;

    //assert(v_refinement <= target_refinement);
    //assert(target_refinement % v_refinement == 0);

    int ratio = target_refinement / v_refinement;

    Position from = ratio * v;
    Position to = from + (ratio - 1) * Position::one();

    return std::tie(from, to);
}




/**
 *
 * @tparam D dimension, template parameter
 * @param i index of neighbor in link
 * @param l AMRLink
 * @param v_glob vertex in global coordinates
 * @param v_refinement refinement of v_glob
 * @return true, if the core of i-th block in the link contains v_glob
 * if v_glob corresponds to multiple points on the level of the i-th block,
 * only one of these points is checked.
 */
template<unsigned D>
bool neighbor_contains(size_t i, diy::AMRLink* l, diy::Point<int, D>& v_glob, int v_refinement)
{
    using Position = diy::Point<int, D>;

    Position from = point_from_dynamic_point<D>(l->core(i).min);
    Position to = point_from_dynamic_point<D>(l->core(i).max);
    // TODO: vector refinement
    int refinement = l->refinement(i)[0];

    bool debug = false;

    if (debug) {
        fmt::print("Enter neighbor_contains, v_glob = {}, v_refinemtn = {}, from = {}, to = {}, ref = {}\n", v_glob,
                   v_refinement, from, to, refinement);
    }

    if (v_refinement < refinement) {
        // vertex is below, its coordinates must be multiplied to pull it to level_
        // In fact, we check only one vertex in the whole set of vertices
        // that correspond to v_glob at our level, but it is OK for our purposes
        //assert(refinement % v_refinement == 0);
        int scale = refinement / v_refinement;

        for (size_t i = 0; i < D; ++i)
            if (scale * v_glob[i] < from[i] || scale * v_glob[i] > to[i]) {
                return false;
            }

    } else if (v_refinement > refinement) {
        // vertex is above, out coords must be multiplied
        int scale = v_refinement / refinement;
        //assert(v_refinement % refinement == 0);
        if (debug) { fmt::print("In NeighborBox.contains, scale = {}\n", scale); }

        Position new_to = scale * (to + Position::one()) - Position::one();
        for (size_t i = 0; i < D; ++i)
            if (v_glob[i] < scale * from[i] || v_glob[i] > new_to[i]) {
                if (debug) { fmt::print("In NeighborBox.contains, coord = {}, does NOT contain, return false\n", i); }
                return false;
            }
    } else {
        // v_glob is at the same level as box, just compare coordinates
        for (size_t i = 0; i < D; ++i)
            if (v_glob[i] < from[i] || v_glob[i] > to[i]) {
                return false;
            }
    }

    if (debug) { fmt::print("In neighbor_contains, DOES contain, return true\n"); }
    return true;
}

template<unsigned D>
diy::DiscreteBounds refine_bounds(const diy::DiscreteBounds& bounds_in, int scale)
{
    //assert(scale >= 1);
    diy::DiscreteBounds result;

    result.min = scale * bounds_in.min;
    result.max = scale * (bounds_in.max + diy::DynamicPoint<int, 4>::one(D)) - diy::DynamicPoint<int, 4>::one(D);

    for(int i = D; i < DIY_MAX_DIM; ++i) {
        result.min[i] = 0;
        result.max[i] = 0;
    }

    return result;
}

