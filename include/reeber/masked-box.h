#ifndef REEBER_MASKED_BOX_H
#define REEBER_MASKED_BOX_H

#include <functional>

#include "amr-vertex.h"
#include "range/filtered.h"
#include "range/transformed.h"
#include "range/utility.h"

#include "diy/vertices.hpp"
#include "format.h"

#include "grid.h"
#include "box.h"
#include "vertices.h"

#include "amr_helper.h"

namespace reeber {

    template<unsigned D>
    class MaskedBox
    {
    public:
        using MaskValue = int;   // type of single cell in mask
        using MaskType = Grid<MaskValue, D>;
        using Position = typename MaskType::Vertex;
        using NewDynamicPoint = diy::DynamicPoint<int, D>;

        // all special mask values must be negative
        // indicates vertex is not active due to the function value
        static constexpr MaskValue LOW = -1;
        // indicates vertex is ghost
        static constexpr MaskValue GHOST = -2;
        // indicates vertex is active, i.e. unmasked and value is higher than threshold
        // (or lower, depends on the negate flag)
        static constexpr MaskValue ACTIVE = -3;
        // indicates vertex was never touched
        static constexpr MaskValue UNINIT = -4;

        class FreudenthalLinkIterator;

        using FreudenthalLinkRange = range::iterator_range<FreudenthalLinkIterator>;

        using VI = VerticesIterator<Position>;

        // Topology interface
        using Vertex = reeber::AmrVertexId;

        MaskedBox() {}

        MaskedBox(const NewDynamicPoint& core_from,
                  const NewDynamicPoint& core_to,
                  const NewDynamicPoint& bounds_from,
                  const NewDynamicPoint& bounds_to,
                  int _ref,
                  int _level,
                  int _gid,
                  bool c_order) :
                core_from_( point_from_dynamic_point<D>(core_from)),
                core_to_(point_from_dynamic_point<D>(core_to)),
                bounds_from_(point_from_dynamic_point<D>(bounds_from)),
                bounds_to_(point_from_dynamic_point<D>(bounds_to)),
                local_box_(nullptr, bounds_to_ - bounds_from_ + Position::one(), c_order),
                mask_from_(core_from_ - Position::one()),
                mask_to_(core_to_ + Position::one()),
                core_shape_(core_to_ - core_from_ + Position::one()),
                bounds_shape_(bounds_to_ - bounds_from_ + Position::one()),
                mask_shape_(mask_to_ - mask_from_ + Position::one()),
                ghost_adjustment_(core_from_ - bounds_from_),
                mask_adjustment_(bounds_from_ - mask_from_),
                mask_(mask_shape_, c_order),
                refinement_(_ref),
                level_(_level),
                gid_(_gid)
        {
            assert(ghost_adjustment_ == bounds_to_ - core_to_);
            diy::for_each(mask_.shape(), [this](const Position& p) { this->set_mask(p, this->UNINIT); });
        }

        /**
         *
         * @param max_gid int: maximal gid
         *
         * for debug, sanity check of mask values
         * ghost vertices cannot be active
         * low vertices cannot be ghost
         * mask value must be between 0 and max_gid
         */
        void check_mask_validity(int max_gid)
        {
            bool debug = false;
            diy::for_each(mask_.shape(), [this, max_gid, debug](const Position& p_mask) {
                MaskValue m = this->mask_(p_mask);
                Position p_core = p_mask - Position::one();

                bool is_valid = (m == ACTIVE and not is_outer(p_mask)) or
                                (is_in_core(p_core) and m == LOW) or
                                (m >= 0 and m <= max_gid and m != this->gid());

                if (not is_valid)
                {
                    fmt::print(
                            "in check validity, p_mask ={}, p_core = {}, m = {}, is_ghost = {}, contains = {}, core_from = {}, core_to = {}\n",
                            p_mask, p_core, pretty_mask_value(m), is_outer(p_mask), is_in_core(p_core), core_from_,
                            core_to_);
                    throw std::runtime_error("bad mask value");
                }

                assert(is_valid);
            });
        }

        static unsigned dimension() { return D; }

        /**
         *
         * @return range of all active vertices as vertex indices w.r.t. bounds
         */
        decltype(auto) vertices() const
        {
            return range::iterator_range<VI>(VI::begin(core_from_, core_to_), VI::end(core_from_, core_to_))
                   | range::filtered(std::bind(&MaskedBox::is_active_global, this, std::placeholders::_1))
                   | range::transformed(std::bind(&MaskedBox::global_position_to_vertex, this, std::placeholders::_1) );
        }

        /**
         *
         * @return range of all active core vertices as points in global coordinates
         */
        decltype(auto) active_global_positions() const
        {
            return range::iterator_range<VI>(VI::begin(core_from_, core_to_), VI::end(core_from_, core_to_))
                   | range::filtered(std::bind(&MaskedBox::is_active_global, this, std::placeholders::_1));
        }

        /**
         *
         * @param v AmrVertexId: index of a cell
         * @return range of all vertex indices of active vertices
         * in the link of v.
         */
        decltype(auto) link(const Vertex& v) const
        {
            return local_position_link_vertices(local_position(v));
        }


        // take global position p_glob
        // return outer positions in global coords
        decltype(auto) outer_edge_link(const Position& p_global) const
        {
            return FreudenthalLinkRange(FreudenthalLinkIterator::begin(p_global),
                                        FreudenthalLinkIterator::end(p_global))
                   | range::filtered(std::bind(&MaskedBox::is_outer_edge_start_glob, this, std::placeholders::_1));
        }


        // for quick-and-dirty fix only
        decltype(auto) inside_link(const Position& p_global) const
        {
            Position p_from = p_global;
            Position p_to = p_global;
            for(size_t i = 0; i < D; ++i)
            {
                p_from[i] -= 1;
                p_to[i] += 1;
            }

            return range::iterator_range<VI>(VI::begin(p_from, p_to), VI::end(p_from, p_to))
                    | range::filtered(std::bind(&MaskedBox::is_strictly_inside, this, std::placeholders::_1));
        }

        void swap(MaskedBox& other)
        {
            mask_.swap(other.mask_);
            std::swap(core_from_, other.core_from_);
            std::swap(core_to_, other.core_to_);
            std::swap(bounds_from_, other.bounds_from_);
            std::swap(bounds_to_, other.bounds_to_);
        }

        bool operator==(const MaskedBox& other) const
        {
            return std::tie(core_from_, core_to_, bounds_from_, bounds_to_, mask_) ==
                   std::tie(other.core_from_, other.core_to_, other.bounds_from_, other.bounds_to_, mask_);
        }

        template<class C_, class T_>
        friend std::basic_ostream<C_, T_>&
        operator<<(std::basic_ostream<C_, T_>& out, const MaskedBox& b)
        {
            out << "MaskedBox: " << b.core_from_ << " - " << b.core_to_ << ", mask_shape: " << b.mask_shape();
            out << ", bounds_from_: " << b.bounds_from_ << ", bounds_to_: " << b.bounds_to_;
            out << ", mask_from_: " << b.mask_from_ << ", mask_to_: " << b.mask_to_;
            out << ", ghost_adjustement = " << b.ghost_adjustment_ << ", mask_adjustment = " << b.mask_adjustment_;
            return out;
        }

        int gid() const { return gid_; }
        /**
         *
         * @param p_mask cell in mask coordinates (w.r.t. bounds_from)
         * @param value new mask value
         */
        void set_mask(const Position& p_mask, MaskValue value)
        {
            mask_(p_mask) = value;
        }

        /**
         *
         * @param p_bounds cell in mask coordinates (w.r.t bounds_from)
         * @return mask value
         */
        MaskValue mask(const Position& p_mask) const
        {
#ifdef REEBER_ENABLE_CHECKS
            if (not is_valid_mask_position(p_mask))
            {
                fmt::print("error in mask, p_mask = {}\n", p_mask);
                throw std::runtime_error("bad p_mask");
            }
#endif
            //assert(is_valid_mask_position(p_bounds));
            return mask_(p_mask);
        }

        const MaskType& mask_grid() const
        {
            return mask_;
        }


        /**
         *
         * @param p_mask cell in mask coordinates
         * @return true, if cell is ghost
         * just checks that cell is not in the core
         */
        bool is_outer(const Position& p_mask) const
        {
            for (size_t i = 0; i < D; ++i) {
                if (p_mask[i] < 1 or p_mask[i] > core_shape_[i]) {
                    return true;
                }
            }
            return false;
        }

        Position bounds_from() const { return bounds_from_; }

        Position mask_from() const { return mask_from_; }

        Position mask_shape() const { return mask_shape_; }

        Position core_shape() const { return core_shape_; }

        int level() const { return level_; }

        int refinement() const { return refinement_; }

        std::string pretty_mask_value(MaskValue a) const
        {
            switch (a) {
                case ACTIVE :
                    return "ACTIVE";
                case GHOST :
                    return "GHOST";
                case LOW :
                    return "LOW";
                case UNINIT:
                    return "UNINIT";
                default:
                    return std::to_string(a);
            };
        }

        // get readable mask value, Position must be in local coordinates w.r.t. bounds
        // (i.e, to get mask at bounds_from supply (0,0,0), not (-1,-1,-1)
        std::string pretty_mask_value(const Position& p_mask) const
        {
            return pretty_mask_value(mask_(p_mask));
        }

        // get readable mask by index w.r.t to core (i,e, does not support ghost vertices)
        // do not rename to pretty_mask_value to avoid confusion: Vertex in non-AMR setting is size_t,
        // MaskValue is int, and the caller would have to be very careful
        std::string pretty_mask_value_by_index(const Vertex& v) const
        {
            return pretty_mask_value(mask_by_index(v));
        }

        MaskValue mask_by_index(const Vertex& v) const
        {
            Position p_mask = mask_position(v);

            //assert(Position::zero().is_less_or_eq(p_mask) && mask_shape().is_greater_or_eq(p_mask));

            return mask_(p_mask);
        }

        Position bounds_shape() const { return bounds_shape_; }
        Position ghost_adjustment() const { return ghost_adjustment_; }

        /**
         *
         * @param p_local cell in local coordinates
         * @return cell in global coordinates
         */

        Position global_position_from_local(const Position& p_local) const
        {
            return p_local + bounds_from_;
        }

        /**
         *
         * @param v index of cell (AmrVertexId)
         * @return cell in global coordinates
         */
        Position global_position(const Vertex& v) const
        {
            Position p = local_box_.vertex(static_cast<size_t>(v));
            return global_position_from_local(p);
        }

        /**
         *
         * @param p_global cell in global coordinates
         * @return cell index (AmrVertexId)
         */

        Vertex get_vertex_from_global_position(Position p_global) const
        {
            //assert(bounds_contains_global(p_global));
            Position p_local = local_position_from_global(p_global);
            return AmrVertexId { gid(), local_box_.index(p_local) };
        }

        Position mask_position_from_local(const Position& p_bounds) const
        {
            return p_bounds + mask_adjustment_;
        }

        Position local_position_from_mask(const Position& p_mask) const
        {
            return p_mask - mask_adjustment_;
        }

        /**
         *
         * @param p_global cell in global coordinates
         * @return cell in local coordinates
         */
        Position local_position_from_global(const Position& p_global) const
        {
            return p_global - bounds_from_;
        }

        Position local_position_from_global(const NewDynamicPoint& p_global) const { return local_position_from_global(point_from_dynamic_point<D>(p_global)); }

        /**
         *
         * @param v index of vertex (w.r.t. bounds)
         * @return
         */
        Position mask_position(Vertex v) const
        {
            Position p_bounds = local_box_.vertex(static_cast<size_t>(v));
            return p_bounds + mask_adjustment_;
        }

        /**
         *
         * @param v index of cell (AmrVertexId)
         * @return cell in local coordinates
         */
        Position local_position(const Vertex& v) const
        {
            //assert(static_cast<size_t>(v) < mask_.size());
            auto result = local_box_.vertex(static_cast<size_t>(v));
            return result;
        }

        /**
         *
         * @param p_global cell in global coordinates
         * @return true, if cell belongs to the core of the box (between core_from and core_to)
         */
        bool core_contains_global(const Position& p_global) const;
        bool core_contains_global(const NewDynamicPoint& p_global) const { return core_contains_global(point_from_dynamic_point<D>(p_global)); }

        bool is_on_boundary(const Position& p_global) const;
        bool is_on_boundary(const Vertex& v) const { return is_on_boundary(global_position(v)); }

        bool is_strictly_inside(const Position& p_global) const;
        /**
         *
         * @param p_global cell in global coordinates
         * @return true, if cell belongs to the box with ghosts (between bounds_from and bounds_to)
         */
        bool bounds_contains_global(const Position& p_global) const;
        bool bounds_contains_global(const NewDynamicPoint& p_global) const { return bounds_contains_global(point_from_dynamic_point<D>(p_global)); }


        /**
         *
         * @param p_core cell in core coordinates
         * @return true, if cell belongs to core
         */
        bool is_in_core(const Position& p_core) const;

        /**
         *
         * @param p Cell in local coordinates
         * @return range of all active local vertices of Freudenthal link of p
         * in local coordinates
         */
        decltype(auto) local_position_link(const Position& p_local) const
        {
            return FreudenthalLinkRange(FreudenthalLinkIterator::begin(p_local), FreudenthalLinkIterator::end(p_local))
                   | range::filtered(std::bind(&MaskedBox::is_active_local, this, std::placeholders::_1));
        }

        // for test only
        decltype(auto) local_position_link(const NewDynamicPoint& p) const { return local_position_link(point_from_dynamic_point<D>(p)); }

        Vertex local_position_to_vertex(const Position& p_local) const
        {
            return AmrVertexId { gid(), local_box_.index(p_local) };
        }


        static void save(const void* mb, diy::BinaryBuffer& bb);

        static void load(void* mb, diy::BinaryBuffer& bb);

    private:
        /**
         *
         * @param p cell in local coordinates
         * @return range of all vertex indices of active vertices
         * in the link of p.
         */
        decltype(auto) local_position_link_vertices(const Position& p) const
        {
            return local_position_link(p)
                   | range::transformed(std::bind(&MaskedBox::local_position_to_vertex, this, std::placeholders::_1) );
        }



        Position mask_position_from_global(const Position& p_global) const
        {
            return p_global - bounds_from_ + mask_adjustment_;
        }

        /**
         *
         * @param p_global cell in global coordinates
         * @return true, if cell is active
         */
        bool is_active_global(const Position& p_global) const
        {
            return mask_(mask_position_from_global(p_global)) == ACTIVE;
        }

        /**
         *
         * @param p_local cell in local coordinates (w.r.t bounds_from)
         * @return true, if cell is active
         */
        bool is_active_local(const Position& p_local) const
        {
            Position p_mask = mask_position_from_local(p_local);
            if (not is_valid_mask_position(p_mask))
                return false;
            else
                return mask(p_mask) == ACTIVE;
        }

        /**
         *
         * @param v index of a cell (AmrVertexId)
         * @return true, if cell is active
         */
        bool is_active_index(const Vertex& v) const
        {
            return mask_(mask_position(v)) == ACTIVE;
        }

        /**
         *
         * @param p_global cell in global coordinates
         * @return true, if TODO
         */
        bool is_outer_edge_start_glob(const Position& p_global) const
        {
            Position p_mask = mask_position_from_global(p_global);

            if (not is_valid_mask_position(p_mask))
                return false;

#ifdef REEBER_ENABLE_CHECKS
            if ((mask_(p_mask) == ACTIVE or mask_(p_mask) == LOW) and is_outer(p_mask))
            {
                fmt::print("Error in is_outer_edge_start, p_global = {}, p_mask = {}, this = {}, is outer\n", p_global, p_mask, *this);
                throw std::runtime_error("Error in is_outer_edge_start-2");
            }
#endif

            return mask_(p_mask) != ACTIVE and mask_(p_mask) != LOW;
        }

        bool is_valid_mask_position(const Position& p_mask) const
        {
            for(int i = 0; i < (int)D; ++i)
            {
                if (p_mask[i] > mask_shape_[i] or
                    p_mask[i] < 0)
                {
                    return false;
                }
            }
            return true;
//            return p_mask.is_less_or_eq(mask_shape_) and p_mask.is_greater_or_eq(Position::zero());
        }

        bool is_valid_local_position(const Position& p_bounds) const
        {
            for(int i = 0; i < (int)D; ++i)
            {
                if (p_bounds[i] > bounds_shape_[i] or
                    p_bounds[i] < 0)
                {
                    return false;
                }
            }
            return true;
        }

        Vertex global_position_to_vertex(const Position& p_global) const
        {
            return local_position_to_vertex(local_position_from_global(p_global));
        }


        // data
        const Position core_from_, core_to_;
        const Position bounds_from_, bounds_to_;
        diy::GridRef<void*, D> local_box_ { nullptr, Position::zero() };
        const Position mask_from_, mask_to_;
        const Position core_shape_;
        const Position bounds_shape_;
        const Position mask_shape_;
        const Position ghost_adjustment_;
        const Position mask_adjustment_;
        MaskType mask_;
        const int refinement_ { 0 };
        const int level_ { -1 };
        const int gid_ { -1 };
    };

}

#include "masked-box.hpp"

#endif
