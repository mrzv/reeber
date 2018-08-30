#ifndef REEBER_MASKED_BOX_H
#define REEBER_MASKED_BOX_H

#include <functional>

#include "amr-vertex.h"
#include "range/filtered.h"
#include "range/transformed.h"
#include "range/utility.h"

#include "diy/vertices.hpp"
#include "diy/fmt/format.h"

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

        MaskedBox() :
                core_from_(),
                core_to_(),
                bounds_from_(),
                bounds_to_(),
                core_shape_(),
                mask_shape_(),
                ghost_adjustment_(),
                mask_(),
                refinement_(0),
                level_(0),
                gid_(-1)
        {}

        MaskedBox(const Position& core_from,
                  const Position& core_to,
                  const Position& bounds_from,
                  const Position& bounds_to,
                  int _ref,
                  int _level,
                  int _gid,
                  bool c_order) :
                core_from_(core_from),
                core_to_(core_to),
                bounds_from_(bounds_from),
                bounds_to_(bounds_to),
                core_shape_(core_to - core_from + Position::one()),
                mask_shape_(bounds_to - bounds_from + Position::one()),
                ghost_adjustment_(core_from - bounds_from),
                mask_(mask_shape_, c_order),
                refinement_(_ref),
                level_(_level),
                gid_(_gid)
        {
            //assert(ghost_adjustment_ == bounds_to - core_to);
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
            diy::for_each(mask_.shape(), [this, max_gid, debug](const Position& p_bounds) {
                MaskValue m = this->mask_(p_bounds);
                Position p_core = p_bounds - Position::one();

                if (debug) fmt::print( "in check validity, p_bounds ={}, p_core = {}, m = {}, is_ghost = {}, contains = {}, core_from = {}, core_to = {}\n", p_bounds, p_core, pretty_mask_value(m), is_ghost(p_core), contains_local(p_bounds), core_from_, core_to_);

                assert((m == ACTIVE and not is_ghost(p_bounds)) or
                       (contains_local(p_bounds) and m == LOW) or
                       (m >= 0 and m <= max_gid and m != this->gid()));
            });
        }

        static unsigned dimension() { return D; }

        /**
         *
         * @return range of all active vertices as vertex indices
         */
        decltype(auto) vertices() const
        {
            return range::iterator_range<VI>(VI::begin(bounds_from_, bounds_to_), VI::end(bounds_from_, bounds_to_))
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
         * @param p cell in local coordinates
         * @return range of all vertex indices of active vertices
         * in the link of p.
         */
        decltype(auto) local_position_link_vertices(const Position& p) const
        {
            return local_position_link(p)
                   | range::transformed(std::bind(&MaskedBox::local_position_to_vertex, this, std::placeholders::_1) );
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

        /**
         *
         * @param p Cell in local coordinates
         * @return range of all active local vertices of Freudenthal link of p
         * in local coordinates
         */
        decltype(auto) local_position_link(const Position& p) const
        {
            return FreudenthalLinkRange(FreudenthalLinkIterator::begin(p), FreudenthalLinkIterator::end(p))
                   | range::filtered(std::bind(&MaskedBox::is_active_local, this, std::placeholders::_1));
        }

        /**
         *
         * @param p_global cell in global coordinates
         * @return true, if cell belongs to the core of the box (between core_from and core_to)
         */
        bool core_contains_global(const Position& p_global) const;

        /**
         *
         * @param p_global cell in global coordinates
         * @return true, if cell belongs to the box with ghosts (between bounds_from and bounds_to)
         */
        bool bounds_contains_global(const Position& p_global) const;


        /**
         *
         * @param p_local cell in local coordinates
         * @return true, if cell belongs to the box
         */
        bool contains_local(const Position& p_local) const;

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
            return out;
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
         * @param p_global cell in global coordinates
         * @return cell index (AmrVertexId)
         */

        Vertex get_vertex_from_global_position(Position p_global) const
        {
            //assert(bounds_contains_global(p_global));
            Position p_local = local_position_from_global(p_global);
            return AmrVertexId { gid(), mask_.index(p_local) };
        }

        /**
         *
         * @param p_global cell in global coordinates
         * @return true, if cell is active
         */
        bool is_active_global(const Position& p_global) const
        {
            return mask_(p_global - bounds_from_) == ACTIVE;
        }

        /**
         *
         * @param p_local cell in local coordinates (w.r.t bounds_from)
         * @return true, if cell is active
         */
        bool is_active_local(const Position& p_local) const
        {
            return mask_(p_local) == ACTIVE;
        }

        /**
         *
         * @param v index of a cell (AmrVertexId)
         * @return true, if cell is active
         */
        bool is_active_index(const Vertex& v) const
        {
            return mask_(local_position(v)) == ACTIVE;
        }

        /**
         *
         * @param p_global cell in global coordinates
         * @return true, if TODO
         */
        bool is_outer_edge_start_glob(const Position& p_global) const
        {
            Position p_bounds = p_global - bounds_from_;
            return mask_(p_bounds) != ACTIVE and mask_(p_bounds) != LOW;
        }

        bool is_valid_mask_position(const Position& p_bounds) const
        {
            return p_bounds.is_less_or_eq(mask_shape_) and p_bounds.is_greater_or_eq(Position::zero());
        }

        /**
         *
         * @param v index of cell (AmrVertexId)
         * @return cell in local coordinates
         */
        Position local_position(const Vertex& v) const
        {
            //assert(static_cast<size_t>(v) < mask_.size());
            auto result = mask_.vertex(static_cast<size_t>(v));
            return result;
            return mask_.vertex(static_cast<size_t>(v));
        }

//        void check_me() const
//        {
//            diy::for_each(mask_.shape(), [this](const Position& p) {
//                assert(this->local_position(this->local_position_to_vertex(p)) == p);
//            });
//            fmt::print("CHECK OK\n");
//        }

        /**
         *
         * @param v index of cell (AmrVertexId)
         * @return cell in global coordinates
         */
         Position global_position(const Vertex& v) const
        {
            Position p = mask_.vertex(static_cast<size_t>(v));
            return p + bounds_from();
        }


        /**
         *
         * @return gid of block to which masked box belongs
         */
        int gid() const { return gid_; }

        /**
         *
         * @param p cell in mask coordinates (w.r.t. bounds_from)
         * @param value new mask value
         */
        void set_mask(const Position& p, MaskValue value)
        {
            mask_(p) = value;
        }

        /**
         *
         * @param p_bounds cell in mask coordinates (w.r.t bounds_from)
         * @return mask value
         */
        MaskValue mask(const Position& p_bounds) const
        {
            //assert(is_valid_mask_position(p_bounds));
            return mask_(p_bounds);
        }

        const MaskType& mask_grid() const
        {
            return mask_;
        }


        /**
         *
         * @param p_local cell in local coordinates
         * @return true, if cell is ghost
         * just checks that cell is not in the core
         */
        bool is_ghost(const Position& p_local)
        {
            Position p_core = p_local - ghost_adjustment_;
            for (size_t i = 0; i < D; ++i) {
                if (p_core[i] < 0 or p_core[i] >= core_shape_[i]) {
                    return true;
                }
            }
            return false;
        }

        Vertex global_position_to_vertex(const Position& p_global) const
        {
            return local_position_to_vertex(local_position_from_global(p_global));
        }

        Vertex local_position_to_vertex(const Position& p_local) const
        {
            return AmrVertexId { gid(), mask_.index(p_local) };
        }

        Position bounds_from() const { return bounds_from_; }

        Position mask_shape() const { return mask_shape_; }

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
        std::string pretty_mask_value(const Position& p_bounds) const
        {
            return pretty_mask_value(mask_(p_bounds));
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
            Position p_bounds = local_position(v);

            //assert(Position::zero().is_less_or_eq(p_bounds) && mask_shape().is_greater_or_eq(p_bounds));

            return mask_(p_bounds);
        }

        Position core_shape() const { return core_shape_; }

        static void save(const void* mb, diy::BinaryBuffer& bb);

        static void load(void* mb, diy::BinaryBuffer& bb);

    private:
        // data
        const Position core_from_, core_to_;
        const Position bounds_from_, bounds_to_;
        const Position core_shape_, mask_shape_;
        const Position ghost_adjustment_;
        MaskType mask_;
        const int refinement_;
        const int level_;
        const int gid_;

    };

}

#include "masked-box.hpp"

#endif
