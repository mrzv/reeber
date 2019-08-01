#include "catch/catch.hpp"

#include <sstream>
#include <iostream>

#include <diy/master.hpp>
#include <diy/io/block.hpp>
#include <opts/opts.h>

#include <reeber/box.h>

#include "fab-block.h"
#include "fab-cc-block.h"
#include "reader-interfaces.h"
#include "diy/vertices.hpp"
#include "reeber/grid.h"

TEST_CASE("Small check", "[masked_box][dim2]")
{
    using MaskedBox = reeber::MaskedBox<2>;
    using Position = MaskedBox::Position;
    using DynPoint = MaskedBox::NewDynamicPoint;
    using AmrVertexId = reeber::AmrVertexId;

    const int refinement = 1;
    const int level = 0;
    const int gid = 0;

    const DynPoint e1{1, 0};
    const DynPoint e2{0, 1};
    const DynPoint one{1, 1};

    const DynPoint core_from{3, 3};
    const DynPoint core_to{4, 4};
    const DynPoint bounds_from = core_from - one;
    const DynPoint bounds_to = core_to + one;

    const DynPoint core_shape = core_to - core_from + one;
    const DynPoint bounds_shape = bounds_to - bounds_from + one;
    const size_t masked_box_size = bounds_shape[0] * bounds_shape[1];
    const int core_size = core_shape[0] * core_shape[1];

    MaskedBox mb(core_from, core_to, bounds_from, bounds_to, refinement, level, gid, false);

    MaskedBox mb_no_ghosts(core_from, core_to, core_from, core_to, refinement, level, gid, false);

    SECTION("check is_active, simple case: all vertices are active")
    {

        // set all inner vertices to be active
        diy::for_each(mb.mask_shape(), [&mb](const Position& p) {
            if (mb.is_outer(p))
                // pretend ghost vertices are masked by another box, does not matter which
                mb.set_mask(p, 2);
            else
                mb.set_mask(p, MaskedBox::ACTIVE);
        });

        auto vs = mb.vertices();
        std::set<AmrVertexId> vertices(std::begin(vs), std::end(vs));

        std::set<AmrVertexId> correct_vertices;

        for(size_t i = 0; i < masked_box_size; ++i)
        {
            if (not mb.is_outer(mb.mask_position(AmrVertexId{gid, i})))
                correct_vertices.emplace(gid, i);
        }

        REQUIRE(vertices.size() == core_size);
        REQUIRE(vertices.size() == correct_vertices.size());
        REQUIRE(vertices == correct_vertices);

        REQUIRE(mb.is_outer(Position{0, 0}));
        REQUIRE(mb.is_outer(Position{0, 1}));
        REQUIRE(mb.is_outer(Position{1, 0}));
        REQUIRE(not mb.is_outer(Position{1, 1}));
        REQUIRE(not mb.is_outer(Position{1, 2}));
        REQUIRE(not mb.is_outer(Position{2, 1}));
        REQUIRE(not mb.is_outer(Position{2, 2}));
        REQUIRE(mb.is_outer(Position{3, 2}));
        REQUIRE(mb.is_outer(Position{2, 3}));
        REQUIRE(mb.is_outer(Position{3, 3}));
        REQUIRE(mb.is_outer(Position{4, 4}));
        REQUIRE(mb.is_outer(Position{4, 5}));
        REQUIRE(mb.is_outer(Position{5, 4}));
        REQUIRE(mb.is_outer(Position{5, 5}));
        REQUIRE(mb.bounds_from() == Position({2, 2}));
        REQUIRE(mb.mask_from() == Position({2, 2}));
        REQUIRE(mb.mask_shape() == Position({4, 4}));
        REQUIRE(mb.bounds_shape() == Position({4, 4}));
        REQUIRE(mb.core_shape() == Position({2, 2}));
        REQUIRE(mb.global_position_from_local(Position({0, 0})) == mb.bounds_from());
        REQUIRE(mb.global_position(AmrVertexId({0, 0})) == mb.bounds_from());
        REQUIRE(mb.get_vertex_from_global_position(mb.bounds_from()) == AmrVertexId({0, 0}));
        REQUIRE(mb.mask_position_from_local(Position({0, 0})) == Position({0, 0}));
        REQUIRE(mb.local_position_from_global(Position({2, 2})) == Position({0, 0}));
        REQUIRE(mb.mask_position(AmrVertexId({0, 0})) == Position({0, 0}));
    }

    SECTION("check no-ghosts case, simple case: all vertices are active")
    {

        // set all inner vertices to be active
        diy::for_each(mb_no_ghosts.mask_shape(), [&mb_no_ghosts](const Position& p) {
            if (mb_no_ghosts.is_outer(p))
                // pretend ghost vertices are masked by another box, does not matter which
                mb_no_ghosts.set_mask(p, 2);
            else
                mb_no_ghosts.set_mask(p, MaskedBox::ACTIVE);
        });

        auto vs = mb_no_ghosts.vertices();
        std::set<AmrVertexId> vertices(std::begin(vs), std::end(vs));

        std::set<AmrVertexId> correct_vertices;

        for(size_t i = 0; i < masked_box_size; ++i)
        {
            if (not mb_no_ghosts.is_outer(mb_no_ghosts.mask_position(AmrVertexId{gid, i})))
                correct_vertices.emplace(gid, i);
        }

        REQUIRE(vertices.size() == core_size);
        REQUIRE(vertices.size() == correct_vertices.size());
        REQUIRE(vertices == correct_vertices);

        REQUIRE(mb_no_ghosts.is_outer(Position{0, 0}));
        REQUIRE(mb_no_ghosts.is_outer(Position{0, 1}));
        REQUIRE(mb_no_ghosts.is_outer(Position{1, 0}));
        REQUIRE(not mb_no_ghosts.is_outer(Position{1, 1}));
        REQUIRE(not mb_no_ghosts.is_outer(Position{1, 2}));
        REQUIRE(not mb_no_ghosts.is_outer(Position{2, 1}));
        REQUIRE(not mb_no_ghosts.is_outer(Position{2, 2}));
        REQUIRE(mb_no_ghosts.is_outer(Position{3, 2}));
        REQUIRE(mb_no_ghosts.is_outer(Position{2, 3}));
        REQUIRE(mb_no_ghosts.is_outer(Position{3, 3}));
        REQUIRE(mb_no_ghosts.is_outer(Position{4, 4}));
        REQUIRE(mb_no_ghosts.is_outer(Position{4, 5}));
        REQUIRE(mb_no_ghosts.is_outer(Position{5, 4}));
        REQUIRE(mb_no_ghosts.is_outer(Position{5, 5}));
        REQUIRE(mb_no_ghosts.bounds_from() == Position({3, 3}));
        REQUIRE(mb_no_ghosts.mask_from() == Position({2, 2}));
        REQUIRE(mb_no_ghosts.mask_shape() == Position({4, 4}));
        REQUIRE(mb_no_ghosts.bounds_shape() == Position({2, 2}));
        REQUIRE(mb_no_ghosts.core_shape() == Position({2, 2}));
        REQUIRE(mb_no_ghosts.global_position_from_local(Position({0, 0})) == mb_no_ghosts.bounds_from());
        REQUIRE(mb_no_ghosts.global_position(AmrVertexId({0, 0})) == mb_no_ghosts.bounds_from());
        REQUIRE(mb_no_ghosts.get_vertex_from_global_position(mb_no_ghosts.bounds_from()) == AmrVertexId({0, 0}));
        REQUIRE(mb_no_ghosts.mask_position_from_local(Position({0, 0})) == Position({1, 1}));
        REQUIRE(mb_no_ghosts.local_position_from_global(Position({2, 2})) == Position({-1, -1}));
        REQUIRE(mb_no_ghosts.mask_position(AmrVertexId({0, 0})) == Position({1, 1}));
    }
}

TEST_CASE("Check masked_box in 2 dimensions", "[masked_box][dim2]")
{
    using MaskedBox = reeber::MaskedBox<2>;
    using Position = MaskedBox::Position;
    using DynPoint = MaskedBox::NewDynamicPoint;
    using AmrVertexId = reeber::AmrVertexId;
    using Vertex = AmrVertexId;

    const DynPoint one{1, 1};
    const DynPoint e1{1, 0};
    const DynPoint e2{0, 1};

    const Position pone{1, 1};
    const Position pe1{1, 0};
    const Position pe2{0, 1};

    const int refinement = 1;
    const int level = 0;
    const int gid = 0;
    const DynPoint core_from{3, 2};
    const DynPoint core_to{10, 15};
    const DynPoint bounds_from = core_from - one;
    const DynPoint bounds_to = core_to + one;

    const DynPoint corner_upper_left = core_from + 13 * e2;
    const DynPoint corner_lower_right = core_from + 7 * e1;

    const DynPoint core_shape = core_to - core_from + one;

    const int core_size = core_shape[0] * core_shape[1];

    MaskedBox mb(core_from, core_to, bounds_from, bounds_to, refinement, level, gid, false);

    const size_t masked_box_size = mb.mask_grid().size();

    SECTION("check is_outer")
    {

        // around  lower left corner
        REQUIRE(mb.is_outer(Position{0, 0}));
        REQUIRE(mb.is_outer(Position{1, 0}));
        REQUIRE(mb.is_outer(Position{0, 1}));
        REQUIRE(mb.is_outer(Position{0, 3}));
        REQUIRE(mb.is_outer(Position{3, 0}));

        REQUIRE(not mb.is_outer(Position{1, 1}));
        REQUIRE(not mb.is_outer(Position{2, 3}));
        REQUIRE(not mb.is_outer(Position{1, 1}));


        // around upper right corner
        REQUIRE(not mb.is_outer(Position{8, 14}));
        REQUIRE(not mb.is_outer(Position{8, 13}));
        REQUIRE(not mb.is_outer(Position{7, 14}));
        REQUIRE(not mb.is_outer(Position{7, 13}));

        REQUIRE(mb.is_outer(Position{9, 15}));
        REQUIRE(mb.is_outer(Position{8, 15}));
        REQUIRE(mb.is_outer(Position{9, 14}));
        REQUIRE(mb.is_outer(Position{7, 15}));
        REQUIRE(mb.is_outer(Position{9, 13}));
    }

    SECTION("check local_position and global position")
    {
        REQUIRE(mb.local_position(AmrVertexId{gid, 0}) == Position{0, 0});
        REQUIRE(mb.local_position(AmrVertexId{gid, 1}) == Position{1, 0});

        REQUIRE(mb.global_position(AmrVertexId{gid, 0}) == point_from_dynamic_point<2>(bounds_from) + Position{0, 0});
        REQUIRE(mb.global_position(AmrVertexId{gid, 1}) == point_from_dynamic_point<2>(bounds_from) + Position{1, 0});
    }

    SECTION("check is_active, simple case: all vertices are active")
    {

        // set all inner vertices to be active
        diy::for_each(mb.mask_shape(), [&mb](const Position& p) {
            if (mb.is_outer(p))
                // pretend ghost vertices are masked by another box, does not matter which
                mb.set_mask(p, 2);
            else
                mb.set_mask(p, MaskedBox::ACTIVE);
        });

        auto vs = mb.vertices();
        std::set<Vertex> vertices(std::begin(vs), std::end(vs));

        std::set<Vertex> correct_vertices;
        for(size_t i = 0; i < masked_box_size; ++i)
        {
            if (not mb.is_outer(mb.local_position(AmrVertexId{gid, i})))
                correct_vertices.emplace(gid, i);
        }

        REQUIRE(vertices.size() == core_size);
        REQUIRE(vertices.size() == correct_vertices.size());
        REQUIRE(vertices == correct_vertices);

        auto aps = mb.active_global_positions();
        std::set<Position> active_global_positions_1(std::begin(aps), std::end(aps));
        std::set<DynPoint> active_global_positions;

        for(const Position& p : active_global_positions_1)
        {
            DynPoint pp{p[0], p[1]};
            active_global_positions.insert(pp);
        }

        REQUIRE((active_global_positions.count(core_from) == 1 and
                active_global_positions.count(core_to) == 1 and
                active_global_positions.count(corner_lower_right) == 1 and
                active_global_positions.count(corner_upper_left) == 1));

        REQUIRE((active_global_positions.count(core_from + e1) == 1 and
                active_global_positions.count(core_to - e2) == 1 and
                active_global_positions.count(corner_lower_right + e2) == 1 and
                active_global_positions.count(corner_upper_left + e1) == 1));

        REQUIRE((active_global_positions.count(core_from - e1) == 0 and
                active_global_positions.count(core_to + e2) == 0 and
                active_global_positions.count(corner_lower_right - e2) == 0 and
                active_global_positions.count(corner_upper_left - e1) == 0));

        auto inner_vertex = e1 + e2 + e1 + e2;

        auto ivl = mb.local_position_link(inner_vertex);
        std::set<Position> inner_vertex_link(std::begin(ivl), std::end(ivl));

        std::set<DynPoint> correct_inner_vertex_link_1{inner_vertex + e1,
                                                       inner_vertex + e2,
                                                       inner_vertex - e1,
                                                       inner_vertex - e2,
                                                       inner_vertex + e1 + e2,
                                                       inner_vertex - e1 - e2};

        std::set<Position> correct_inner_vertex_link;
        for(const DynPoint& p : correct_inner_vertex_link_1)
        {
            Position pp{p[0], p[1]};
            correct_inner_vertex_link.insert(pp);
        }

        REQUIRE(inner_vertex_link == correct_inner_vertex_link);

        auto cfl = mb.local_position_link(Position::one());
        std::set<Position> core_from_link(std::begin(cfl), std::end(cfl));

        std::set<DynPoint> correct_core_from_link_1{one + e1, one + e2, one + e1 + e2};
        std::set<Position> correct_core_from_link;
        for(auto p : correct_core_from_link_1)
        {
            Position pp{p[0], p[1]};
            correct_core_from_link.insert(pp);
        }

        REQUIRE(core_from_link == correct_core_from_link);

        auto ull = mb.local_position_link(corner_upper_left - bounds_from);
        std::set<Position> upper_left_corner_link(std::begin(ull), std::end(ull));

        Position corner_upper_left_local = mb.local_position_from_global(corner_upper_left);
        std::set<Position> correct_upper_left_corner_link{corner_upper_left_local + pe1,
                                                          corner_upper_left_local - pe2};

        REQUIRE(upper_left_corner_link == correct_upper_left_corner_link);
    }
}

TEST_CASE("Ghosts and no ghosts", "[masked_box][dim2]")
{
    using MaskedBox = reeber::MaskedBox<2>; using Position = MaskedBox::Position;
    using DynPoint = MaskedBox::NewDynamicPoint; using AmrVertexId = reeber::AmrVertexId;

    const int refinement = 1, level = 0, gid = 0;

    const DynPoint e1{1, 0}, e2{0, 1}, one{1, 1}, core_from{3, 3}, core_to{4, 4};
    const DynPoint bounds_from = core_from - one, bounds_to = core_to + one;
    const DynPoint core_shape = core_to - core_from + one, bounds_shape = bounds_to - bounds_from + one;

    MaskedBox mb(core_from, core_to, bounds_from, bounds_to, refinement, level, gid, false);
    MaskedBox mb_no_ghosts(core_from, core_to, core_from, core_to, refinement, level, gid, false);
    MaskedBox mb_2_ghosts(core_from, core_to, bounds_from - one, bounds_to + one, refinement, level, gid, false);

    SECTION("check vertices, simple case: all vertices are active")
    {
        int n_outer = 0;
        // set all inner vertices to be active
        diy::for_each(mb.mask_shape(), [&mb, &n_outer](const Position& p) {
            if (mb.is_outer(p))
            {
                // pretend ghost vertices are masked by another box, does not matter which
                mb.set_mask(p, 2);
                n_outer++;
            } else
                mb.set_mask(p, MaskedBox::ACTIVE);
        });
        REQUIRE(n_outer ==12);

        n_outer = 0;
        diy::for_each(mb_no_ghosts.mask_shape(), [&mb_no_ghosts, &n_outer](const Position& p) {
            if (mb_no_ghosts.is_outer(p))
            {
                // pretend ghost vertices are masked by another box, does not matter which
                mb_no_ghosts.set_mask(p, 2);
                n_outer++;
            } else
                mb_no_ghosts.set_mask(p, MaskedBox::ACTIVE);
        });
        REQUIRE(n_outer == 12);

        n_outer = 0;
        diy::for_each(mb_2_ghosts.mask_shape(), [&mb_2_ghosts, &n_outer](const Position& p) {
            if (mb_2_ghosts.is_outer(p))
            {
                // pretend ghost vertices are masked by another box, does not matter which
                mb_2_ghosts.set_mask(p, 2);
                n_outer++;
            } else
                mb_2_ghosts.set_mask(p, MaskedBox::ACTIVE);
        });
        REQUIRE(n_outer == 12);

        auto vs_mb = mb.vertices();
        std::set<AmrVertexId> vertices_mb(std::begin(vs_mb), std::end(vs_mb));

        auto vs_mb_ng = mb_no_ghosts.vertices();
        std::set<AmrVertexId> vertices_mb_ng(std::begin(vs_mb_ng), std::end(vs_mb_ng));

        auto vs_mb_2g = mb_2_ghosts.vertices();
        std::set<AmrVertexId> vertices_mb_2g(std::begin(vs_mb_2g), std::end(vs_mb_2g));

        std::set<AmrVertexId> correct_vertices { AmrVertexId{gid, 5}, AmrVertexId{gid, 6},
                                                 AmrVertexId{gid, 9}, AmrVertexId{gid, 10}};

        std::set<AmrVertexId> correct_vertices_ng{ AmrVertexId{gid, 0}, AmrVertexId{gid, 1},
                                                   AmrVertexId{gid, 2}, AmrVertexId{gid, 3}};

        std::set<AmrVertexId> correct_vertices_2g{ AmrVertexId{gid, 14}, AmrVertexId{gid, 15},
                                                   AmrVertexId{gid, 20}, AmrVertexId{gid, 21}};

        REQUIRE(vertices_mb == correct_vertices);
        REQUIRE(vertices_mb_ng == correct_vertices_ng);
        REQUIRE(vertices_mb_2g == correct_vertices_2g);

        // outer edges
        auto v_glob = point_from_dynamic_point<2>(core_from);

        auto range_core_from_ng_outer_edge_link = mb_no_ghosts.outer_edge_link(v_glob);
        auto range_core_from_2g_outer_edge_link = mb_2_ghosts.outer_edge_link(v_glob);
        auto range_core_from_outer_edge_link = mb.outer_edge_link(v_glob);

        std::set<Position> core_from_ng_outer_edge_link { std::begin(range_core_from_ng_outer_edge_link),
                                                          std::end(range_core_from_ng_outer_edge_link) };
        std::set<Position> core_from_2g_outer_edge_link { std::begin(range_core_from_2g_outer_edge_link),
                                                          std::end(range_core_from_2g_outer_edge_link) };

        std::set<Position> core_from_outer_edge_link { std::begin(range_core_from_outer_edge_link),
                                                       std::end(range_core_from_outer_edge_link) };

        std::set<Position> core_from_outer_edge_link_correct { point_from_dynamic_point<2>(core_from - e1),
                point_from_dynamic_point<2>(core_from - e2), point_from_dynamic_point<2>(core_from - one) };

        REQUIRE(core_from_ng_outer_edge_link == core_from_outer_edge_link_correct);
        REQUIRE(core_from_2g_outer_edge_link == core_from_outer_edge_link_correct);
        REQUIRE(core_from_outer_edge_link == core_from_outer_edge_link_correct);
    };
}




TEST_CASE("Check blocks constructor in simplest case", "[FabTmtBlock][dim2]")
{
    using Point = diy::DynamicPoint<int, 4>;

    using DynPoint = diy::DynamicPoint<int, 2>;
    using DynPoint4 = diy::DynamicPoint<int, 4>;
    using DBounds = diy::DiscreteBounds;
    using Grid = diy::Grid<double, 2>;
    using IntGrid = diy::Grid<int, 2>;
    using Block = FabComponentBlock<Real, 2>;

    using Vertex = Block::Vertex;
    using MaskedBox = Block::MaskedBox;
    using AMRLink = diy::AMRLink;
    using Bounds = diy::DiscreteBounds;

    using AmrVertexId = reeber::AmrVertexId;
    using AmrEdge = reeber::AmrEdge;
    using Position = MaskedBox::Position;

    //using ConnComponent = Block::Component;

    std::vector<int> levels;

    std::vector<DynPoint> refinements;

    std::vector<Bounds> bounds, cores;
    std::vector<Grid> grids;
    std::vector<Vertex> grid_shapes;
    std::vector<diy::BlockID> bids;
    std::vector<Block> blocks;

    int dim = 2;
    bool negate = false;
    bool absolute = true;
    int dom_size = 3;
    Point min_domain{0, 0, 0, 0};
    Point max_domain{dom_size, dom_size, 0, 0};
    Bounds domain{min_domain, max_domain};

    SECTION("2x2 cell structure, simplest case")
    {

        int blocks_size_1 = 4;
        int n_blocks = 4;

        for(int i = 0; i < n_blocks; ++i)
        {
            levels.push_back(0);
            refinements.push_back(DynPoint({1, 1}));
        }

        Vertex grid_shape{blocks_size_1, blocks_size_1};

        double rho = 20000.0;

        // fill in function values
        for(int i = 0; i < n_blocks; ++i)
        {
            grid_shapes.push_back(grid_shape);
            grids.emplace_back(grid_shape);
        }

        for(int i = 0; i < (int) blocks_size_1; ++i)
            for(int j = 0; j < (int) blocks_size_1; ++j)
                for(int grid_idx = 0; grid_idx < (int) grids.size(); ++grid_idx)
                {
                    Vertex v{i, j};
                    grids[grid_idx](v) = 1.0;
                }


        // 3d one, 4th coord is 0
        DynPoint4 zero{0, 0, 0, 0};
        DynPoint4 one{1, 1, 0, 0};
        DynPoint4 two_x{2, 0, 0, 0};
        DynPoint4 two_y{0, 2, 0, 0};

        cores.emplace_back(diy::DiscreteBounds{zero, one});
        cores.emplace_back(diy::DiscreteBounds{two_x, two_x + one});
        cores.emplace_back(diy::DiscreteBounds{two_y, two_y + one});
        cores.emplace_back(diy::DiscreteBounds{two_x + two_y, two_x + two_y + one});

        for(int i = 0; i < n_blocks; ++i)
        {
            bounds.emplace_back(diy::DiscreteBounds{cores[i].min - one, cores[i].max + one});
        }

        int proc = 0;
        for(int gid = 0; gid < (int) grids.size(); ++gid)
        {
            bids.push_back(diy::BlockID{gid, proc});
        }

        std::vector<AMRLink> links;
        for(int i = 0; i < n_blocks; ++i)
        {
            links.emplace_back(dim, levels[i], refinements[i][0], cores[i], bounds[i]);
            for(int j = 0; j < (int) grids.size(); ++j)
            {
                if (j == i)
                    continue;
                links.back().add_neighbor(bids[j]);
                links.back().add_bounds(levels[j], refinements[j][0], cores[j], bounds[j]);
            }
        }

        for(int i = 0; i < n_blocks; ++i)
        {
            blocks.emplace_back(grids[i], refinements[i][0], levels[i], domain, bounds[i], cores[i], bids[i].gid,
                    &links[i], rho, negate, absolute);
        }

        std::vector<IntGrid> correct_masks;
        for(int i = 0; i < n_blocks; ++i)
        {
            correct_masks.emplace_back(blocks[i].local_.mask_shape());
        }

        for(int i = 0; i < n_blocks; ++i)
        {
            diy::for_each(blocks[i].local_.mask_shape(), [&](const MaskedBox::Position& p) {
                if (blocks[i].local_.is_outer(p))
                    correct_masks[i](p) = MaskedBox::GHOST;
                else
                    correct_masks[i](p) = MaskedBox::ACTIVE;
            });
        }

        int i = 0;

        correct_masks[i](MaskedBox::Position({0, 0})) = 3;
        correct_masks[i](MaskedBox::Position({0, 1})) = 1;
        correct_masks[i](MaskedBox::Position({0, 2})) = 1;
        correct_masks[i](MaskedBox::Position({0, 3})) = 3;
        correct_masks[i](MaskedBox::Position({1, 3})) = 2;
        correct_masks[i](MaskedBox::Position({2, 3})) = 2;
        correct_masks[i](MaskedBox::Position({3, 3})) = 3;
        correct_masks[i](MaskedBox::Position({3, 2})) = 1;
        correct_masks[i](MaskedBox::Position({3, 1})) = 1;
        correct_masks[i](MaskedBox::Position({3, 0})) = 3;
        correct_masks[i](MaskedBox::Position({2, 0})) = 2;
        correct_masks[i](MaskedBox::Position({1, 0})) = 2;

        i++; // i == 1

        correct_masks[i](MaskedBox::Position({0, 0})) = 2;
        correct_masks[i](MaskedBox::Position({0, 1})) = 0;
        correct_masks[i](MaskedBox::Position({0, 2})) = 0;
        correct_masks[i](MaskedBox::Position({0, 3})) = 2;
        correct_masks[i](MaskedBox::Position({1, 3})) = 3;
        correct_masks[i](MaskedBox::Position({2, 3})) = 3;
        correct_masks[i](MaskedBox::Position({3, 3})) = 2;
        correct_masks[i](MaskedBox::Position({3, 2})) = 0;
        correct_masks[i](MaskedBox::Position({3, 1})) = 0;
        correct_masks[i](MaskedBox::Position({3, 0})) = 2;
        correct_masks[i](MaskedBox::Position({2, 0})) = 3;
        correct_masks[i](MaskedBox::Position({1, 0})) = 3;

        i++; // i == 2
        correct_masks[i](MaskedBox::Position({0, 0})) = 1;
        correct_masks[i](MaskedBox::Position({0, 1})) = 3;
        correct_masks[i](MaskedBox::Position({0, 2})) = 3;
        correct_masks[i](MaskedBox::Position({0, 3})) = 1;
        correct_masks[i](MaskedBox::Position({1, 3})) = 0;
        correct_masks[i](MaskedBox::Position({2, 3})) = 0;
        correct_masks[i](MaskedBox::Position({3, 3})) = 1;
        correct_masks[i](MaskedBox::Position({3, 2})) = 3;
        correct_masks[i](MaskedBox::Position({3, 1})) = 3;
        correct_masks[i](MaskedBox::Position({3, 0})) = 1;
        correct_masks[i](MaskedBox::Position({2, 0})) = 0;
        correct_masks[i](MaskedBox::Position({1, 0})) = 0;

        i++; // i == 3
        correct_masks[i](MaskedBox::Position({0, 0})) = 0;
        correct_masks[i](MaskedBox::Position({0, 1})) = 2;
        correct_masks[i](MaskedBox::Position({0, 2})) = 2;
        correct_masks[i](MaskedBox::Position({0, 3})) = 0;
        correct_masks[i](MaskedBox::Position({1, 3})) = 1;
        correct_masks[i](MaskedBox::Position({2, 3})) = 1;
        correct_masks[i](MaskedBox::Position({3, 3})) = 0;
        correct_masks[i](MaskedBox::Position({3, 2})) = 2;
        correct_masks[i](MaskedBox::Position({3, 1})) = 2;
        correct_masks[i](MaskedBox::Position({3, 0})) = 0;
        correct_masks[i](MaskedBox::Position({2, 0})) = 1;
        correct_masks[i](MaskedBox::Position({1, 0})) = 1;

        for(int i = 0; i < n_blocks; ++i)
        {
            REQUIRE(blocks[i].components_.size() == 1);
            if (correct_masks[i] != blocks[i].local_.mask_grid())
            {
                fmt::print("ACHTUNG!, i = {}\n", i);
                for(int kkk = 0; kkk < correct_masks[i].size(); ++kkk)
                {
                    AmrVertexId v;
                    v.gid = i;
                    v.vertex = kkk;
                    MaskedBox::Position local_pos = blocks[i].local_.local_position(v);
                    MaskedBox::Position global_pos = blocks[i].local_.global_position(v);
                    auto my_value = blocks[i].local_.mask_grid()(v);
                    auto correct_value = correct_masks[i](v);
                    if (my_value != correct_value)
                        fmt::print("i = {}, local pos = {}, global pos = {}, correct = {}, my = {}\n", i, local_pos,
                                global_pos, correct_value, my_value);
                }
            }
            REQUIRE(correct_masks[i] == blocks[i].local_.mask_grid());
        }


        // fill this map with edges given by gid and local coordinates
        std::map<std::pair<int, Position>, std::vector<std::pair<int, Position>>> correct_outgoing_edges_pos;

        // will be filled from correct_outgoing_edges_pos with proper AmrEdge-s
        std::map<AmrVertexId, std::set<AmrEdge>> correct_outgoing_edges;

        correct_outgoing_edges_pos[{0, {1, 1}}] = {{2, {1, 2}},
                                                   {3, {2, 2}},
                                                   {1, {2, 1}}};

        for(const auto& vert_vector_pair : correct_outgoing_edges_pos)
        {
            auto inner_v = vert_vector_pair.first;
            int inner_gid = inner_v.first;
            AmrVertexId my_vertex{inner_gid, blocks[inner_gid].local_.local_position_to_vertex(inner_v.second)};
            for(const auto& outer_vertex_pair : vert_vector_pair.second)
            {
                int outer_gid = outer_vertex_pair.first;
                AmrVertexId outer_vertex{outer_gid,
                                         blocks[outer_gid].local_.local_position_to_vertex(outer_vertex_pair.second)};
                correct_outgoing_edges[my_vertex].emplace(AmrEdge{my_vertex, outer_vertex});
            }
        }

        for(const auto& vertex_edge_set_pair : correct_outgoing_edges)
        {
            int gid = vertex_edge_set_pair.first.gid;
            auto v_glob = blocks[gid].local_.global_position(vertex_edge_set_pair.first);
            auto result = get_vertex_edges(v_glob, blocks[gid].local_, &links[gid], domain);
            std::set<AmrEdge> result_set(result.begin(), result.end());
            REQUIRE(vertex_edge_set_pair.second == result_set);
        }
    }
}

TEST_CASE("Check blocks constructor in simplest case-no ghosts", "[FabTmtBlock][dim2]")
{
    using Point = diy::DynamicPoint<int, 4>;

    using DynPoint = diy::DynamicPoint<int, 2>;
    using DynPoint4 = diy::DynamicPoint<int, 4>;
    using DBounds = diy::DiscreteBounds;
    using Grid = diy::Grid<double, 2>;
    using IntGrid = diy::Grid<int, 2>;
    using Block = FabTmtBlock<double, 2>;

    using Vertex = Block::Vertex;
    using MaskedBox = Block::MaskedBox;
    using AMRLink = diy::AMRLink;
    using Bounds = diy::DiscreteBounds;

    using AmrVertexId = reeber::AmrVertexId;
    using AmrEdge = reeber::AmrEdge;
    using Position = MaskedBox::Position;

    //using ConnComponent = Block::Component;

    std::vector<int> levels;

    std::vector<DynPoint> refinements;

    std::vector<Bounds> cores;
    std::vector<Grid> grids;
    std::vector<Vertex> grid_shapes;
    std::vector<diy::BlockID> bids;
    std::vector<Block> blocks;

    int dim = 2;
    bool negate = false;
    bool absolute = true;
    int dom_size = 3;
    Point min_domain{0, 0, 0, 0};
    Point max_domain{dom_size, dom_size, 0, 0};
    Bounds domain{min_domain, max_domain};

    SECTION("2x2 cell structure, simplest case")
    {

        int blocks_size_1 = 4;
        int n_blocks = 4;

        for(int i = 0; i < n_blocks; ++i)
        {
            levels.push_back(0);
            refinements.push_back(DynPoint({1, 1}));
        }

        Vertex grid_shape{blocks_size_1, blocks_size_1};

        double rho = 20000.0;

        // fill in function values
        for(int i = 0; i < n_blocks; ++i)
        {
            grid_shapes.push_back(grid_shape);
            grids.emplace_back(grid_shape);
        }

        for(int i = 0; i < (int) blocks_size_1; ++i)
            for(int j = 0; j < (int) blocks_size_1; ++j)
                for(int grid_idx = 0; grid_idx < (int) grids.size(); ++grid_idx)
                {
                    Vertex v{i, j};
                    grids[grid_idx](v) = 1.0;
                }


        // 3d one, 4th coord is 0
        DynPoint4 zero{0, 0, 0, 0};
        DynPoint4 one{1, 1, 0, 0};
        DynPoint4 two_x{2, 0, 0, 0};
        DynPoint4 two_y{0, 2, 0, 0};

        cores.emplace_back(diy::DiscreteBounds{zero, one});
        cores.emplace_back(diy::DiscreteBounds{two_x, two_x + one});
        cores.emplace_back(diy::DiscreteBounds{two_y, two_y + one});
        cores.emplace_back(diy::DiscreteBounds{two_x + two_y, two_x + two_y + one});

        int proc = 0;
        for(int gid = 0; gid < (int) grids.size(); ++gid)
        {
            bids.push_back(diy::BlockID{gid, proc});
        }

        std::vector<AMRLink> links;
        for(int i = 0; i < n_blocks; ++i)
        {
            links.emplace_back(dim, levels[i], refinements[i][0], cores[i], cores[i]);
            for(int j = 0; j < (int) grids.size(); ++j)
            {
                if (j == i)
                    continue;
                links.back().add_neighbor(bids[j]);
                links.back().add_bounds(levels[j], refinements[j][0], cores[j], cores[j]);
            }
        }

        for(int i = 0; i < n_blocks; ++i)
        {
            blocks.emplace_back(grids[i], refinements[i][0], levels[i], domain, cores[i], cores[i], bids[i].gid,
                    &links[i], rho, negate, absolute);
        }

        std::vector<IntGrid> correct_masks;
        for(int i = 0; i < n_blocks; ++i)
        {
            correct_masks.emplace_back(blocks[i].local_.mask_shape());
        }

        for(int i = 0; i < n_blocks; ++i)
        {
            diy::for_each(blocks[i].local_.mask_shape(), [&](const MaskedBox::Position& p) {
                if (blocks[i].local_.is_outer(p))
                    correct_masks[i](p) = MaskedBox::GHOST;
                else
                    correct_masks[i](p) = MaskedBox::ACTIVE;
            });
        }

        int i = 0;

        correct_masks[i](MaskedBox::Position({0, 0})) = 3;
        correct_masks[i](MaskedBox::Position({0, 1})) = 1;
        correct_masks[i](MaskedBox::Position({0, 2})) = 1;
        correct_masks[i](MaskedBox::Position({0, 3})) = 3;
        correct_masks[i](MaskedBox::Position({1, 3})) = 2;
        correct_masks[i](MaskedBox::Position({2, 3})) = 2;
        correct_masks[i](MaskedBox::Position({3, 3})) = 3;
        correct_masks[i](MaskedBox::Position({3, 2})) = 1;
        correct_masks[i](MaskedBox::Position({3, 1})) = 1;
        correct_masks[i](MaskedBox::Position({3, 0})) = 3;
        correct_masks[i](MaskedBox::Position({2, 0})) = 2;
        correct_masks[i](MaskedBox::Position({1, 0})) = 2;

        i++; // i == 1

        correct_masks[i](MaskedBox::Position({0, 0})) = 2;
        correct_masks[i](MaskedBox::Position({0, 1})) = 0;
        correct_masks[i](MaskedBox::Position({0, 2})) = 0;
        correct_masks[i](MaskedBox::Position({0, 3})) = 2;
        correct_masks[i](MaskedBox::Position({1, 3})) = 3;
        correct_masks[i](MaskedBox::Position({2, 3})) = 3;
        correct_masks[i](MaskedBox::Position({3, 3})) = 2;
        correct_masks[i](MaskedBox::Position({3, 2})) = 0;
        correct_masks[i](MaskedBox::Position({3, 1})) = 0;
        correct_masks[i](MaskedBox::Position({3, 0})) = 2;
        correct_masks[i](MaskedBox::Position({2, 0})) = 3;
        correct_masks[i](MaskedBox::Position({1, 0})) = 3;

        i++; // i == 2
        correct_masks[i](MaskedBox::Position({0, 0})) = 1;
        correct_masks[i](MaskedBox::Position({0, 1})) = 3;
        correct_masks[i](MaskedBox::Position({0, 2})) = 3;
        correct_masks[i](MaskedBox::Position({0, 3})) = 1;
        correct_masks[i](MaskedBox::Position({1, 3})) = 0;
        correct_masks[i](MaskedBox::Position({2, 3})) = 0;
        correct_masks[i](MaskedBox::Position({3, 3})) = 1;
        correct_masks[i](MaskedBox::Position({3, 2})) = 3;
        correct_masks[i](MaskedBox::Position({3, 1})) = 3;
        correct_masks[i](MaskedBox::Position({3, 0})) = 1;
        correct_masks[i](MaskedBox::Position({2, 0})) = 0;
        correct_masks[i](MaskedBox::Position({1, 0})) = 0;

        i++; // i == 3
        correct_masks[i](MaskedBox::Position({0, 0})) = 0;
        correct_masks[i](MaskedBox::Position({0, 1})) = 2;
        correct_masks[i](MaskedBox::Position({0, 2})) = 2;
        correct_masks[i](MaskedBox::Position({0, 3})) = 0;
        correct_masks[i](MaskedBox::Position({1, 3})) = 1;
        correct_masks[i](MaskedBox::Position({2, 3})) = 1;
        correct_masks[i](MaskedBox::Position({3, 3})) = 0;
        correct_masks[i](MaskedBox::Position({3, 2})) = 2;
        correct_masks[i](MaskedBox::Position({3, 1})) = 2;
        correct_masks[i](MaskedBox::Position({3, 0})) = 0;
        correct_masks[i](MaskedBox::Position({2, 0})) = 1;
        correct_masks[i](MaskedBox::Position({1, 0})) = 1;

        for(int i = 0; i < n_blocks; ++i)
        {
            REQUIRE(blocks[i].components_.size() == 1);
            if (correct_masks[i] != blocks[i].local_.mask_grid())
            {
                fmt::print("ACHTUNG!, i = {}\n", i);
                for(int i1 = 0; i1 < correct_masks[i].size(); ++i1)
                {
                    AmrVertexId v;
                    v.gid = i;
                    v.vertex = i1;
                    MaskedBox::Position local_pos = blocks[i].local_.local_position(v);
                    MaskedBox::Position global_pos = blocks[i].local_.global_position(v);
                    auto my_value = blocks[i].local_.mask_grid()(v);
                    auto correct_value = correct_masks[i](v);
                    if (my_value != correct_value)
                        fmt::print("i = {}, local pos = {}, global pos = {}, correct = {}, my = {}\n", i, local_pos,
                                global_pos, correct_value, my_value);
                }
            }
            REQUIRE(correct_masks[i] == blocks[i].local_.mask_grid());
        }


        // fill this map with edges given by gid and local coordinates
        std::map<std::pair<int, Position>, std::vector<std::pair<int, Position>>> correct_outgoing_edges_pos;

        // will be filled from correct_outgoing_edges_pos with proper AmrEdge-s
        std::map<AmrVertexId, std::set<AmrEdge>> correct_outgoing_edges;

        correct_outgoing_edges_pos[{0, {1, 1}}] = {{2, {1, 0}},
                                                   {3, {0, 0}},
                                                   {1, {0, 1}}};

        for(const auto& vert_vector_pair : correct_outgoing_edges_pos)
        {
            auto inner_v = vert_vector_pair.first;
            int inner_gid = inner_v.first;
            AmrVertexId my_vertex{inner_gid, blocks[inner_gid].local_.local_position_to_vertex(inner_v.second)};
            for(const auto& outer_vertex_pair : vert_vector_pair.second)
            {
                int outer_gid = outer_vertex_pair.first;
                AmrVertexId outer_vertex{outer_gid,
                                         blocks[outer_gid].local_.local_position_to_vertex(outer_vertex_pair.second)};
                correct_outgoing_edges[my_vertex].emplace(AmrEdge{my_vertex, outer_vertex});
            }
        }

//        for(size_t i = 0; i < 4; ++i)
//        {
//            AmrVertexId v { 0, i };
//            std::cout << "v = " << v << ", position = " << blocks[v.gid].local_.global_position(v) << std::endl;
//        };

        for(const auto& vertex_edge_set_pair : correct_outgoing_edges)
        {
            std::cout << "HERE, vertex = " << vertex_edge_set_pair.first << std::endl;
            int gid = vertex_edge_set_pair.first.gid;
            auto v_glob = blocks[gid].local_.global_position(vertex_edge_set_pair.first);
            auto result = get_vertex_edges(v_glob, blocks[gid].local_, &links[gid], domain);
            std::set<AmrEdge> result_set(result.begin(), result.end());
            REQUIRE(vertex_edge_set_pair.second == result_set);
        }
    }
}


/*

TEST_CASE("Check blocks constructor in simple masked case", "[FabTmtBlock][dim2][mask]")
{
    using Point = diy::DynamicPoint<int, 4>;
    using Grid = diy::Grid<double, 2>;
    using IntGrid = diy::Grid<int, 2>;
    using Block = FabTmtBlock<double, 2>;

    using Vertex = Block::Vertex;
    using MaskedBox = Block::MaskedBox;
    using AMRLink = diy::AMRLink;
    using Bounds = diy::DiscreteBounds;

    using AmrVertexId = reeber::AmrVertexId;
    using AmrEdge = reeber::AmrEdge;
    using Position = MaskedBox::Position;

    //using ConnComponent = Block::Component;

    std::vector<int> levels, refinements;
    std::vector<Bounds> bounds, cores;
    std::vector<Grid> grids;
    std::vector<Vertex> grid_shapes;
    std::vector<diy::BlockID> bids;
    std::vector<Block> blocks;

    int dim = 2;
    bool negate = false;
    bool absolute = true;
    int dom_size = 3;
    Point min_domain { 0, 0, 0, 0 };
    Point max_domain { dom_size, dom_size, 0, 0 };
    Bounds domain { min_domain, max_domain };

    SECTION("2x2 cell structure, simplest case") {

        int blocks_size_1 = 2;
        int n_blocks = 5;

        for (int i = 0; i < n_blocks; ++i) {
            if (i < n_blocks - 1) {
                levels.push_back(0);
                refinements.push_back(1);
            } else {
                levels.push_back(1);
                refinements.push_back(2);
            }
        }

        Vertex grid_shape { blocks_size_1, blocks_size_1 };
        Vertex mask_shape { grid_shape + Vertex::one() + Vertex::one() };

        double rho = 20000.0;

        // fill in function values
        for (int i = 0; i < n_blocks; ++i) {
            grid_shapes.push_back(mask_shape);
            grids.emplace_back(mask_shape);
        }

        for(int i = 0; i < n_blocks; ++i) {
            diy::for_each(grid_shapes[i], [&](const Vertex v) {
                bool is_outer = false;
                for(int j = 0; j < dim; ++j) {
                    if (v[j] == 0 || v[j] == grid_shapes[i][j] - 1) {
                        is_outer = true;
                        break;
                    }
                }
                grids[i](v) = is_outer ? 2.0 : 1.0;
            });
        }


        // 3d one, 4th coord is 0
        Point zero { 0, 0, 0, 0 };
        Point one { 1, 1, 0, 0 };
        Point two_x { 2, 0, 0, 0 };
        Point two_y { 0, 2, 0, 0 };

        cores.emplace_back(diy::DiscreteBounds { zero, one });
        cores.emplace_back(diy::DiscreteBounds { two_x, two_x + one });
        cores.emplace_back(diy::DiscreteBounds { two_y, two_y + one });
        cores.emplace_back(diy::DiscreteBounds { two_x + two_y, two_x + two_y + one });
        cores.emplace_back(diy::DiscreteBounds { two_x + two_y, two_x + two_y + one });


        for (int i = 0; i < n_blocks; ++i) {
            bounds.emplace_back(diy::DiscreteBounds { cores[i].min - one, cores[i].max + one });
        }


        int proc = 0;
        for (int gid = 0; gid < (int) grids.size(); ++gid) {
            bids.push_back(diy::BlockID { gid, proc });
        }

        std::vector<AMRLink> links;
        for (int i = 0; i < n_blocks; ++i) {
            links.emplace_back(dim, levels[i], refinements[i], cores[i], bounds[i]);
            for (int j = 0; j < (int)grids.size(); ++j) {
                if (j == i)
                    continue;
                links.back().add_neighbor(bids[j]);
                links.back().add_bounds(levels[j], refinements[j], cores[j], bounds[j]);
            }
        }

        for (int i = 0; i < n_blocks; ++i) {
            blocks.emplace_back(grids[i], refinements[i], levels[i], domain, bounds[i], cores[i], bids[i].gid,
                                &links[i], rho, negate, absolute);
        }

        std::vector<IntGrid> correct_masks;
        for (int i = 0; i < n_blocks; ++i) {
            correct_masks.emplace_back(blocks[i].local_.mask_shape());
        }

        for (int i = 0; i < n_blocks; ++i) {
            diy::for_each(blocks[i].local_.mask_shape(), [&](const MaskedBox::Position& p) {
                if (blocks[i].local_.is_outer(p))
                   correct_masks[i](p) = MaskedBox::GHOST;
                else
                   correct_masks[i](p) = MaskedBox::ACTIVE;
            });
        }

        int i = 0;

        correct_masks[i](MaskedBox::Position({0,0})) = 3;
        correct_masks[i](MaskedBox::Position({0,1})) = 1;
        correct_masks[i](MaskedBox::Position({0,2})) = 1;
        correct_masks[i](MaskedBox::Position({0,3})) = 3;
        correct_masks[i](MaskedBox::Position({1,3})) = 2;
        correct_masks[i](MaskedBox::Position({2,3})) = 2;
        correct_masks[i](MaskedBox::Position({3,3})) = 3;
        correct_masks[i](MaskedBox::Position({3,2})) = 1;
        correct_masks[i](MaskedBox::Position({3,1})) = 1;
        correct_masks[i](MaskedBox::Position({3,0})) = 3;
        correct_masks[i](MaskedBox::Position({2,0})) = 2;
        correct_masks[i](MaskedBox::Position({1,0})) = 2;

        correct_masks[i](MaskedBox::Position({2,2})) = 4;

        i++; // i == 1

        correct_masks[i](MaskedBox::Position({0,0})) = 2;
        correct_masks[i](MaskedBox::Position({0,1})) = 0;
        correct_masks[i](MaskedBox::Position({0,2})) = 4;
        correct_masks[i](MaskedBox::Position({0,3})) = 2;
        correct_masks[i](MaskedBox::Position({1,3})) = 3;
        correct_masks[i](MaskedBox::Position({2,3})) = 3;
        correct_masks[i](MaskedBox::Position({3,3})) = 2;
        correct_masks[i](MaskedBox::Position({3,2})) = 0;
        correct_masks[i](MaskedBox::Position({3,1})) = 0;
        correct_masks[i](MaskedBox::Position({3,0})) = 2;
        correct_masks[i](MaskedBox::Position({2,0})) = 3;
        correct_masks[i](MaskedBox::Position({1,0})) = 3;


        i++; // i == 2
        correct_masks[i](MaskedBox::Position({0,0})) = 1;
        correct_masks[i](MaskedBox::Position({0,1})) = 3;
        correct_masks[i](MaskedBox::Position({0,2})) = 3;
        correct_masks[i](MaskedBox::Position({0,3})) = 1;
        correct_masks[i](MaskedBox::Position({1,3})) = 0;
        correct_masks[i](MaskedBox::Position({2,3})) = 0;
        correct_masks[i](MaskedBox::Position({3,3})) = 1;
        correct_masks[i](MaskedBox::Position({3,2})) = 3;
        correct_masks[i](MaskedBox::Position({3,1})) = 3;
        correct_masks[i](MaskedBox::Position({3,0})) = 1;
        correct_masks[i](MaskedBox::Position({2,0})) = 4;
        correct_masks[i](MaskedBox::Position({1,0})) = 0;

        i++; // i == 3
        correct_masks[i](MaskedBox::Position({0,0})) = 4;
        correct_masks[i](MaskedBox::Position({0,1})) = 2;
        correct_masks[i](MaskedBox::Position({0,2})) = 2;
        correct_masks[i](MaskedBox::Position({0,3})) = 0;
        correct_masks[i](MaskedBox::Position({1,3})) = 1;
        correct_masks[i](MaskedBox::Position({2,3})) = 1;
        correct_masks[i](MaskedBox::Position({3,3})) = 0;
        correct_masks[i](MaskedBox::Position({3,2})) = 2;
        correct_masks[i](MaskedBox::Position({3,1})) = 2;
        correct_masks[i](MaskedBox::Position({3,0})) = 0;
        correct_masks[i](MaskedBox::Position({2,0})) = 1;
        correct_masks[i](MaskedBox::Position({1,0})) = 1;

        i++; // i == 4, block on level 1 masking (1,1)
        correct_masks[i](MaskedBox::Position({0,0})) = 0;
        correct_masks[i](MaskedBox::Position({0,1})) = 0;
        correct_masks[i](MaskedBox::Position({0,2})) = 0;
        correct_masks[i](MaskedBox::Position({0,3})) = 2;
        correct_masks[i](MaskedBox::Position({1,3})) = 2;
        correct_masks[i](MaskedBox::Position({2,3})) = 2;
        correct_masks[i](MaskedBox::Position({3,3})) = 3;
        correct_masks[i](MaskedBox::Position({3,2})) = 1;
        correct_masks[i](MaskedBox::Position({3,1})) = 1;
        correct_masks[i](MaskedBox::Position({3,0})) = 1;
        correct_masks[i](MaskedBox::Position({2,0})) = 0;
        correct_masks[i](MaskedBox::Position({1,0})) = 0;


        for (int i = 0; i < n_blocks; ++i) {
            REQUIRE(blocks[i].components_.size() == 1);
            REQUIRE(correct_masks[i] == blocks[i].local_.mask_grid());
//            fmt::print("i = {} OK\n", i);
        }


        // fill this map with edges given by gid and local coordinates
        std::map<std::pair<int, Position>, std::vector<std::pair<int, Position>>> correct_outgoing_edges_pos;

        // will be filled from correct_outgoing_edges_pos with proper AmrEdge-s
        std::map<AmrVertexId, std::set<AmrEdge>> correct_outgoing_edges;

        correct_outgoing_edges_pos[{0, {1,1}}] = { {2, {1, 2}}, { 3, {2, 2} }, { 1, {2, 1}}, {4, {1,1}} };
        correct_outgoing_edges_pos[{0, {1,2}}] = { {1, {2, 1}}, { 1, {2, 2} }, { 2, {1, 1}}, {2, {2,1}}, {4, {1,1}}, {4, {1, 2}} };
        correct_outgoing_edges_pos[{0, {2,1}}] = { {2, {1, 2}}, { 4, {1, 1} }, { 4, {2, 1}}, {1, {1,2}}, {1, {1,1}}, {2, {2,2}} };
        correct_outgoing_edges_pos[{0, {2,2}}] = {};


        correct_outgoing_edges_pos[{2, {1,1}}] = { {1, {2, 2}}, {3, {2,1}}, {0, {1,2}} };
        correct_outgoing_edges_pos[{2, {2,1}}] = { {0, {1, 2}}, {3, {1,2}}, {3, {1,1}}, {4, {1,2}}, {4, {2,2}}};
        correct_outgoing_edges_pos[{2, {1,2}}] = { {3, {2, 1}}, {3, {2,2}}, {0, {1,1}}, {0, {2,1}}};
        correct_outgoing_edges_pos[{2, {2,2}}] = { {0, {2, 1}}, {1, {1,1}}, {3, {1,2}}};


        correct_outgoing_edges_pos[{3, {1,1}}] = { {4, {2,2}}, {2, {2, 1}}, {1, {1, 2}} };

        correct_outgoing_edges_pos[{4, {1,1}}] = { {0, {2,1}}, {0, {1, 1}}, {0, {1, 2}}};
        correct_outgoing_edges_pos[{4, {1,2}}] = { {0, {1,2}}, {2, {2, 1}}};
        correct_outgoing_edges_pos[{4, {2,1}}] = { {0, {2,1}}, {1, {1, 2}}};
        correct_outgoing_edges_pos[{4, {2,2}}] = { {2, {2, 1}}, {3, {1,1}}, {1, {1, 2}}};

        for(const auto& vert_vector_pair : correct_outgoing_edges_pos) {
            auto inner_v = vert_vector_pair.first;
            int inner_gid = inner_v.first;
            AmrVertexId my_vertex { inner_gid, blocks[inner_gid].local_.local_position_to_vertex(inner_v.second) };
            for(const auto& outer_vertex_pair : vert_vector_pair.second) {
                int outer_gid = outer_vertex_pair.first;
                AmrVertexId outer_vertex { outer_gid, blocks[outer_gid].local_.local_position_to_vertex(outer_vertex_pair.second) };
                correct_outgoing_edges[my_vertex].emplace(AmrEdge{my_vertex, outer_vertex});
            }
        }

        for(const auto& vertex_edge_set_pair : correct_outgoing_edges) {
            int gid = vertex_edge_set_pair.first.gid;
            auto v_glob = blocks[gid].local_.global_position(vertex_edge_set_pair.first);
            auto result = get_vertex_edges(v_glob, blocks[gid].local_, &links[gid], domain);
            std::set<AmrEdge> result_set(result.begin(), result.end());
            REQUIRE(vertex_edge_set_pair.second == result_set);
        }
    }
}
*/
