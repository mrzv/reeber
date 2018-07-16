#include <diy/master.hpp>
#include <diy/io/block.hpp>
#include <opts/opts.h>

#include <reeber/box.h>

#include "fab-block.h"
#include "fab-tmt-block.h"
#include "reader-interfaces.h"
#include "diy/vertices.hpp"
#include "reeber/grid.h"


// block-independent types
using Bounds = diy::DiscreteBounds;
using AmrVertexId = r::AmrVertexId;
using AmrEdge = reeber::AmrEdge;
using AmrEdgeContainer = reeber::AmrEdgeContainer;
using AmrLink = diy::AMRLink;
using BlockID = diy::BlockID;
using DiscreteBounds = diy::DiscreteBounds;


using Block = FabTmtBlock<Real, 3>;
using Vertex = Block::Vertex;
using AmrVertexContainer = Block::AmrVertexContainer;
using GidContainer = Block::GidContainer;
using Component = Block::Component;
using VertexNeighborMap = Block::TripletMergeTree::VertexNeighborMap;
using AmrTripletMergeTree = Block::TripletMergeTree;
using MaskedBox = Block::MaskedBox;


void f1()
{
    using MaskedBox = reeber::MaskedBox<2>;
    using Position = MaskedBox::Position;
    using AmrVertexId = reeber::AmrVertexId;

    const int refinement = 1;
    const int level = 0;
    const int gid = 0;
    const Position core_from { 3, 3 };
    const Position core_to { 4, 4 };
    const Position bounds_from = core_from - Position::one();
    const Position bounds_to = core_to + Position::one();

    const Position e1 { 1, 0 };
    const Position e2 { 0, 1 };

    const Position core_shape = core_to - core_from + Position::one();
    const Position bounds_shape = bounds_to - bounds_from + Position::one();
    const size_t masked_box_size = bounds_shape[0] * bounds_shape[1];
    const int core_size = core_shape[0] * core_shape[1];
    (void)core_size;

    MaskedBox mb(core_from, core_to, bounds_from, bounds_to, refinement, level, gid, false);

    {

        // set all inner vertices to be active
        diy::for_each(mb.mask_shape(), [&mb](const Position& p) {
            if (mb.is_ghost(p))
                // pretend ghost vertices are masked by another box, does not matter which
                mb.set_mask(p, 2);
            else
                mb.set_mask(p, MaskedBox::ACTIVE);
        });

        auto vs = mb.vertices();
        std::set<AmrVertexId> vertices(std::begin(vs), std::end(vs));

        std::set<AmrVertexId> correct_vertices;

        for (size_t i = 0; i < masked_box_size; ++i) {
            AmrVertexId v { gid, i };
            if (not mb.is_ghost(mb.local_position(AmrVertexId { gid, i })))
                correct_vertices.emplace(gid, i);
        }

        assert(vertices.size() - core_size == 0);
        assert(vertices.size() == correct_vertices.size());
        assert(vertices == correct_vertices);
    }

}

void f2()
{
    using Point = diy::Point<int, 4>;
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
    int dom_size = 3;
    Point min_domain { 0, 0, 0, 0 };
    Point max_domain { dom_size, dom_size, 0, 0 };
    Bounds domain { min_domain, max_domain };

    {

        int blocks_size_1 = 2;
        int n_blocks = 4;

        for (int i = 0; i < n_blocks; ++i) {
            levels.push_back(0);
            refinements.push_back(1);
        }

        Vertex grid_shape { blocks_size_1, blocks_size_1 };
        Vertex mask_shape { grid_shape + Vertex::one() + Vertex::one() };

        double rho = 20000.0;

        // fill in function values
        for (int i = 0; i < n_blocks; ++i) {
            grid_shapes.push_back(mask_shape);
            grids.emplace_back(mask_shape);
        }

        for (int grid_idx = 0; grid_idx < n_blocks; ++grid_idx) {
            diy::for_each(grid_shapes[grid_idx], [&](const Vertex v) {
                bool is_ghost = false;
                for (int j = 0; j < dim; ++j) {
                    if (v[j] == 0 || v[j] == grid_shapes[grid_idx][j] - 1) {
                        is_ghost = true;
                        break;
                    }
                }
                grids[grid_idx](v) = is_ghost ? 2.0 : 1.0;
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
            for (int j = 0; j < (int) grids.size(); ++j) {
                if (j == i)
                    continue;
                links.back().add_neighbor(bids[j]);
                links.back().add_bounds(levels[j], refinements[j], cores[j], bounds[j]);
            }
        }

        for (int i = 0; i < n_blocks; ++i) {
            blocks.emplace_back(grids[i], refinements[i], levels[i], domain, bounds[i], cores[i], bids[i].gid,
                                &links[i], rho, negate);
        }

        std::vector<IntGrid> correct_masks;
        for (int i = 0; i < n_blocks; ++i) {
            correct_masks.emplace_back(blocks[i].local_.mask_shape());
        }

        for (int i = 0; i < n_blocks; ++i) {
            diy::for_each(blocks[i].local_.mask_shape(), [&](const MaskedBox::Position& p) {
                if (blocks[i].local_.is_ghost(p))
                    correct_masks[i](p) = MaskedBox::GHOST;
                else
                    correct_masks[i](p) = MaskedBox::ACTIVE;
            });
        }

        int i = 0;

        correct_masks[i](MaskedBox::Position({ 0, 0 })) = 3;
        correct_masks[i](MaskedBox::Position({ 0, 1 })) = 1;
        correct_masks[i](MaskedBox::Position({ 0, 2 })) = 1;
        correct_masks[i](MaskedBox::Position({ 0, 3 })) = 3;
        correct_masks[i](MaskedBox::Position({ 1, 3 })) = 2;
        correct_masks[i](MaskedBox::Position({ 2, 3 })) = 2;
        correct_masks[i](MaskedBox::Position({ 3, 3 })) = 3;
        correct_masks[i](MaskedBox::Position({ 3, 2 })) = 1;
        correct_masks[i](MaskedBox::Position({ 3, 1 })) = 1;
        correct_masks[i](MaskedBox::Position({ 3, 0 })) = 3;
        correct_masks[i](MaskedBox::Position({ 2, 0 })) = 2;
        correct_masks[i](MaskedBox::Position({ 1, 0 })) = 2;

        i++; // i == 1

        correct_masks[i](MaskedBox::Position({ 0, 0 })) = 2;
        correct_masks[i](MaskedBox::Position({ 0, 1 })) = 0;
        correct_masks[i](MaskedBox::Position({ 0, 2 })) = 0;
        correct_masks[i](MaskedBox::Position({ 0, 3 })) = 2;
        correct_masks[i](MaskedBox::Position({ 1, 3 })) = 3;
        correct_masks[i](MaskedBox::Position({ 2, 3 })) = 3;
        correct_masks[i](MaskedBox::Position({ 3, 3 })) = 2;
        correct_masks[i](MaskedBox::Position({ 3, 2 })) = 0;
        correct_masks[i](MaskedBox::Position({ 3, 1 })) = 0;
        correct_masks[i](MaskedBox::Position({ 3, 0 })) = 2;
        correct_masks[i](MaskedBox::Position({ 2, 0 })) = 3;
        correct_masks[i](MaskedBox::Position({ 1, 0 })) = 3;


        i++; // i == 2
        correct_masks[i](MaskedBox::Position({ 0, 0 })) = 1;
        correct_masks[i](MaskedBox::Position({ 0, 1 })) = 3;
        correct_masks[i](MaskedBox::Position({ 0, 2 })) = 3;
        correct_masks[i](MaskedBox::Position({ 0, 3 })) = 1;
        correct_masks[i](MaskedBox::Position({ 1, 3 })) = 0;
        correct_masks[i](MaskedBox::Position({ 2, 3 })) = 0;
        correct_masks[i](MaskedBox::Position({ 3, 3 })) = 1;
        correct_masks[i](MaskedBox::Position({ 3, 2 })) = 3;
        correct_masks[i](MaskedBox::Position({ 3, 1 })) = 3;
        correct_masks[i](MaskedBox::Position({ 3, 0 })) = 1;
        correct_masks[i](MaskedBox::Position({ 2, 0 })) = 0;
        correct_masks[i](MaskedBox::Position({ 1, 0 })) = 0;

        i++; // i == 3
        correct_masks[i](MaskedBox::Position({ 0, 0 })) = 0;
        correct_masks[i](MaskedBox::Position({ 0, 1 })) = 2;
        correct_masks[i](MaskedBox::Position({ 0, 2 })) = 2;
        correct_masks[i](MaskedBox::Position({ 0, 3 })) = 0;
        correct_masks[i](MaskedBox::Position({ 1, 3 })) = 1;
        correct_masks[i](MaskedBox::Position({ 2, 3 })) = 1;
        correct_masks[i](MaskedBox::Position({ 3, 3 })) = 0;
        correct_masks[i](MaskedBox::Position({ 3, 2 })) = 2;
        correct_masks[i](MaskedBox::Position({ 3, 1 })) = 2;
        correct_masks[i](MaskedBox::Position({ 3, 0 })) = 0;
        correct_masks[i](MaskedBox::Position({ 2, 0 })) = 1;
        correct_masks[i](MaskedBox::Position({ 1, 0 })) = 1;


        for (int i = 0; i < n_blocks; ++i) {
            assert(blocks[i].components_.size() == 1);
            assert(correct_masks[i] == blocks[i].local_.mask_grid());
        }

        for (int block_idx = 0; block_idx < n_blocks; ++block_idx) {
            auto& curr_block = blocks[block_idx];
            for (size_t i = 0; i < curr_block.local_.mask_grid().size(); ++i) {
                auto p_loc = curr_block.local_.local_position(reeber::AmrVertexId { curr_block.gid, i });
                auto p_glob = curr_block.local_.global_position(reeber::AmrVertexId { curr_block.gid, i });
                fmt::print("in block bound_from = {} bound_to = {}, i = {}, local pos = {}, global pos = {}\n",
                           bounds[block_idx].min, bounds[block_idx].max, i,
                           p_loc, p_glob);
            }
        }


        // fill this map with edges given by gid and local coordinates
        std::map<std::pair<int, Position>, std::vector<std::pair<int, Position>>> correct_outgoing_edges_pos;

        // will be filled from correct_outgoing_edges_pos with proper AmrEdge-s
        std::map<AmrVertexId, std::set<AmrEdge>> correct_outgoing_edges;

        correct_outgoing_edges_pos[{ 0, { 1, 1 }}] = {{ 2, { 1, 2 }},
                                                      { 3, { 2, 2 }},
                                                      { 1, { 2, 1 }}};

        for (const auto& vert_vector_pair : correct_outgoing_edges_pos) {
            auto inner_v = vert_vector_pair.first;
            int inner_gid = inner_v.first;
            AmrVertexId my_vertex { inner_gid, blocks[inner_gid].local_.local_position_to_vertex(inner_v.second) };
            for (const auto& outer_vertex_pair : vert_vector_pair.second) {
                int outer_gid = outer_vertex_pair.first;
                AmrVertexId outer_vertex { outer_gid, blocks[outer_gid].local_.local_position_to_vertex(
                        outer_vertex_pair.second) };
                correct_outgoing_edges[my_vertex].emplace(AmrEdge { my_vertex, outer_vertex });
            }
        }

        for (const auto& vertex_edge_set_pair : correct_outgoing_edges) {
            int gid = vertex_edge_set_pair.first.gid;
            auto v_glob = blocks[gid].local_.global_position(vertex_edge_set_pair.first);
            auto result = get_vertex_edges(v_glob, blocks[gid].local_, &links[gid], domain);
            std::set<AmrEdge> result_set(result.begin(), result.end());
            assert(vertex_edge_set_pair.second == result_set);
        }
    }
}
void f3()
{
    using Point = diy::Point<int, 4>;
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
    int dom_size = 3;
    Point min_domain { 0, 0, 0, 0 };
    Point max_domain { dom_size, dom_size, 0, 0 };
    Bounds domain { min_domain, max_domain };

    {

        int blocks_size_1 = 2;
        int n_blocks = 4;

        for (int i = 0; i < n_blocks; ++i) {
            levels.push_back(0);
            refinements.push_back(1);
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
                bool is_ghost = false;
                for(int j = 0; j < dim; ++j) {
                    if (v[j] == 0 || v[j] == grid_shapes[i][j] - 1) {
                        is_ghost = true;
                        break;
                    }
                }
                grids[i](v) = is_ghost ? 2.0 : 1.0;
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
                                &links[i], rho, negate);
        }

        std::vector<IntGrid> correct_masks;
        for (int i = 0; i < n_blocks; ++i) {
            correct_masks.emplace_back(blocks[i].local_.mask_shape());
        }

        for (int i = 0; i < n_blocks; ++i) {
            diy::for_each(blocks[i].local_.mask_shape(), [&](const MaskedBox::Position& p) {
                if (blocks[i].local_.is_ghost(p))
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

        i++; // i == 1

        correct_masks[i](MaskedBox::Position({0,0})) = 2;
        correct_masks[i](MaskedBox::Position({0,1})) = 0;
        correct_masks[i](MaskedBox::Position({0,2})) = 0;
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
        correct_masks[i](MaskedBox::Position({2,0})) = 0;
        correct_masks[i](MaskedBox::Position({1,0})) = 0;

        i++; // i == 3
        correct_masks[i](MaskedBox::Position({0,0})) = 0;
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


        for (int i = 0; i < n_blocks; ++i) {
            assert(blocks[i].components_.size() == 1);
            assert(correct_masks[i] == blocks[i].local_.mask_grid());
        }


        // fill this map with edges given by gid and local coordinates
        std::map<std::pair<int, Position>, std::vector<std::pair<int, Position>>> correct_outgoing_edges_pos;

        // will be filled from correct_outgoing_edges_pos with proper AmrEdge-s
        std::map<AmrVertexId, std::set<AmrEdge>> correct_outgoing_edges;

        correct_outgoing_edges_pos[{0, {1,1}}] = { {2, {1, 2}}, { 3, {2, 2} }, { 1, {2, 1}} };

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
            assert(vertex_edge_set_pair.second == result_set);
        }
    }
}

void f4()
{
    using Point = diy::Point<int, 4>;
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
    int dom_size = 3;
    Point min_domain { 0, 0, 0, 0 };
    Point max_domain { dom_size, dom_size, 0, 0 };
    Bounds domain { min_domain, max_domain };

    {

        int blocks_size_1 = 2;
        int n_blocks = 4;

        for (int i = 0; i < n_blocks; ++i) {
            levels.push_back(0);
            refinements.push_back(1);
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
                bool is_ghost = false;
                for(int j = 0; j < dim; ++j) {
                    if (v[j] == 0 || v[j] == grid_shapes[i][j] - 1) {
                        is_ghost = true;
                        break;
                    }
                }
                grids[i](v) = is_ghost ? 2.0 : 1.0;
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
                                &links[i], rho, negate);
        }

        std::vector<IntGrid> correct_masks;
        for (int i = 0; i < n_blocks; ++i) {
            correct_masks.emplace_back(blocks[i].local_.mask_shape());
        }

        for (int i = 0; i < n_blocks; ++i) {
            diy::for_each(blocks[i].local_.mask_shape(), [&](const MaskedBox::Position& p) {
                if (blocks[i].local_.is_ghost(p))
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

        i++; // i == 1

        correct_masks[i](MaskedBox::Position({0,0})) = 2;
        correct_masks[i](MaskedBox::Position({0,1})) = 0;
        correct_masks[i](MaskedBox::Position({0,2})) = 0;
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
        correct_masks[i](MaskedBox::Position({2,0})) = 0;
        correct_masks[i](MaskedBox::Position({1,0})) = 0;

        i++; // i == 3
        correct_masks[i](MaskedBox::Position({0,0})) = 0;
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


        for (int i = 0; i < n_blocks; ++i) {
            assert(blocks[i].components_.size() == 1);
            assert(correct_masks[i] == blocks[i].local_.mask_grid());
        }


        // fill this map with edges given by gid and local coordinates
        std::map<std::pair<int, Position>, std::vector<std::pair<int, Position>>> correct_outgoing_edges_pos;

        // will be filled from correct_outgoing_edges_pos with proper AmrEdge-s
        std::map<AmrVertexId, std::set<AmrEdge>> correct_outgoing_edges;

        correct_outgoing_edges_pos[{0, {1,1}}] = { {2, {1, 2}}, { 3, {2, 2} }, { 1, {2, 1}} };

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
            assert(vertex_edge_set_pair.second == result_set);
        }
    }
}

void f5()
{
    using Point = diy::Point<int, 4>;
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
    int dom_size = 3;
    Point min_domain { 0, 0, 0, 0 };
    Point max_domain { dom_size, dom_size, 0, 0 };
    Bounds domain { min_domain, max_domain };

    {

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
                bool is_ghost = false;
                for(int j = 0; j < dim; ++j) {
                    if (v[j] == 0 || v[j] == grid_shapes[i][j] - 1) {
                        is_ghost = true;
                        break;
                    }
                }
                grids[i](v) = is_ghost ? 2.0 : 1.0;
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
                                &links[i], rho, negate);
        }

        std::vector<IntGrid> correct_masks;
        for (int i = 0; i < n_blocks; ++i) {
            correct_masks.emplace_back(blocks[i].local_.mask_shape());
        }

        for (int i = 0; i < n_blocks; ++i) {
            diy::for_each(blocks[i].local_.mask_shape(), [&](const MaskedBox::Position& p) {
                if (blocks[i].local_.is_ghost(p))
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
            assert(blocks[i].components_.size() == 1);
            assert(correct_masks[i] == blocks[i].local_.mask_grid());
            fmt::print("i = {} OK\n", i);
        }


        // fill this map with edges given by gid and local coordinates
        std::map<std::pair<int, Position>, std::vector<std::pair<int, Position>>> correct_outgoing_edges_pos;

        // will be filled from correct_outgoing_edges_pos with proper AmrEdge-s
        std::map<AmrVertexId, std::set<AmrEdge>> correct_outgoing_edges;

        correct_outgoing_edges_pos[{0, {1,1}}] = { {2, {1, 2}}, { 3, {2, 2} }, { 1, {2, 1}}, {4, {1,1}} };
        correct_outgoing_edges_pos[{0, {1,2}}] = { {1, {2, 1}}, { 1, {2, 2} }, { 2, {1, 1}}, {3, {1,1}} };
        correct_outgoing_edges_pos[{0, {2,1}}] = { {3, {1, 2}}, { 4, {1, 1} }, { 4, {2, 1}}, {1, {1,2}}, {1, {1,1}} };
        correct_outgoing_edges_pos[{0, {2,2}}] = {};

        correct_outgoing_edges_pos[{4, {1,1}}] = { {2, {1, 2}}, { 3, {2, 2} }, { 1, {2, 1}}, {4, {1,1}} };

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
            assert(vertex_edge_set_pair.second == result_set);
        }
    }
}

int main()
{
    f5();
}

template<unsigned D>
void check_edge_symmetry(const typename std::vector<FabTmtBlock<double, D>>& blocks)
{
    using BID = diy::BlockID;
    std::map<BID, std::map<BID, std::set<AmrEdge>>> edge_matr;

    size_t n_blocks = blocks.size();

    for (unsigned i = 0; i < n_blocks; ++i) {
        BID from_bid = blocks[i].bid;
        for (const auto& cc : blocks[i].comps_) {
            for (const auto& e : cc.outgoing_edges_) {
                BID to_bid = std::get<1>(e).bid;
                if (from_bid < to_bid) {
                    edge_matr[from_bid][to_bid].insert(e);
                } else {
                    assert(to_bid < from_bid);
                    edge_matr[from_bid][to_bid].emplace(std::get<1>(e), std::get<0>(e));
                }
            }
        }
    }

    for (size_t i = 0; i < n_blocks; ++i)
        for (size_t j = i + 1; j < n_blocks; ++j) {
            BID bid_i = blocks[i].bid;
            BID bid_j = blocks[j].bid;
            fmt::print("Checking edges from {} to {}: {} vs {}\n", bid_i.gid, bid_j.gid, edge_matr[bid_i][bid_j].size(),
                       edge_matr[bid_j][bid_i].size());
            assert(edge_matr[bid_i][bid_j] == edge_matr[bid_j][bid_i]);
        }

    fmt::print("Outgoing edges are symmetric\n");
}
