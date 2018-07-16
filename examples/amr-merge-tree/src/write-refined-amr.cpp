#include <diy/master.hpp>
#include <diy/assigner.hpp>
#include <diy/link.hpp>
#include <diy/grid.hpp>
#include <diy/fmt/format.h>
#include <diy/io/block.hpp>         // for saving blocks in DIY format
#include <diy/io/numpy.hpp>

#include "reeber/amr_helper.h"
#include "fab-block.h"

using Real  = double;
using Block = FabBlock<Real, 3>;

template<class R, unsigned D>
diy::Grid<R, D> refine_grid(const diy::GridRef<R, D>& grid_in, int grid_refinement, int target_refinement)
{
    assert(target_refinement >= grid_refinement and target_refinement % grid_refinement == 0);

    int scale = target_refinement / grid_refinement;

    using Vertex = typename diy::Grid<R, D>::Vertex;

    Vertex in_shape = grid_in.shape();
    Vertex out_shape = scale * in_shape;

    diy::Grid<R, D> result(out_shape);

    R maxval = std::numeric_limits<R>::lowest();
    R minval = std::numeric_limits<R>::max();

    diy::for_each(in_shape, [&result, grid_refinement, target_refinement, grid_in, &maxval, &minval](const Vertex& v) {
        Vertex p_from, p_to;
        std::tie(p_from, p_to) = refine_vertex(v, grid_refinement, target_refinement);
        diy::for_each(p_from, p_to, [&result, &grid_in, &v, &maxval, &minval](const Vertex& p) {
            R val = grid_in(v);
            if (std::isnan(val))
                val = 0;
            maxval = std::max(maxval, val);
            minval = std::min(minval, val);
            result(p) = val; });
    });

//    fmt::print("in refined grid, maxval = {}, minval = {}\n", maxval, minval);
    return result;
}


int main(int argc, char** argv)
{
    if (argc < 2) {
        fmt::print(std::cerr, "Usage: {} IN.amr OUT.npy\n", argv[0]);
        return 1;
    }
    std::string infn = argv[1];
    std::string outfn = argv[2];

    diy::mpi::environment env(argc, argv);
    diy::mpi::communicator world;
    diy::Master master(world,
                       1, -1,
                       Block::create,
            //                       Block::destroy,
                       0,
                       0,
                       Block::save,
                       Block::load);

    diy::ContiguousAssigner assigner(world.size(), 0);
    diy::MemoryBuffer header;

    diy::io::read_blocks(infn, world, assigner, master, header);

    master.foreach([](Block* b, const diy::Master::ProxyWithLink& cp) {
        change_to_c_order<Real, 3>(b);
    });

    diy::DiscreteBounds domain;
    diy::load(header, domain);
    fmt::print("Domain: {} - {}\n", domain.min, domain.max);

    int max_level = 0;
    int max_refinement = 1;

    master.foreach([&max_level, &max_refinement](Block* b, const diy::Master::ProxyWithLink& cp) {
        auto* l = static_cast<diy::AMRLink*>(cp.link());

        max_level = std::max(max_level, l->level());
        max_refinement = std::max(max_refinement, l->refinement());

        fmt::print("{}: level = {}, refinement = {}, shape = {}, core = {} - {}, bounds = {} - {}\n",
                   cp.gid(), l->level(), l->refinement(), b->fab.shape(),
                   l->core().min, l->core().max,
                   l->bounds().min, l->bounds().max);
    });

    diy::mpi::io::file out(world, outfn, diy::mpi::io::file::wronly | diy::mpi::io::file::create);

    diy::io::NumPy writer(out);

    diy::DiscreteBounds refined_domain = refine_bounds<3>(domain, max_refinement);

    writer.write_header<double>(3, refined_domain);

    for (int current_level = 0; current_level <= max_level; ++current_level) {
        master.foreach([&writer, max_refinement, current_level](Block* b, const diy::Master::ProxyWithLink& cp) {
            auto* l = static_cast<diy::AMRLink*>(cp.link());
            if (l->level() == current_level) {
                fmt::print("processing block = {}\n", cp.gid());
                auto refined_grid = refine_grid(b->fab_storage_, l->refinement(), max_refinement);
                int scale = max_refinement / l->refinement();
                diy::DiscreteBounds refined_bounds = refine_bounds<3>(l->bounds(), scale);
                diy::DiscreteBounds refined_core = refine_bounds<3>(l->core(), scale);

                assert(project_point<3>(refined_bounds.max - refined_bounds.min + diy::Point<int, 4>::one()) == refined_grid.shape());

                writer.write(refined_bounds, refined_grid.data(), refined_core);
            }
        });
    }
}
