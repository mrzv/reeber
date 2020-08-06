#include "../include/read-hdf5.h"

#include <reeber-real.h>

#include <diy/master.hpp>
#include <diy/vertices.hpp>
#include <diy/assigner.hpp>
#include <diy/decomposition.hpp>
#include <diy/link.hpp>
#include <reeber/format.h>

#include <highfive/H5File.hpp>
namespace h5 = HighFive;

#include <dlog/stats.h>
#include <dlog/log.h>

#include "fab-block.h"

void read_from_hdf5_file(std::string infn,
                         std::vector<std::string> all_var_names, // HDF5 only: all fields that will be read from plotfile
                         int n_mt_vars,                          // sum of first n_mt_vars in all_var_names will be stored in fab of FabBlock,
                                                                 // for each variable listed in all_var_names FabBlock will have an extra GridRef
                         diy::mpi::communicator& world,
                         int nblocks,
                         diy::Master& master_reader,
                         diy::ContiguousAssigner& assigner,
                         diy::MemoryBuffer& header,
                         diy::DiscreteBounds& domain)
{
    constexpr unsigned D = 3;
    using FabBlockR = FabBlock<Real, D>;

    h5::File in(infn);

    // get shape
    std::vector<h5::DataSet> datasets;
    for (auto name : all_var_names)
        datasets.emplace_back(in.getDataSet(name));
    auto dimensions = datasets[0].getDimensions();

    using Decomposer = diy::RegularDecomposer<diy::DiscreteBounds>;
    using Point   = diy::DynamicPoint<int, DIY_MAX_DIM>;

    Point one = Point::one(D);
    for(unsigned i = 0; i < D; ++i)
    {
        domain.min[i] = 0;
        domain.max[i] = dimensions[i] - 1;
    }

    Decomposer::BoolVector wrap { true, true, true };               // TODO
    Decomposer decomposer(D, domain, nblocks,
                          Decomposer::BoolVector { false, false, false },   // share_face
                          wrap,
                          Decomposer::CoordinateVector { 0, 0, 0 });        // ghosts

    decomposer.decompose(world.rank(), assigner, [&master_reader, &wrap, &datasets, &all_var_names, n_mt_vars, one]
                                                                      (int gid,
                                                                       const Decomposer::Bounds& core,
                                                                       const Decomposer::Bounds& bounds,
                                                                       const Decomposer::Bounds& domain,
                                                                       const Decomposer::Link& link) {
        auto* b = new FabBlockR;

        // we never want ghosts
        auto my_bounds = core;

        auto shape_4d = my_bounds.max - my_bounds.min + one;
        typename FabBlockR::Shape shape(&shape_4d[0]);       // quick and hacky
        bool c_order = true;
        b->fab_storage_ = decltype(b->fab_storage_)(shape, c_order);
        b->fab = decltype(b->fab)(b->fab_storage_.data(), shape, c_order);

//#ifdef ZARIJA
//        // pretend that the 1st and only field is particle_mass_density
//        // for compatability with integral printing code
//        b->extra_names_.push_back("particle_mass_density");
//        b->extra_fabs_.push_back(b->fab);
//#endif

        diy::Grid<Real, D> core_grid(shape);

        std::vector<size_t> from(D), size(D);
        for (unsigned i = 0; i < D; ++i)
        {
            from[i] = core.min[i];
            size[i] = core.max[i] - core.min[i] + 1;
        }
        for (size_t i = 0; i < all_var_names.size(); ++i)
        {
            Real* extra_fab = new Real[b->fab.size()];
            b->extra_fabs_.emplace_back(extra_fab, shape, c_order);
            b->extra_names_.push_back(all_var_names[i]);

            datasets[i].select(from, size).read(core_grid.data());

            auto& g = b->extra_fabs_.back();
            diy::for_each(shape, [&](const typename FabBlockR::Shape& p) {
                g(p) = core_grid(p);
            });
        }

        // fill the field for MT computation
        diy::for_each(shape, [&](const typename FabBlockR::Shape& p) {
            b->fab(p) = 0;
        });
        for (int i = 0; i < n_mt_vars; ++i)
            diy::for_each(shape, [&](const typename FabBlockR::Shape& p) {
                b->fab(p) += b->extra_fabs_[i](p);
            });

        // copy link
        diy::AMRLink* amr_link = new diy::AMRLink(D, 0, 1, link.core(), my_bounds);
        for (int i = 0; i < link.size(); ++i)
        {
            amr_link->add_neighbor(link.target(i));
            // shrink core from bounds, since it's not stored in the RegularLink explicitly
            auto nbr_core = link.bounds(i);
            auto nbr_bounds = link.bounds(i);
            amr_link->add_bounds(0, 1, nbr_core, nbr_bounds);
        }

        // record wrap
        for (int dir_x : { -1, 0, 1 })
        {
            if (!wrap[0] && dir_x) continue;
            if (dir_x < 0 && core.min[0] != domain.min[0]) continue;
            if (dir_x > 0 && core.max[0] != domain.max[0]) continue;

            for (int dir_y : { -1, 0, 1 })
            {
                if (!wrap[1] && dir_y) continue;
                if (dir_y < 0 && core.min[1] != domain.min[1]) continue;
                if (dir_y > 0 && core.max[1] != domain.max[1]) continue;

                for (int dir_z : { -1, 0, 1 })
                {
                    if (dir_x == 0 and dir_y == 0 and dir_z == 0)
                        continue;

                    if (!wrap[2] && dir_z) continue;
                    if (dir_z < 0 && core.min[2] != domain.min[2]) continue;
                    if (dir_z > 0 && core.max[2] != domain.max[2]) continue;
        // TODO! fix wrap
                    amr_link->add_wrap(diy::Direction{ dir_x, dir_y, dir_z });
                }
            }
        }
        master_reader.add(gid, b, amr_link);
    });
}
