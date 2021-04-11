#pragma once

#include <string>

#include <reeber-real.h>

#include <diy/master.hpp>
#include <diy/vertices.hpp>
#include <diy/assigner.hpp>
#include <diy/decomposition.hpp>
#include <diy/link.hpp>
#include <diy/io/block.hpp>         // for saving blocks in DIY format
#include <diy/io/numpy.hpp>

#include <reeber/format.h>

#include <dlog/stats.h>
#include <dlog/log.h>

#include "fab-block.h"

template<unsigned D>
void read_from_npy_file(std::string infn,
                        diy::mpi::communicator& world,
                        int nblocks,
                        diy::Master& master_reader,
                        diy::ContiguousAssigner& assigner,
                        diy::MemoryBuffer& header,
                        diy::DiscreteBounds& domain,
                        diy::RegularDecomposer<diy::DiscreteBounds>::BoolVector wrap)
{
    using FabBlockR = FabBlock<Real, D>;

    diy::mpi::io::file in(world, infn, diy::mpi::io::file::rdonly);
    diy::io::NumPy reader(in);
    reader.read_header();

    if (reader.word_size() != sizeof(Real))
    {
        LOG_SEV_IF(world.rank() == 0, fatal) << "Type mismatch: in numpy " << reader.word_size() << ", expect floating point type of size " << sizeof(Real);
        throw std::runtime_error("Type mismatch");
    }

    using Decomposer = diy::RegularDecomposer<diy::DiscreteBounds>;
    using Point   = diy::DynamicPoint<int, DIY_MAX_DIM>;

    domain.min = { 0, 0, 0, 0 };
    domain.max = { 0, 0, 0, 0 };
    Point one = Point::one(D);

//    for(unsigned i = D; i < DIY_MAX_DIM; ++i)
//        one[i] = 0;

    for(unsigned i = 0; i < D; ++i)
    {
        domain.max[i] = reader.shape()[i] - 1;
    }

    Decomposer decomposer(D, domain, nblocks,
                          Decomposer::BoolVector { false, false, false },   // share_face
                          wrap,
                          Decomposer::CoordinateVector { 0, 0, 0 });        // no ghosts -- cannot read ghost data correctly

    decomposer.decompose(world.rank(), assigner, [&master_reader, &wrap, &reader, one](int gid,
                                                                       const Decomposer::Bounds& core,
                                                                       const Decomposer::Bounds& bounds,
                                                                       const Decomposer::Bounds& domain,
                                                                       const Decomposer::Link& link) {
        auto* b = new FabBlockR;

        auto shape_4d = bounds.max - bounds.min + one;
        typename FabBlockR::Shape shape(&shape_4d[0]);       // quick and hacky
        bool c_order = true;
        b->fab_storage_ = decltype(b->fab_storage_)(shape, c_order);
        b->fab = decltype(b->fab)(b->fab_storage_.data(), shape, c_order);

#ifdef ZARIJA
        // pretend that the 1st and only field is particle_mass_density
        // for compatability with integral printing code
        b->extra_names_.push_back("particle_mass_density");
        b->extra_fabs_.push_back(b->fab);
#endif

        reader.read(core, b->fab.data());

        // copy link
        diy::AMRLink* amr_link = new diy::AMRLink(D, 0, 1, link.core(), bounds);
        for (int i = 0; i < link.size(); ++i)
        {
            amr_link->add_neighbor(link.target(i));
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
