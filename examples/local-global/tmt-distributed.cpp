#ifdef REEBER_USE_TBB
#define DIY_NO_THREADS
#endif

#include <iostream>
#include <string>

#include <dlog/stats.h>
#include <dlog/log.h>
#include <opts/opts.h>

#include <diy/master.hpp>
#include <diy/assigner.hpp>
#include <diy/storage.hpp>
#include <diy/decomposition.hpp>
#include <diy/reduce.hpp>
#include <diy/partners/swap.hpp>
#include <diy/io/block.hpp>

#include "memory.h"

#include "reader-interfaces.h"
#include "triplet-merge-tree-block.h"
#include <reeber/distributed-tmt.h>

// debug only
void save_grids(void* b_, const diy::Master::ProxyWithLink& cp, void*)
{
    TripletMergeTreeBlock* b = static_cast<TripletMergeTreeBlock*>(b_);

    std::string outfn = fmt::format("block-debug-b{}.npy", b->gid);
    diy::mpi::io::file  f(cp.master()->communicator(), outfn, diy::mpi::io::file::create | diy::mpi::io::file::wronly);
    diy::io::NumPy writer(f);
    writer.write_header<Real,TripletMergeTreeBlock::Vertex>(b->grid.shape());
    diy::DiscreteBounds bounds {3};
    for (unsigned i = 0; i < 3; ++i)
    {
        bounds.min[i] = 0;
        bounds.max[i] = b->grid.shape()[i] - 1;
    }
    writer.write(bounds, b->grid.data());
}

void test_link(void* b_, const diy::Master::ProxyWithLink& cp, void*)
{
    TripletMergeTreeBlock* b = static_cast<TripletMergeTreeBlock*>(b_);

    fmt::print("Block {}: {} - {}\n", b->gid, b->local.from(), b->local.to());
    fmt::print("Link of {} -> {}\n", 0, b->local.position(0));
    for(TripletMergeTreeBlock::Index u : b->local.link(0))
        fmt::print("  {} -> {}\n", u, b->local.position(u));
}

int main(int argc, char** argv)
{
    diy::mpi::environment   env(argc, argv);
    diy::mpi::communicator  world;
#ifdef REEBER_USE_BOXLIB_READER
    reeber::io::BoxLib::environment boxlib_env(argc, argv, world);
#endif

    using namespace opts;

    int         nblocks    = world.size();
    std::string prefix     = "./DIY.XXXXXX";
    int         in_memory  = -1;
    int         jobs       = 1;

    std::string profile_path;
    std::string log_level = "info";
    int         threads = r::task_scheduler_init::automatic;

    Options ops(argc, argv);
    ops
        >> Option('b', "blocks",    nblocks,      "number of blocks to use")
        >> Option('m', "memory",    in_memory,    "maximum blocks to store in memory")
        >> Option('j', "jobs",      jobs,         "threads to use during the computation")
        >> Option('s', "storage",   prefix,       "storage prefix")
        >> Option('p', "profile",   profile_path, "path to keep the execution profile")
        >> Option('l', "log",       log_level,    "log level")
        >> Option('t', "threads",   threads,      "number of threads to use (with TBB)")
    ;
    bool        negate      = ops >> Present('n', "negate", "sweep superlevel sets");
    bool        wrap_       = ops >> Present('w', "wrap",   "periodic boundary conditions");
    bool        split       = ops >> Present(     "split",  "use split IO");

    std::string infn, outfn;
    if (  ops >> Present('h', "help", "show help message") ||
        !(ops >> PosOption(infn) >> PosOption(outfn)))
    {
        if (world.rank() == 0)
        {
            fmt::print("Usage: {} INPUT OUT.mt\n", argv[0]);
            fmt::print("Compute local-global tree from NumPy");
#ifdef REEBER_USE_BOXLIB_READER
            fmt::print(" or BoxLib");
#endif
            fmt::print(" input.\n");
            fmt::print("{}", ops);
        }
        return 1;
    }

    r::task_scheduler_init init(threads);

    dlog::add_stream(std::cerr, dlog::severity(log_level))
        << dlog::stamp() << dlog::aux_reporter(world.rank()) << dlog::color_pre() << dlog::level() << dlog::color_post() >> dlog::flush();

    std::ofstream   profile_stream;
    if (profile_path == "-")
        dlog::prof.add_stream(std::cerr);
    else if (!profile_path.empty())
    {
        std::string profile_fn = fmt::format("{}-r{}.prf", profile_path, world.rank());
        profile_stream.open(profile_fn.c_str());
        dlog::prof.add_stream(profile_stream);

        std::string stats_fn = fmt::format("{}-r{}.txt", profile_path, world.rank());
        dlog::stats.open(stats_fn.c_str(), std::ios::out);
    }

    world.barrier();
    dlog::Timer timer;
    LOG_SEV_IF(world.rank() == 0, info) << "Starting computation";

    diy::FileStorage            storage(prefix);

    diy::Master                 master(world,
                                       jobs,
                                       in_memory,
                                       &TripletMergeTreeBlock::create,
                                       &TripletMergeTreeBlock::destroy,
                                       &storage,
                                       &TripletMergeTreeBlock::save,
                                       &TripletMergeTreeBlock::load);

    diy::ContiguousAssigner     assigner(world.size(), nblocks);

    // set up the reader
    Reader* reader_ptr = Reader::create(infn, world);
    Reader& reader  = *reader_ptr;

    diy::DiscreteBounds domain {3};
    domain.min[0] = domain.min[1] = domain.min[2] = 0;
    for (unsigned i = 0; i < reader.shape().size(); ++i)
      domain.max[i] = reader.shape()[i] - 1;

    diy::RegularDecomposer<diy::DiscreteBounds>::BoolVector       share_face(3, false);
    diy::RegularDecomposer<diy::DiscreteBounds>::BoolVector       wrap(3, wrap_);
    diy::RegularDecomposer<diy::DiscreteBounds>                   decomposer(3, domain, assigner.nblocks(), share_face, wrap);
    if (wrap_)
    {
        for (unsigned i = 0; i < 3; ++i)
            if (decomposer.divisions[i] < 3)
            {
                LOG_SEV(fatal) << "Can't have fewer than three divisions per side, when wrap is on";
                return -1;
            }
    }

    decomposer.decompose(world.rank(), assigner,
                         [&](int                          gid,
                            const diy::DiscreteBounds&   core,
                                  diy::DiscreteBounds    bounds,
                            const diy::DiscreteBounds&   domain,
                            const diy::RegularGridLink&  link)
    {
        using Vertex     = TripletMergeTreeBlock::Vertex;
        using Box        = TripletMergeTreeBlock::Box;
        using OffsetGrid = TripletMergeTreeBlock::OffsetGrid;

        TripletMergeTreeBlock*  b  = new TripletMergeTreeBlock;
        diy::RegularGridLink*   l  = new diy::RegularGridLink(link);

        Vertex                  full_shape = Vertex(&domain.max[0]) - Vertex(&domain.min[0]) + Vertex::one();

        LOG_SEV(debug) << "[" << b->gid << "] " << "Bounds: " << Vertex(&bounds.min[0]) << " - " << Vertex(&bounds.max[0]);
        LOG_SEV(debug) << "[" << b->gid << "] " << "Core:   " << Vertex(&core.min[0])   << " - " << Vertex(&core.max[0]);

        OffsetGrid(full_shape, &bounds.min[0], &bounds.max[0]).swap(b->grid);
        reader.read(bounds, b->grid.data(), true);      // collective; implicitly assumes same number of blocks on every processor

        b->gid = gid;
        b->cell_size = reader.cell_size();
        b->mt.set_negate(negate);
        b->local  = Box(full_shape, &core.min[0], &core.max[0]);
        b->global = Box(full_shape, &domain.min[0], &domain.max[0]);
        LOG_SEV(debug) << "[" << b->gid << "] Local box:  " << b->local.from()  << " - " << b->local.to();
        LOG_SEV(debug) << "[" << b->gid << "] Global box: " << b->global.from() << " - " << b->global.to();
        LOG_SEV(debug) << "[" << b->gid << "] Grid shape: " << b->grid.shape();

        int lid   = master.add(gid, b, l);
        static_cast<void>(lid);     // shut up the compiler about lid
    });
    LOG_SEV_IF(world.rank() == 0, info) << "Domain decomposed: " << master.size();
    LOG_SEV_IF(world.rank() == 0, info) << "  (data read)";
    delete reader_ptr;

    world.barrier();
    LOG_SEV_IF(world.rank() == 0, info) << "Time to read data:       " << dlog::clock_to_string(timer.elapsed());
    timer.restart();

    // debug only
    //master.foreach(&save_grids);
    //master.foreach(&test_link);

    auto expand = [wrap_](const TripletMergeTreeBlock::Box& box)
                  {
                      auto expanded = box;
                      if (!wrap_)
                      {
                          for (unsigned i = 0; i < expanded.dimension(); ++i)
                          {
                              if (expanded.from()[i] > 0)
                                  expanded.from()[i]--;
                              if (expanded.to()[i] < expanded.grid_shape()[i] - 1)
                                  expanded.to()[i]++;
                          }
                      }
                      else
                      {
                          expanded.from() -= TripletMergeTreeBlock::Vertex::one();
                          expanded.to()   += TripletMergeTreeBlock::Vertex::one();
                      }

                      return expanded;
                  };

    reeber::compute_merge_tree(master, assigner,
                               &TripletMergeTreeBlock::mt,
                               &TripletMergeTreeBlock::edge_maps,
                               [&expand](TripletMergeTreeBlock* b)                          // topology
                               {
                                   return expand(b->local);
                               },
                               [](TripletMergeTreeBlock* b) -> const decltype(b->grid)&     // function
                               { return b->grid; },
                               [&](TripletMergeTreeBlock* b)                                // gid generator
                               {
                                   auto expanded = expand(b->local);
                                   return [&decomposer,expanded](TripletMergeTreeBlock::Index i)
                                          {
                                            auto p = expanded.position(i);
                                            int gid = decomposer.point_to_gid(p);
                                            return gid;
                                          };
                               });

    // save the result
    timer.restart();
    if (outfn != "none")
    {
    if (!split)
        diy::io::write_blocks(outfn, world, master);
    else
        diy::io::split::write_blocks(outfn, world, master);
    }

    world.barrier();
    LOG_SEV_IF(world.rank() == 0, info) << "Time to output trees:    " << dlog::clock_to_string(timer.elapsed());
    timer.restart();

    dlog::prof.flush();     // TODO: this is necessary because the profile file will close before
                            //       the global dlog::prof goes out of scope and flushes the events.
                            //       Need to eventually fix this.
    dlog::stats.flush();
}

