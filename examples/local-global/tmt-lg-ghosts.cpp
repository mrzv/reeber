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

#include "format.h"
#include "memory.h"

#include "reader-interfaces.h"
#include "ghosts.h"
#include "triplet-merge-tree-block.h"

// Load the specified chunk of data, compute local merge tree, add block to diy::Master
struct LoadAdd
{
    typedef         TripletMergeTreeBlock::Vertex                   Vertex;
    typedef         TripletMergeTreeBlock::TripletMergeTree         TripletMergeTree;
    typedef         TripletMergeTreeBlock::Box                      Box;
    typedef         TripletMergeTreeBlock::OffsetGrid               OffsetGrid;

                    LoadAdd(diy::Master& master_, const Reader& reader_, bool negate_, bool wrap_):
                        master(&master_), reader(reader_), negate(negate_), wrap(wrap_)      {}

    void            operator()(int                          gid,
                               const diy::DiscreteBounds&   core,
                                     diy::DiscreteBounds    bounds,
                               const diy::DiscreteBounds&   domain,
                               const diy::RegularGridLink&  link) const
    {
        TripletMergeTreeBlock*  b  = new TripletMergeTreeBlock;
        diy::RegularGridLink*   l  = new diy::RegularGridLink(link);

        Vertex                  full_shape = Vertex(domain.max) - Vertex(domain.min) + Vertex::one();

        LOG_SEV(debug) << "[" << b->gid << "] " << "Bounds: " << Vertex(bounds.min) << " - " << Vertex(bounds.max);
        LOG_SEV(debug) << "[" << b->gid << "] " << "Core:   " << Vertex(core.min)   << " - " << Vertex(core.max);

        OffsetGrid(full_shape, bounds.min, bounds.max).swap(b->grid);
        reader.read(bounds, b->grid.data(), true);      // collective; implicitly assumes same number of blocks on every processor

        b->gid = gid;
        b->cell_size = reader.cell_size();
        b->mt.set_negate(negate);
        b->core  = Box(full_shape, core.min, core.max);
        b->local = b->global = Box(full_shape, bounds.min, bounds.max);
        LOG_SEV(debug) << "[" << b->gid << "] Local box:  " << b->local.from()  << " - " << b->local.to();
        LOG_SEV(debug) << "[" << b->gid << "] Global box: " << b->global.from() << " - " << b->global.to();
        LOG_SEV(debug) << "[" << b->gid << "] Grid shape: " << b->grid.shape();

        int lid   = master->add(gid, b, l);
        static_cast<void>(lid);     // shut up the compiler about lid
    }

    diy::Master*            master;
    const Reader&           reader;
    bool                    negate;
    bool                    wrap;
};

void record_stats(const char* message, const char* format, fmt::ArgList args)
{
    fmt::print(dlog::stats, "{:<25} ", message);
    fmt::print(dlog::stats, format, args);
    fmt::print(dlog::stats, " (hwm = {})", proc_status_value("VmHWM"));
    fmt::print(dlog::stats, "\n");
    dlog::prof.flush();
    dlog::stats.flush();
}
FMT_VARIADIC(void, record_stats, const char*, const char*)

void compute_tree(void* b_, const diy::Master::ProxyWithLink& cp, void* aux)
{
    TripletMergeTreeBlock* b = static_cast<TripletMergeTreeBlock*>(b_);

    record_stats("Local box:", "{}", b->local);
    r::compute_merge_tree(b->mt, b->local, b->grid);

    LOG_SEV(debug) << "[" << b->gid << "] " << "Initial tree size: " << b->mt.size();
    record_stats("Initial tree size:", "{}", b->mt.size());

    TripletMergeTreeBlock::OffsetGrid().swap(b->grid);     // clear out the grid, we don't need it anymore
}

void save_no_vertices(diy::BinaryBuffer& bb, const TripletMergeTreeBlock::TripletMergeTree& mt)
{
    reeber::Serialization<TripletMergeTreeBlock::TripletMergeTree>::save(bb, mt);
}

struct GlobalBoundary
{
    typedef     TripletMergeTreeBlock::Box             Box;

                GlobalBoundary(const Box& global_, bool wrap_):
                    global(global_), wrap(wrap_)                                    {}

    bool        operator()(TripletMergeTreeBlock::Index v) const
    {
        Box::Position       vp          = global.position(v);
        const Box::Position full_shape  = global.grid_shape();
        for (size_t i = 0; i < global.dimension(); ++i)
        {
            if      ( wrap &&  global.from()[i] == 0 && global.to()[i] == full_shape[i])
                continue;
            else if (!wrap && (   (global.from()[i] == 0               && vp[i] == global.from()[i])
                               || (global.to()[i] == full_shape[i] - 1 && vp[i] == global.to()[i])))
                continue;

            if (vp[i] == global.from()[i] || vp[i] == global.to()[i])
                return true;
        }

        return false;
    }

    const Box&      global;
    bool            wrap;
};

struct LocalOrGlobalBoundary
{
    typedef     TripletMergeTreeBlock::Box             Box;

                LocalOrGlobalBoundary(const Box& local_, const Box& global_, bool wrap):
                    local_test(local_), global_test(global_, wrap)                  {}

    bool        operator()(TripletMergeTreeBlock::Index v) const                    { return local_test(v) || global_test(v); }

    Box::BoundsTest         local_test;
    GlobalBoundary          global_test;
};

struct MergeSparsify
{
            MergeSparsify(bool wrap_):
                wrap(wrap_)                                                         {}

    void    operator()(void* b_, const diy::ReduceProxy& srp, const diy::RegularSwapPartners& partners) const
    {
        typedef              TripletMergeTreeBlock::TripletMergeTree    TripletMergeTree;
        typedef              TripletMergeTreeBlock::Box                 Box;
        typedef              TripletMergeTreeBlock::Vertex              Vertex;
        typedef              TripletMergeTreeBlock::Index               Index;

        LOG_SEV(debug) << "Entered merge_sparsify()";

        TripletMergeTreeBlock*             b        = static_cast<TripletMergeTreeBlock*>(b_);
        unsigned                    round    = srp.round();
        LOG_SEV(debug) << "Round: " << round;
        LOG_SEV_IF(srp.master()->communicator().rank() == 0, info) << "round = " << srp.round();

        // receive trees, merge, and sparsify
        int in_size = srp.in_link().size();
        LOG_SEV(debug) << "  incoming link size: " << in_size;
        if (in_size)
        {
            std::vector<Box>        bounds(in_size, b->global.grid_shape());
            std::vector<TripletMergeTree>  trees(in_size, b->mt.negate());
            for (int i = 0; i < in_size; ++i)
            {
              int nbr_gid = srp.in_link().target(i).gid;
              if (nbr_gid == srp.gid())
              {
                  bounds[i].swap(b->global);
                  trees[i].swap(b->mt);
                  LOG_SEV(debug) << "  swapped in tree of size: " << trees[i].size();
              } else
              {
                  srp.dequeue(nbr_gid, bounds[i]);
                  srp.dequeue(nbr_gid, trees[i]);
                  LOG_SEV(debug) << "  received tree of size: " << trees[i].size();
                  srp.incoming(nbr_gid).wipe();
              }
            }
            LOG_SEV(debug) << "  trees and bounds received";

            // merge boxes
            b->global.from() = bounds.front().from();
            b->global.to()   = bounds.back().to();
            LOG_SEV(debug) << "  boxes merged: " << b->global.from() << " - " << b->global.to() << " (" << b->global.grid_shape() << ')';

            // merge trees and move vertices
            record_stats("Merging trees:", "{} and {}", trees[0].size(), trees[1].size());

            int d = -1;
            for (int i = 0; i < Box::dimension(); ++i) if (bounds.front().from()[i] != bounds.back().from()[i]) d = i;

            Vertex u, v;
            u = bounds.front().from();
            u[d] = bounds.front().to()[d];
            v = bounds.back().to();
            v[d] = bounds.back().from()[d];
            Box edges_domain(b->global.grid_shape(), u, v);

            Box side = edges_domain.side(d, false);
            std::vector<std::tuple<Index, Index>> edges;
            r::VerticesIterator<Vertex> it = r::VerticesIterator<Vertex>::begin(side.from(), side.to()),
                                        end = r::VerticesIterator<Vertex>::end(side.from(), side.to());

            while (it != end)
            {
                Vertex u = *it;
                for (const Index& v : edges_domain.link(u))
                {
                    edges.push_back(std::make_tuple(b->global.position_to_vertex()(u), v));
                }
                ++it;
            }

            if (wrap)
            {
                for (int i = 0; i < Box::dimension(); ++i)
                {
                    if (bounds.front().from()[i] == 0 && bounds.back().to()[i] == b->global.grid_shape()[i] - 1 && bounds.front().to()[i] + 1 == bounds.back().from()[i])
                    {
                        Box side = bounds.back().side(i, true);
                        Vertex side_to = side.to();
                        for (int j = 0; j < Box::dimension(); ++j)
                        {
                            if (side.to()[j] == b->global.grid_shape()[j] - 1) side.to()[j]++;
                        }

                        it = r::VerticesIterator<Vertex>::begin(side.from(), side_to);
                        end = r::VerticesIterator<Vertex>::end(side.from(), side_to);
                        while (it != end)
                        {
                            Vertex u = *it;
                            if (bounds.back().contains(u)) {
                                for (const Index& v : side.link(u))
                                {
                                    if (bounds.front().contains(v)) edges.push_back(std::make_tuple(b->global.position_to_vertex()(u), v));
                                }
                            }
                            ++it;
                        }
                    }
                }
            }

            trees[0].swap(b->mt);
            r::merge(b->mt, trees[1], edges);

            trees.clear();
            LOG_SEV(debug) << "  trees merged: " << b->mt.size();
            record_stats("Trees merged:", "{}", b->mt.size());

            // sparsify
            r::sparsify(b->mt, LocalOrGlobalBoundary(b->local, b->global, wrap));
            record_stats("Trees sparsified:", "{}", b->mt.size());
        }

        // send (without the vertices) to the neighbors
        int out_size = srp.out_link().size();
        if (out_size == 0)        // final round: create the final local-global tree, nothing needs to be sent
        {
            LOG_SEV(debug) << "Sparsifying final tree of size: " << b->mt.size();
            r::sparsify(b->mt, b->local.bounds_test());
            record_stats("Final sparsified:", "{}", b->mt.size());
            LOG_SEV(debug) << "[" << b->gid << "] " << "Final tree size: " << b->mt.size();
            return;
        }

        TripletMergeTree mt_out(b->mt.negate());       // tree sparsified w.r.t. global boundary (dropping internal nodes)
        r::sparsify(mt_out, b->mt, GlobalBoundary(b->global, wrap));
        record_stats("Outgoing tree:", "{}", mt_out.size());

        for (int i = 0; i < out_size; ++i)
        {
          diy::BlockID nbr_bid = srp.out_link().target(i);
          if (nbr_bid.gid != srp.gid())
          {
            srp.enqueue(nbr_bid, b->global);
            srp.enqueue(nbr_bid, mt_out, &save_no_vertices);
          }
        }
    }

    bool                    wrap;
};

// debug only
void save_grids(void* b_, const diy::Master::ProxyWithLink& cp, void*)
{
    TripletMergeTreeBlock* b = static_cast<TripletMergeTreeBlock*>(b_);

    std::string outfn = fmt::format("block-debug-b{}.npy", b->gid);
    diy::mpi::io::file  f(cp.master()->communicator(), outfn, diy::mpi::io::file::create | diy::mpi::io::file::wronly);
    diy::io::NumPy writer(f);
    writer.write_header<Real,TripletMergeTreeBlock::Vertex>(b->grid.shape());
    diy::DiscreteBounds bounds;
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
    BOOST_FOREACH(TripletMergeTreeBlock::Index u, b->local.link(0))
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
    int         threads    = 1;

    std::string profile_path;
    std::string log_level = "info";

    Options ops(argc, argv);
    ops
        >> Option('b', "blocks",    nblocks,      "number of blocks to use")
        >> Option('m', "memory",    in_memory,    "maximum blocks to store in memory")
        >> Option('j', "jobs",      threads,      "threads to use during the computation")
        >> Option('s', "storage",   prefix,       "storage prefix")
        >> Option('p', "profile",   profile_path, "path to keep the execution profile")
        >> Option('l', "log",       log_level,    "log level")
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
                                       threads,
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

    diy::DiscreteBounds domain;
    domain.min[0] = domain.min[1] = domain.min[2] = 0;
    for (unsigned i = 0; i < reader.shape().size(); ++i)
      domain.max[i] = reader.shape()[i] - 1;

    diy::RegularDecomposer<diy::DiscreteBounds>::BoolVector       share_face(3, false);
    diy::RegularDecomposer<diy::DiscreteBounds>::BoolVector       wrap(3, wrap_);
    diy::RegularDecomposer<diy::DiscreteBounds>                   decomposer(3, domain, assigner.nblocks(), share_face, wrap);
    if (wrap_)
    {
        for (unsigned i = 0; i < 3; ++i)
            if (decomposer.divisions[i] < 2)
            {
                LOG_SEV(fatal) << "Can't have fewer than two divisions per side, when wrap is on";
                return -1;
            }
    }

    LoadAdd create(master, reader, negate, wrap_);
    decomposer.decompose(world.rank(), assigner, create);
    LOG_SEV_IF(world.rank() == 0, info) << "Domain decomposed: " << master.size();
    LOG_SEV_IF(world.rank() == 0, info) << "  (data read)";
    delete reader_ptr;

    world.barrier();
    LOG_SEV_IF(world.rank() == 0, info) << "Time to read data:       " << dlog::clock_to_string(timer.elapsed());
    timer.restart();

    // debug only
    //master.foreach(&save_grids);
    //master.foreach(&test_link);

    master.foreach(&compute_tree);

    world.barrier();
    LOG_SEV_IF(world.rank() == 0, info) << "Time to compute tree:    " << dlog::clock_to_string(timer.elapsed());
    timer.restart();

    // perform the global swap-reduce
    int k = 2;
    diy::RegularSwapPartners  partners(decomposer, k, true);
    diy::reduce(master, assigner, partners, MergeSparsify(wrap_));

    world.barrier();
    LOG_SEV_IF(world.rank() == 0, info) << "Time for the reduction:  " << dlog::clock_to_string(timer.elapsed());
    timer.restart();

    // save the result
    if (!split)
        diy::io::write_blocks(outfn, world, master);
    else
        diy::io::split::write_blocks(outfn, world, master);

    world.barrier();
    LOG_SEV_IF(world.rank() == 0, info) << "Time to output trees:    " << dlog::clock_to_string(timer.elapsed());
    timer.restart();

    dlog::prof.flush();     // TODO: this is necessary because the profile file will close before
                            //       the global dlog::prof goes out of scope and flushes the events.
                            //       Need to eventually fix this.
    dlog::stats.flush();
}