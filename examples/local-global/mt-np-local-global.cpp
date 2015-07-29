#include <iostream>
#include <string>
#include <array>

#include <dlog/stats.h>
#include <dlog/log.h>
#include <opts/opts.h>

#include <diy/master.hpp>
#include <diy/assigner.hpp>
#include <diy/storage.hpp>
#include <diy/decomposition.hpp>
#include <diy/reduce.hpp>
#include <diy/partners/swap.hpp>
#include <diy/io/numpy.hpp>
#include <diy/io/block.hpp>

#ifdef REEBER_USE_BOXLIB_READER
#include <algorithm>
#include <reeber/io/boxlib.h>
#endif

#include "format.h"

#include "merge-tree-block.h"

#ifdef REEBER_USE_BOXLIB_READER
typedef reeber::io::BoxLib::reader Reader;
#else
typedef diy::io::NumPy Reader;
#endif

// Load the specified chunk of data, compute local merge tree, add block to diy::Master
struct LoadComputeAdd
{
    typedef         MergeTreeBlock::Vertex                          Vertex;
    typedef         MergeTreeBlock::MergeTree                       MergeTree;
    typedef         MergeTreeBlock::Box                             Box;
    typedef         MergeTreeBlock::OffsetGrid                      OffsetGrid;

                    LoadComputeAdd(diy::Master& master_, const Reader& reader_, bool negate_):
                        master(&master_), reader(reader_), negate(negate_)      {}

    void            operator()(int                          gid,
                               const diy::DiscreteBounds&   core,
                               const diy::DiscreteBounds&   bounds,
                               const diy::DiscreteBounds&   domain,
                               const diy::RegularGridLink&  link) const
    {
        MergeTreeBlock*         b  = new MergeTreeBlock;
        diy::RegularGridLink*   l  = new diy::RegularGridLink(link);

        Vertex                  full_shape = Vertex(domain.max) - Vertex(domain.min) + Vertex::one();

        OffsetGrid g(full_shape, bounds.min, bounds.max);
        reader.read(bounds, g.data(), true);      // collective; implicitly assumes same number of blocks on every processor

        b->gid = gid;
        b->mt.set_negate(negate);
        b->local = b->global = Box(full_shape, bounds.min, bounds.max);
        LOG_SEV(debug) << "Local box:  " << b->local.from()  << " - " << b->local.to();
        LOG_SEV(debug) << "Global box: " << b->global.from() << " - " << b->global.to();

        r::compute_merge_tree(b->mt, b->local, g,
                          [b,&g](MergeTreeBlock::Index v)
                          {
                            if (b->local.internal_test()(v))
                                return true;

                            MergeTreeBlock::Vertex  p    = b->local.position(v);
                            std::array<int, 3>      side = {{ 0, 0, 0 }};
                            for (int i = 0; i < 3; ++i)
                                if (p[i] == b->local.from()[i])
                                    side[i] = -1;
                                else if (p[i] == b->local.to()[i])
                                    side[i] = 1;

                            int zeroes = 0;
                            MergeTreeBlock::Box     side_box = b->local;
                            for (int i = 0; i < 3; ++i)
                            {
                                if (side[i] == -1)
                                    side_box.to()[i] = side_box.from()[i];
                                else if (side[i] == 1)
                                    side_box.from()[i] = side_box.to()[i];
                                else // (side[i] == 0)
                                    ++zeroes;
                            }
                            if (zeroes < 2)     // corner
                                return false;

                            typedef     MergeTreeBlock::MergeTree::Node::ValueVertex    ValueVertex;
                            ValueVertex vval = { g(v), v };
                            BOOST_FOREACH(MergeTreeBlock::Index u, side_box.link(v))
                            {
                                ValueVertex uval = { g(u), u };
                                if (b->mt.cmp(uval, vval))      // v is not a minimum
                                    return true;
                            }

                            return false;
                          });
        AssertMsg(b->mt.count_roots() == 1, "The tree can have only one root, not " << b->mt.count_roots());

        LOG_SEV(info) << "Initial tree size (" << b->gid << "): " << b->mt.size();

        int lid   = master->add(gid, b, l);
        static_cast<void>(lid);     // shut up the compiler about lid
    }

    diy::Master*            master;
    const Reader&           reader;
    bool                    negate;
};

void save_no_vertices(diy::BinaryBuffer& bb, const MergeTreeBlock::MergeTree& mt)
{
    reeber::Serialization<MergeTreeBlock::MergeTree>::save(bb, mt, false);
}

struct GlobalBoundary
{
    typedef     MergeTreeBlock::Box             Box;

                GlobalBoundary(const Box& global_):
                    global_test(global_)                        {}

    bool        operator()(MergeTreeBlock::Index v) const                           { return global_test(v); }

    Box::BoundaryTest       global_test;
};

struct LocalOrGlobalBoundary
{
    typedef     MergeTreeBlock::Box             Box;

                LocalOrGlobalBoundary(const Box& local_, const Box& global_):
                    local_test(local_), global_test(global_)                        {}

    bool        operator()(MergeTreeBlock::Index v) const                           { return local_test(v) || global_test(v); }

    Box::BoundsTest         local_test;
    Box::BoundaryTest       global_test;
};

void merge_sparsify(void* b_, const diy::ReduceProxy& srp, const diy::RegularSwapPartners& partners)
{
    typedef                     MergeTreeBlock::MergeTree           MergeTree;
    typedef                     MergeTreeBlock::Box                 Box;
    typedef                     MergeTree::Neighbor                 Neighbor;

    LOG_SEV(debug) << "Entered merge_sparsify()";

    MergeTreeBlock*             b        = static_cast<MergeTreeBlock*>(b_);
    unsigned                    round    = srp.round();
    LOG_SEV(debug) << "Round: " << round;

    // receive trees, merge, and sparsify
    int in_size = srp.in_link().size();
    LOG_SEV(debug) << "  incoming link size: " << in_size;
    if (in_size)
    {
        std::vector<Box>        bounds(in_size, b->global.grid_shape());
        std::vector<MergeTree>  trees(in_size, b->mt.negate());
        int in_pos = -1;
        for (int i = 0; i < in_size; ++i)
        {
          int nbr_gid = srp.in_link().target(i).gid;
          if (nbr_gid == srp.gid())
          {
              in_pos = i;
              bounds[i].swap(b->global);
              trees[i].swap(b->mt);
              LOG_SEV(debug) << "  swapped in tree of size: " << trees[i].size();
          } else
          {
              srp.dequeue(nbr_gid, bounds[i]);
              srp.dequeue(nbr_gid, trees[i]);
              LOG_SEV(debug) << "  received tree of size: " << trees[i].size();
          }
        }
        LOG_SEV(debug) << "  trees and bounds received";

        // merge boxes
        b->global.from() = bounds.front().from();
        b->global.to()   = bounds.back().to();
        LOG_SEV(debug) << "  boxes merged: " << b->global.from() << " - " << b->global.to() << " (" << b->global.grid_shape() << ')';

        // merge trees and move vertices
        r::merge(b->mt, trees);
        BOOST_FOREACH(Neighbor n, static_cast<const MergeTree&>(trees[in_pos]).nodes() | r::ba::map_values)
            if (!n->vertices.empty())
                b->mt[n->vertex]->vertices.swap(n->vertices);
        trees.clear();
        LOG_SEV(debug) << "  trees merged: " << b->mt.size();

        // sparsify
        sparsify(b->mt, LocalOrGlobalBoundary(b->local, b->global));
        remove_degree2(b->mt, b->local.bounds_test(), GlobalBoundary(b->global));
    }

    // send (without the vertices) to the neighbors
    int out_size = srp.out_link().size();
    if (out_size == 0)        // final round: create the final local-global tree, nothing needs to be sent
    {
        //LOG_SEV(info) << "Sparsifying final tree of size: " << b->mt.size();
        sparsify(b->mt, b->local.bounds_test());
        remove_degree2(b->mt, b->local.bounds_test());
        LOG_SEV(info) << "Final tree size: " << b->mt.size();
        return;
    }

    MergeTree mt_out(b->mt.negate());       // tree sparsified w.r.t. global boundary (dropping internal nodes)
    sparsify(mt_out, b->mt, b->global.boundary_test());

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

    std::string infn, outfn;
    if (  ops >> Present('h', "help", "show help message") ||
        !(ops >> PosOption(infn) >> PosOption(outfn)))
    {
        if (world.rank() == 0)
            fmt::print("Usage: {} IN.npy OUT.mt\n{}", argv[0], ops);
        return 1;
    }

    dlog::add_stream(std::cerr, dlog::severity(log_level))
        << dlog::stamp() << dlog::aux_reporter(world.rank()) << dlog::color_pre() << dlog::level() << dlog::color_post() >> dlog::flush();

#ifdef PROFILE
    if (threads != 1)
    {
        LOG_SEV_IF(world.rank() == 0, fatal) << "Cannot use profiling with more than one thread";
        return 1;
    }
#endif

    std::ofstream   profile_stream;
    if (profile_path == "-")
        dlog::prof.add_stream(std::cerr);
    else if (!profile_path.empty())
    {
        std::string profile_fn = fmt::format("{}-r{}.prf", profile_path, world.rank());
        profile_stream.open(profile_fn.c_str());
        dlog::prof.add_stream(profile_stream);
    }

    LOG_SEV(info) << "Starting computation";
    diy::FileStorage            storage(prefix);

    diy::Master                 master(world,
                                       threads,
                                       in_memory,
                                       &MergeTreeBlock::create,
                                       &MergeTreeBlock::destroy,
                                       &storage,
                                       &MergeTreeBlock::save,
                                       &MergeTreeBlock::load);

    diy::ContiguousAssigner     assigner(world.size(), nblocks);

#ifdef REEBER_USE_BOXLIB_READER
    Reader                      reader(infn, world);
#else
    diy::mpi::io::file          in(world, infn, diy::mpi::io::file::rdonly);
    Reader                      reader(in);
    reader.read_header();
#endif

    diy::DiscreteBounds domain;
    domain.min[0] = domain.min[1] = domain.min[2] = 0;
    for (unsigned i = 0; i < reader.shape().size(); ++i)
      domain.max[i] = reader.shape()[i] - 1;

    diy::RegularDecomposer<diy::DiscreteBounds>::BoolVector       share_face(3, true);

    LoadComputeAdd create(master, reader, negate);
    diy::decompose(3, world.rank(), domain, assigner, create, share_face);
    LOG_SEV(info) << "Domain decomposed: " << master.size();
    LOG_SEV(info) << "  (data read + local trees computed)";

    // perform the global swap-reduce
    int k = 2;
    diy::RegularSwapPartners  partners(3, nblocks, k, true);
    diy::reduce(master, assigner, partners, &merge_sparsify);

    // save the result
    diy::io::write_blocks(outfn, world, master);

    dlog::prof.flush();     // TODO: this is necessary because the profile file will close before
                            //       the global dlog::prof goes out of scope and flushes the events.
                            //       Need to eventually fix this.
}
