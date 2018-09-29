#include <iostream>
#include <string>

#include <opts/opts.h>
#include <dlog/log.h>
#include <dlog/stats.h>

#include <diy/decomposition.hpp>
#include <diy/io/block.hpp>

#include "format.h"

#include "triplet-merge-tree-block.h"
#include "output-persistence.h"

typedef diy::RegularDecomposer<diy::DiscreteBounds>                 Decomposer;

struct IsLocalTest
{
    using Position = TripletMergeTreeBlock::Box::Position;
    using Neighbor = TripletMergeTreeBlock::TripletMergeTree::Neighbor;

    const Decomposer& decomposer;

    IsLocalTest(const Decomposer& decomposer_) : decomposer(decomposer_) {}

    bool operator()(const TripletMergeTreeBlock& block, Neighbor from) const
    {
        Position from_position = block.global.position(from->vertex);
        return decomposer.lowest_gid(from_position) == block.gid;
    }
};


using OutputPairsR = OutputPairs<TripletMergeTreeBlock, IsLocalTest>;

int main(int argc, char** argv)
{
    diy::mpi::environment   env(argc, argv);
    diy::mpi::communicator  world;

    using namespace opts;

    std::string prefix     = "./DIY.XXXXXX";
    int         in_memory  = -1;
    int         threads    = 1;

    std::string profile_path;
    std::string log_level = "info";
    Real rho;
    Options ops(argc, argv);
    ops
        >> Option('m', "memory",    in_memory,    "maximum blocks to store in memory")
        >> Option('j', "jobs",      threads,      "threads to use during the computation")
        >> Option('t', "threshold", rho, "threshold")
        >> Option('s', "storage",   prefix,       "storage prefix")
        >> Option('p', "profile",   profile_path, "path to keep the execution profile")
        >> Option('l', "log",       log_level,    "log level")
    ;
    bool verbose = ops >> Present('v', "verbose", "verbose output");
    bool split   = ops >> Present(     "split",   "use split IO");

    std::string infn, outfn;
    if (  ops >> Present('h', "help", "show help message") ||
        !(ops >> PosOption(infn) >> PosOption(outfn)))
    {
        fmt::print("Usage: {} IN.lgt OUT.dgm\n{}", argv[0], ops);
        return 1;
    }

    dlog::add_stream(std::cerr, dlog::severity(log_level))
        << dlog::stamp() << dlog::aux_reporter(world.rank()) << dlog::color_pre() << dlog::level() << dlog::color_post() >> dlog::flush();

    std::ofstream   profile_stream;
    if (profile_path == "-")
        dlog::prof.add_stream(std::cerr);
    else if (!profile_path.empty())
    {
        std::string profile_fn = fmt::format("{}.prf", profile_path);
        profile_stream.open(profile_fn.c_str());
        dlog::prof.add_stream(profile_stream);
    }

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

    diy::ContiguousAssigner     assigner(world.size(), 0);

    // load the trees
    if (!split)
        diy::io::read_blocks(infn, world, assigner, master);
    else
        diy::io::split::read_blocks(infn, world, assigner, master);
    LOG_SEV_IF(world.rank() == 0, info) << "Blocks read: " << master.size();

    // get the domain bounds from any block that's in memory (they are all the same) and set up a decomposer
    TripletMergeTreeBlock::Box global = static_cast<TripletMergeTreeBlock*>(((const diy::Master&) master).block(master.loaded_block()))->global;
    diy::DiscreteBounds domain;
    for (unsigned i = 0; i < 3; ++i)
    {
        domain.min[i] = global.from()[i];
        domain.max[i] = global.to()[i];
    }

    // output persistence
    diy::RegularDecomposer<diy::DiscreteBounds>     decomposer(3, domain, assigner.nblocks());
    IsLocalTest test_local(decomposer);

    bool ignore_zero_persistence = false;
    OutputPairsR::ExtraInfo extra(outfn, verbose, world);
    master.foreach([&extra, &test_local, rho, ignore_zero_persistence] (TripletMergeTreeBlock* b, const diy::Master::ProxyWithLink& cp) {
        output_persistence(b, cp, extra, test_local, rho, ignore_zero_persistence); });

    dlog::prof.flush();
}
