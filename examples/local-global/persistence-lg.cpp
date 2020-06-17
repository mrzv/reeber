#include <iostream>
#include <string>

#include <opts/opts.h>
#include <dlog/log.h>
#include <dlog/stats.h>

#include <diy/decomposition.hpp>
#include <diy/io/block.hpp>

#include <reeber/format.h>

#include "merge-tree-block.h"

typedef diy::RegularDecomposer<diy::DiscreteBounds>                 Decomposer;

struct OutputPairs
{
    struct ExtraInfo
    {
                            ExtraInfo(const std::string& outfn_, const Decomposer& decomposer_, bool verbose_):
                                outfn(outfn_), decomposer(decomposer_), verbose(verbose_)   {}
        std::string         outfn;
        const Decomposer&   decomposer;
        bool                verbose;
    };

    typedef MergeTreeBlock::MergeTree::Neighbor                         Neighbor;
    typedef MergeTreeBlock::Box                                         Box;

            OutputPairs(const MergeTreeBlock& block_, const ExtraInfo& extra_):
                block(block_),
                extra(extra_)
    {
        std::string   dgm_fn = fmt::format("{}-b{}.dgm", extra.outfn, block.gid);
        ofs.open(dgm_fn.c_str());
    }

    void    operator()(Neighbor from, Neighbor through, Neighbor to) const
    {
        Box::Position from_position = block.global.position(from->vertex);
        if (extra.decomposer.lowest_gid(from_position) != block.gid)
            return;

        if (extra.verbose)
        {
            if (from != to)
                fmt::print(ofs, "{} {} {} {} {} {}\n", from->vertex, from->value, through->vertex, through->value, to->vertex, to->value);
            else
                fmt::print(ofs, "{} {} {} --\n",       from->vertex,  from->value, (block.mt.negate() ? "-inf" : "inf"));
        } else
            fmt::print(ofs, "{} {}\n", from->value, through->value);
    }

    const MergeTreeBlock&       block;
    const ExtraInfo&            extra;
    mutable std::ofstream       ofs;
};

void output_persistence(MergeTreeBlock* b, const diy::Master::ProxyWithLink& cp, const OutputPairs::ExtraInfo& extra)
{
    LOG_SEV(debug) << "Block:   " << cp.gid();
    LOG_SEV(debug) << " Tree:   " << b->mt.size() << " with " << b->mt.count_roots() << " roots";
    LOG_SEV(debug) << " Local:  " << b->local.from()  << " - " << b->local.to();
    LOG_SEV(debug) << " Global: " << b->global.from() << " - " << b->global.to();

    r::traverse_persistence(b->mt, OutputPairs(*b, extra));
}

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
    Options ops(argc, argv);
    ops
        >> Option('m', "memory",    in_memory,    "maximum blocks to store in memory")
        >> Option('j', "jobs",      threads,      "threads to use during the computation")
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
                                       &MergeTreeBlock::create,
                                       &MergeTreeBlock::destroy,
                                       &storage,
                                       &MergeTreeBlock::save,
                                       &MergeTreeBlock::load);

    diy::ContiguousAssigner     assigner(world.size(), 0);

    // load the trees
    if (!split)
        diy::io::read_blocks(infn, world, assigner, master);
    else
        diy::io::split::read_blocks(infn, world, assigner, master);
    LOG_SEV_IF(world.rank() == 0, info) << "Blocks read: " << master.size();

    // get the domain bounds from any block that's in memory (they are all the same) and set up a decomposer
    MergeTreeBlock::Box global = static_cast<MergeTreeBlock*>(((const diy::Master&) master).block(master.loaded_block()))->global;
    diy::DiscreteBounds domain {3};
    for (unsigned i = 0; i < 3; ++i)
    {
        domain.min[i] = global.from()[i];
        domain.max[i] = global.to()[i];
    }
    diy::RegularDecomposer<diy::DiscreteBounds>     decomposer(3, domain, assigner.nblocks(), Decomposer::BoolVector(3, true));

    // output persistence
    OutputPairs::ExtraInfo extra(outfn, decomposer, verbose);
    master.foreach([&extra](MergeTreeBlock* b, const diy::Master::ProxyWithLink& cp) { output_persistence(b,cp,extra); });

    dlog::prof.flush();
}
