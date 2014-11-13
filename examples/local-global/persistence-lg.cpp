#include <iostream>

#include <opts/opts.h>
#include <dlog/log.h>
#include <dlog/stats.h>

#include <diy/io/block.hpp>

#include "format.h"

#include "merge-tree-block.h"

struct OutputPairs
{
    typedef MergeTreeBlock::MergeTree::Neighbor                         Neighbor;
    typedef MergeTreeBlock::Box                                         Box;

            OutputPairs(std::ostream& out_, const Box& local_, bool negate_):
                out(out_), local(local_), negate(negate_)               {}

    void    operator()(Neighbor from, Neighbor through, Neighbor to) const
    {
        // FIXME: need to check for the overlap; in principle should use the decomposer, although that seems heavy
        if (!local.bounds_test()(from->vertex))
            return;

        if (from != to)
            fmt::print(out, "{} {} {} {} {} {}\n", from->vertex, from->value, through->vertex, through->value, to->vertex, to->value);
        else
            fmt::print(out, "{} {} {} --\n",    from->vertex,  from->value, (negate ? "-inf" : "inf"));
    }

    std::ostream&       out;
    const Box&          local;
    bool                negate;
};

void output_persistence(void* b_, const diy::Master::ProxyWithLink& cp, void* ofs_)
{
    MergeTreeBlock*     b       = static_cast<MergeTreeBlock*>(b_);
    std::ofstream&      ofs     = *static_cast<std::ofstream*>(ofs_);

    LOG_SEV(debug) << "Block:   " << cp.gid();
    LOG_SEV(debug) << " Tree:   " << b->mt.size() << " with " << b->mt.count_roots() << " roots";
    LOG_SEV(debug) << " Local:  " << b->local.from()  << " - " << b->local.to();
    LOG_SEV(debug) << " Global: " << b->global.from() << " - " << b->global.to();

    r::traverse_persistence(b->mt, OutputPairs(ofs, b->local, b->mt.negate()));
}

int main(int argc, char** argv)
{
    diy::mpi::environment   env(argc, argv);
    diy::mpi::communicator  world;

    using namespace opts;

    std::string prefix     = "./DIY.XXXXXX";
    int         in_memory  = -1;
    int         threads    = -1;

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

    LOG_SEV(info) << "Starting computation";
    diy::FileStorage            storage(prefix);

    diy::Communicator           comm(world);
    diy::Master                 master(comm,
                                       &MergeTreeBlock::create,
                                       &MergeTreeBlock::destroy,
                                       in_memory,
                                       threads,
                                       &storage,
                                       &MergeTreeBlock::save,
                                       &MergeTreeBlock::load);

    diy::ContiguousAssigner     assigner(world.size(), 0);

    // load the trees
    diy::io::read_blocks(infn, world, assigner, master);
    LOG_SEV(info) << "Blocks read: " << master.size();

    // output persistence
    std::ofstream ofs(outfn.c_str());
    master.foreach(&output_persistence, &ofs);

    dlog::prof.flush();
}
