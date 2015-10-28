#include <iostream>
#include <string>

#include <opts/opts.h>
#include <dlog/log.h>
#include <dlog/stats.h>

#include <diy/decomposition.hpp>
#include <diy/io/block.hpp>

#include "format.h"

#include "merge-tree-block.h"

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

    std::string infn, outfn;
    if (  ops >> Present('h', "help", "show help message") ||
        !(ops >> PosOption(infn) >> PosOption(outfn)))
    {
        fmt::print("Usage: {} IN.lgt OUT.txt\n{}", argv[0], ops);
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
    diy::io::read_blocks(infn, world, assigner, master);
    LOG_SEV(info) << "Blocks read: " << master.size();

    std::ofstream out(outfn.c_str());
    master.foreach<MergeTreeBlock>([&](MergeTreeBlock* b, const diy::Master::ProxyWithLink& cp, void*)
                                   {
                                    const MergeTreeBlock::MergeTree& mt = b->mt;

                                    fmt::print(out, "Block {}", b->gid);
                                    if (verbose)
                                        fmt::print(out, " (local: {}, core: {})", b->local, b->core);
                                    fmt::print(out, "\n");
                                    BOOST_FOREACH(MergeTreeBlock::MergeTree::Neighbor n, mt.nodes() | reeber::ba::map_values)
                                    {
                                        fmt::print(out, "Node {}", n->vertex);
                                        if (verbose)
                                            fmt::print(out, " ({})", b->global.position(n->vertex));
                                        if (n->parent)
                                            fmt::print(out, " -> {}:\n", n->parent->vertex);
                                        else
                                            fmt::print(out, " (root):\n");
                                        if (!n->vertices.empty())
                                        {
                                            std::sort(n->vertices.begin(), n->vertices.end(),
                                                      [](const MergeTreeBlock::MergeTree::Node::ValueVertex& x,
                                                         const MergeTreeBlock::MergeTree::Node::ValueVertex& y)
                                                      { return x.second < y.second; });
                                            BOOST_FOREACH(const MergeTreeBlock::MergeTree::Node::ValueVertex& vv, n->vertices)
                                            {
                                                fmt::print(out, " {}", vv.second);
                                                if (verbose)
                                                    fmt::print(out, " ({})", b->global.position(vv.second));
                                            }
                                            fmt::print(out, "\n");
                                        }
                                    }
                                   });

    dlog::prof.flush();
}
