#include <iostream>
#include <string>

#include <boost/range/adaptor/map.hpp>

#include <dlog/stats.h>
#include <dlog/log.h>
#include <opts/opts.h>

#include <diy/io/numpy.hpp>

#include <reeber/merge-tree.h>
#include <reeber/grid.h>
#include <reeber/box.h>
namespace r = reeber;

#include "format.h"

typedef     REEBER_REAL                     Real;

typedef     r::Grid<Real, 3>                Grid;
typedef     Grid::Index                     Index;
typedef     Grid::Vertex                    Vertex;
typedef     Grid::Value                     Value;
typedef     r::Box<3>                       Box;
typedef     r::MergeTree<Index, Value>      MergeTree;

int main(int argc, char** argv)
{
    using namespace opts;

    std::string profile_path;
    std::string log_level = "info";
    Options ops(argc, argv);
    ops
        >> Option('p', "profile", profile_path, "path to keep the execution profile")
        >> Option('l', "log",     log_level,    "log level")
    ;
    bool        negate      = ops >> Present('n', "negate", "sweep superlevel sets");

    std::string infn;
    if (  ops >> Present('h', "help", "show help message") ||
        !(ops >> PosOption(infn)))
    {
        fmt::print("Usage: {} IN.npy\n{}", argv[0], ops);
        return 1;
    }

    dlog::add_stream(std::cerr, dlog::severity(log_level))
        << dlog::stamp() << dlog::color_pre() << dlog::level() << dlog::color_post() >> dlog::flush();

    // need some basic MPI to use for file IO
    diy::mpi::environment   env(argc, argv);
    diy::mpi::communicator  world;

    std::ofstream   profile_stream;
    if (profile_path.empty())
        dlog::prof.add_stream(std::cerr);
    else
    {
        std::string profile_fn = fmt::format("{}-r{}.prf", profile_path, world.rank());
        profile_stream.open(profile_fn.c_str());
        dlog::prof.add_stream(profile_stream);
    }

    // read the full infn
    diy::mpi::io::file      in(world, infn, diy::mpi::io::file::rdonly);
    diy::io::NumPy          reader(in);
    reader.read_header();

    if (reader.word_size() != sizeof(Grid::Value))
    {
        fmt::print("Incompatible value types: {} {}\n", reader.word_size(), sizeof(Grid::Value));
        return 1;
    }

    diy::DiscreteBounds box;
    for (unsigned i = 0; i < 3; ++i)
    {
        box.min[i] = 0;
        box.max[i] = reader.shape()[i] - 1;
    }
    Grid    g(Vertex(reader.shape()));
    reader.read(box, g.data());
    fmt::print("Grid shape: {}\n", g.shape());

    MergeTree mt(negate);

    Box domain(g.shape());
    r::compute_merge_tree(mt, domain, g);
    fmt::print("Tree constructed: {}\n", mt.size());

#ifdef COUNTERS
    fmt::print("{}\n", COUNTER(MergeTree::CollapseEvent));
    fmt::print("{}\n", COUNTER(MergeTree::EraseEvent));
    fmt::print("{}\n", COUNTER(MergeTree::FindStepEvent));
#endif

    size_t leaves = 0;
    BOOST_FOREACH(MergeTree::Neighbor n, ((const MergeTree&) mt).nodes() | boost::adaptors::map_values)
    {
        if (n->parent && n->children.size() == 0)
            ++leaves;
    }
    fmt::print("Leaves: {}\n", leaves);
}
