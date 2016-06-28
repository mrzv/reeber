#include <iostream>
#include <string>

#include <boost/range/adaptor/map.hpp>

#include <dlog/stats.h>
#include <dlog/log.h>
#include <opts/opts.h>

#include <diy/io/numpy.hpp>

#include <reeber/triplet-merge-tree.h>
#include <reeber/grid.h>
#include <reeber/box.h>
#include <reeber/merge-tree-serialization.h>
namespace r = reeber;

#include "format.h"

typedef     REEBER_REAL                     Real;

typedef     r::Grid<Real, 3>                Grid;
typedef     r::OffsetGrid<Real, 3>          OffsetGrid;
typedef     Grid::Index                     Index;
typedef     Grid::Vertex                    Vertex;
typedef     Grid::Value                     Value;
typedef     r::Box<3>                       Box;
typedef     r::MergeTree<Index, Value>      MergeTree;

struct OutputPairs
{
            OutputPairs(std::ostream& out_, MergeTree &mt_):
                out(out_), negate(mt_.negate()), mt(mt_)           {}

    void    operator()(const Index from, const Index through, const Index to) const
    {
        if (from != to)
            fmt::print(out, "{} {} {} {} {} {}\n", from, mt.node(from)->value, through, mt.node(through)->value, to, mt.node(to)->value);
            // fmt::print(out, "{} {} {} {} {} {}\n", from, mt.value(from), through, mt.value(through), to, mt.value(to));
        else
            fmt::print(out, "{} {} {} --\n",    from, mt.node(from)->value, (negate ? "-inf" : "inf"));
            // fmt::print(out, "{} {} {} --\n",    from, mt.value(from), (negate ? "-inf" : "inf"));
    }

    std::ostream&       out;
    bool                negate;
    MergeTree&          mt;
};

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

    std::string infn, outfn;
    if (  ops >> Present('h', "help", "show help message") ||
        !(ops >> PosOption(infn) >> PosOption(outfn)))
    {
        fmt::print("Usage: {} IN.npy OUT.mt\n{}", argv[0], ops);
        return 1;
    }

    dlog::add_stream(std::cerr, dlog::severity(log_level))
        << dlog::stamp() << dlog::color_pre() << dlog::level() << dlog::color_post() >> dlog::flush();

    // need some basic MPI to use for file IO
    diy::mpi::environment   env(argc, argv);
    diy::mpi::communicator  world;

    std::ofstream   profile_stream;
    if (profile_path == "-")
        dlog::prof.add_stream(std::cerr);
    else if (!profile_path.empty())
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

    Vertex v = g.shape() - Vertex::one();
    v[2] /= 2;
    Vertex u = Vertex::zero();
    u[2] = v[2] + 1;
    Vertex w = Vertex::zero();
    w[2] = v[2];
    Vertex x = g.shape() - Vertex::one();
    x[2] = u[2];

    Box domain(g.shape()),
        domain1(g.shape(), Vertex::zero(), v),
        domain2(g.shape(), u, g.shape() - Vertex::one()),
        edges_domain(g.shape(), w, x);

    std::vector<std::tuple<Index, Index>> edges;
    r::VerticesIterator<Vertex> it = r::VerticesIterator<Vertex>::begin(edges_domain.from(), edges_domain.to()),
                                end = r::VerticesIterator<Vertex>::end(edges_domain.from(), edges_domain.to());
    while (it != end)
    {
        Vertex u = *it;
        if (domain1.contains(u))
        {
            for (const Index& v : edges_domain.link(u))
            {
                if (domain2.contains(v)) edges.push_back(std::make_tuple(domain.position_to_vertex()(u), v));
            }
        }
        ++it;
    }

    OffsetGrid g1(g.shape(), domain1.from(), domain1.to()),
               g2(g.shape(), domain2.from(), domain2.to());

    it = r::VerticesIterator<Vertex>::begin(domain1.from(), domain1.to());
    end = r::VerticesIterator<Vertex>::end(domain1.from(), domain1.to());
    while (it != end)
    {
        g1(*it) = g(*it);
        ++it;
    }
    it = r::VerticesIterator<Vertex>::begin(domain2.from(), domain2.to());
    end = r::VerticesIterator<Vertex>::end(domain2.from(), domain2.to());
    while (it != end)
    {
        g2(*it) = g(*it);
        ++it;
    }

    MergeTree mt1(negate);
    // MergeTree mt2(negate);

    // r::compute_merge_tree(mt1, domain1, g1);
    // r::compute_merge_tree(mt2, domain2, g2);
    // dlog::Timer t;
    // r::merge(mt1, mt2, edges);
    // fmt::print("Time: {}\n", t.elapsed());

    dlog::Timer t;
    r::compute_merge_tree2(mt1, domain, g);
    fmt::print("Time: {}\n", t.elapsed());
    fmt::print("Tree constructed: {}\n", mt1.size());

    std::ofstream ofs(outfn.c_str());
    r::traverse_persistence(mt1, OutputPairs(ofs, mt1));


    dlog::prof.flush();     // TODO: this is necessary because the profile file will close before
                            //       the global dlog::prof goes out of scope and flushes the events.
                            //       Need to eventually fix this.
}
