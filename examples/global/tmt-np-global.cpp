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

typedef     REEBER_REAL                       Real;

typedef     r::Grid<Real, 3>                  Grid;
typedef     r::OffsetGrid<Real, 3>            OffsetGrid;
typedef     Grid::Index                       Index;
typedef     Grid::Vertex                      Vertex;
typedef     Grid::Value                       Value;
typedef     r::Box<3>                         Box;
typedef     r::TripletMergeTree<Index, Value> TripletMergeTree;

struct OutputPairs
{
            OutputPairs(std::ostream& out_, TripletMergeTree &mt_):
                out(out_), negate(mt_.negate()), mt(mt_)           {}

    typedef       TripletMergeTree::Neighbor    Neighbor;

    void    operator()(const Neighbor from, const Neighbor through, const Neighbor to) const
    {
        if (from != to)
            fmt::print(out, "{} {} {} {} {} {}\n", from->vertex, from->value, through->vertex, through->value, to->vertex, to->value);
        else
            fmt::print(out, "{} {} {} --\n",    from->vertex, from->value, (negate ? "-inf" : "inf"));
    }

    std::ostream&       out;
    bool                negate;
    TripletMergeTree&          mt;
};

int main(int argc, char** argv)
{
    using namespace opts;

    std::string profile_path;
    std::string log_level = "info";
    int         jobs = r::task_scheduler_init::automatic;
    Options ops(argc, argv);
    ops
        >> Option('p', "profile", profile_path, "path to keep the execution profile")
        >> Option('l', "log",     log_level,    "log level")
        >> Option('j', "jobs",    jobs,         "number of threads to use (with TBB)")
    ;
    bool        negate      = ops >> Present('n', "negate", "sweep superlevel sets");

    std::string infn, outfn;
    if (  ops >> Present('h', "help", "show help message") ||
        !(ops >> PosOption(infn) >> PosOption(outfn)))
    {
        fmt::print("Usage: {} IN.npy OUT.mt\n{}", argv[0], ops);
        return 1;
    }

    r::task_scheduler_init init(jobs);

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
    v[0] /= 2;
    Vertex u = Vertex::zero();
    u[0] = v[0] + 1;
    Vertex w = Vertex::zero();
    w[0] = v[0];
    Vertex x = g.shape() - Vertex::one();
    x[0] = u[0];

    Box domain(g.shape()),
        domain1(g.shape(), Vertex::zero(), v),
        domain2(g.shape(), u, g.shape() - Vertex::one()),
        edges_domain(g.shape(), w, x);

    OffsetGrid g1(g.shape(), domain1.from(), domain1.to()),
               g2(g.shape(), domain2.from(), domain2.to());

    r::VerticesIterator<Vertex> it = r::VerticesIterator<Vertex>::begin(domain1.from(), domain1.to()),
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

    TripletMergeTree mt1(negate);
    TripletMergeTree mt2(negate);

#if 0
    r::compute_merge_tree(mt1, domain1, g1);
    r::compute_merge_tree(mt2, domain2, g2);

    std::vector<std::tuple<Index, Index, Index>> edges;
    it = r::VerticesIterator<Vertex>::begin(edges_domain.from(), edges_domain.to());
    end = r::VerticesIterator<Vertex>::end(edges_domain.from(), edges_domain.to());
    while (it != end)
    {
        Index u = domain.position_to_vertex()(*it);
        if (domain1.contains(u))
        {
            for (const Index& v : edges_domain.link(u))
            {
                if (domain2.contains(v))
                {
                    auto nu = mt1.node(u);
                    auto nv = mt2.node(v);
                    if (mt1.cmp(nv, nu)) edges.push_back(std::make_tuple(u, u, v));
                    else edges.push_back(std::make_tuple(v, v, u)); 
                }
            }
        }
        ++it;
    }

    it = r::VerticesIterator<Vertex>::begin(domain1.from(), domain1.to()),
    end = r::VerticesIterator<Vertex>::end(domain1.from(), domain1.to());
    dlog::Timer t;
    r::merge(mt1, mt2, edges);
    fmt::print("Time: {}\n", t.elapsed());
#else
    dlog::Timer t;
    r::compute_merge_tree2(mt1, domain, g);
    fmt::print("Time: {}\n", t.elapsed());
    fmt::print("Tree constructed: {}\n", mt1.size());
#endif

    std::ofstream ofs(outfn.c_str());
    r::traverse_persistence(mt1, OutputPairs(ofs, mt1));


    dlog::prof.flush();     // TODO: this is necessary because the profile file will close before
                            //       the global dlog::prof goes out of scope and flushes the events.
                            //       Need to eventually fix this.
}
