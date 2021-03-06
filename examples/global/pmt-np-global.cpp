#include <iostream>
#include <string>

#include <dlog/stats.h>
#include <dlog/log.h>
#include <opts/opts.h>

#include <diy/mpi.hpp>
#include <diy/io/numpy.hpp>

#include <reeber/path-merge-tree.h>
#include <reeber/merge-tree.h>
#include <reeber/grid.h>
#include <reeber/box.h>
#include <reeber/path-merge-tree-serialization.h>
namespace r = reeber;

#include <reeber/format.h>

#include "reader-interfaces.h"

typedef     REEBER_REAL                       Real;

typedef     r::Grid<Real, 3>                  Grid;
typedef     r::OffsetGrid<Real, 3>            OffsetGrid;
typedef     Grid::Index                       Index;
typedef     Grid::Vertex                      Vertex;
typedef     Grid::Value                       Value;
//typedef     r::Box<3>                         Box;
typedef     r::PathMergeTree<Index, Value>    PathMergeTree;
typedef     r::MergeTree<Index, Value>        MergeTree;

struct OutputPairs
{
            OutputPairs(std::ostream& out_, MergeTree &mt_):
                out(out_), negate(mt_.negate()), mt(mt_)           {}

    typedef       MergeTree::Neighbor    Neighbor;

    void    operator()(const Neighbor from, const Neighbor through, const Neighbor to) const
    {
        if (from != to)
            fmt::print(out, "{} {} {} {} {} {}\n", from->vertex, from->value, through->vertex, through->value, to->vertex, to->value);
        else
            fmt::print(out, "{} {} {} --\n",    from->vertex, from->value, (negate ? "-inf" : "inf"));
    }

    std::ostream&       out;
    bool                negate;
    MergeTree&          mt;
};

int main(int argc, char** argv)
{
    // need some basic MPI to use for file IO
    diy::mpi::environment   env(argc, argv);
    diy::mpi::communicator  world;

#ifdef REEBER_USE_BOXLIB_READER
    reeber::io::BoxLib::environment boxlib_env(argc, argv, world);
#endif

    using namespace opts;

    std::string profile_path;
    std::string log_level = "info";
    int         jobs = r::task_scheduler_init::automatic;
    int         d = 1;
    std::string tree_fn;

    Options ops(argc, argv);
    ops
        >> Option('p', "profile", profile_path, "path to keep the execution profile")
        >> Option('l', "log",     log_level,    "log level")
        >> Option('j', "jobs",    jobs,         "number of threads to use (with TBB)")
        >> Option('d', "scale",   d,            "downsampling factor")
        >> Option('t', "tree",    tree_fn,      "file to save the tree");
    ;
    bool        negate      = ops >> Present('n', "negate", "sweep superlevel sets");
    bool        split       = ops >> Present('s', "split",  "split domain and merge");

    std::string infn, outfn;
    if (  ops >> Present('h', "help", "show help message") ||
        !(ops >> PosOption(infn) >> PosOption(outfn)))
    {
        fmt::print("Usage: {} IN.npy OUT.dgm\n{}", argv[0], ops);
        return 1;
    }

    fmt::print(std::cerr, "File: {}\nNumber of threads: {}\n", infn, jobs);

    r::task_scheduler_init init(jobs);

    dlog::add_stream(std::cerr, dlog::severity(log_level))
        << dlog::stamp() << dlog::color_pre() << dlog::level() << dlog::color_post() >> dlog::flush();

    std::ofstream   profile_stream;
    if (profile_path == "-")
        dlog::prof.add_stream(std::cerr);
    else if (!profile_path.empty())
    {
        std::string profile_fn = fmt::format("{}-r{}.prf", profile_path, world.rank());
        profile_stream.open(profile_fn.c_str());
        dlog::prof.add_stream(profile_stream);
    }


    // set up the reader
    Reader* reader_ptr = Reader::create(infn, world);
    Reader& reader  = *reader_ptr;

    diy::DiscreteBounds box {3};
    for (unsigned i = 0; i < 3; ++i)
    {
        box.min[i] = 0;
        box.max[i] = reader.shape()[i] - 1;
    }
    Grid    g_(Vertex(reader.shape()));
    reader.read(box, g_.data());
    fmt::print(std::cerr, "Grid shape: {}\n", g_.shape());

    delete reader_ptr;

    Vertex shape = g_.shape();
    shape /= d;
    r::Box<3> domain(shape);
    Grid g(shape);
    r::VerticesIterator<Vertex> it = r::VerticesIterator<Vertex>::begin(domain.from(), domain.to()),
                                end = r::VerticesIterator<Vertex>::end(domain.from(), domain.to());
    while (it != end)
    {
        Vertex a = *it;
        a *= d;
        g(*it) = g_(a);
        ++it;
    }

    Grid().swap(g_);

    fmt::print(std::cerr, "Downsampled grid shape: {}\n", g.shape());

    Vertex v = g.shape() - Vertex::one();
    v[0] /= 2;
    Vertex u = Vertex::zero();
    u[0] = v[0] + 1;
    Vertex w = Vertex::zero();
    w[0] = v[0];
    Vertex x = g.shape() - Vertex::one();
    x[0] = u[0];

    r::Box<3>
        domain1(g.shape(), Vertex::zero(), v),
        domain2(g.shape(), u, g.shape() - Vertex::one()),
        edges_domain(g.shape(), w, x);

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

    PathMergeTree mt1(negate);
    PathMergeTree mt2(negate);

    if (split)
    {
        r::compute_merge_tree2(mt1, domain1, g1);
        r::compute_merge_tree2(mt2, domain2, g2);

        std::vector<std::tuple<Index, Index>> edges;
        it = r::VerticesIterator<Vertex>::begin(edges_domain.from(), edges_domain.to());
        end = r::VerticesIterator<Vertex>::end(edges_domain.from(), edges_domain.to());
        while (it != end)
        {
            Index u = domain.position_to_vertex()(*it);
            if (domain1.contains(u))
            {
                for (const Index& v : edges_domain.link(u))
                {
                    if (domain2.contains(v)) edges.push_back(std::make_tuple(u, v));
                }
            }
            ++it;
        }

        it = r::VerticesIterator<Vertex>::begin(domain1.from(), domain1.to()),
        end = r::VerticesIterator<Vertex>::end(domain1.from(), domain1.to());
        dlog::Timer t;
        r::merge(mt1, mt2, edges);
        dlog::Timer::duration elapsed = t.elapsed();
        fmt::print(std::cerr, "Time to merge: {}\n", elapsed);
        fmt::print("pmt-merge {} {}\n", jobs, elapsed);
    }
    else
    {
        dlog::Timer t;
        r::compute_merge_tree2(mt1, domain, g);
        dlog::Timer::duration elapsed = t.elapsed();
        fmt::print(std::cerr, "Time for compute_merge_tree2: {}\n", elapsed);
        fmt::print("pmt {} {}\n", jobs, elapsed);
    }

    if (outfn != "none")
    {
        std::ofstream ofs(outfn.c_str());

        // Lazy way out: copy mt1 into ordinary merge tree and traverse persistence there
        size_t roots = 0;
        MergeTree mt(negate);
        const auto& mt1_ = mt1;
        for (auto& x : mt1_.nodes())
        {
            auto n  = x.second;
            if (x.first != n->vertex) continue;
            PathMergeTree::Neighbor np = n->parent;

            auto on = mt.find_or_add(n->vertex, n->value);
            auto op = mt.find_or_add(np->vertex, np->value);
            if (on != op)       // mt uses parent == 0 to signify root
            {
                on->parent = op;
                op->children.push_back(on);
            } else
                ++roots;
        }
        fmt::print(std::cerr, "Found {} roots\n", roots);
        fmt::print(std::cerr, "Transferred {} roots\n", mt.count_roots());
        r::remove_degree2(mt, boost::lambda::constant(false));
        r::traverse_persistence(mt, OutputPairs(ofs, mt));
    }

    if (!tree_fn.empty())
    {
        diy::MemoryBuffer bb;
        diy::save(bb, mt1);
        bb.write(tree_fn);
    }

    dlog::prof.flush();     // TODO: this is necessary because the profile file will close before
                            //       the global dlog::prof goes out of scope and flushes the events.
                            //       Need to eventually fix this.
}
