#include <iostream>

#include <opts/opts.h>
#include <dlog/log.h>
#include <dlog/stats.h>

#include <reeber/grid.h>
#include <reeber/merge-tree.h>
#include <reeber/merge-tree-serialization.h>
namespace r = reeber;

#include "format.h"

typedef     REEBER_REAL                     Real;
typedef     r::GridRef<void*, 3>            GridRef;
typedef     GridRef::Index                  Index;
typedef     GridRef::Vertex                 Vertex;
typedef     r::MergeTree<Index, Real>       MergeTree;
typedef     MergeTree::Neighbor             Neighbor;

struct OutputPairs
{
            OutputPairs(std::ostream& out_, bool negate_):
                out(out_), negate(negate_)              {}

    void    operator()(Neighbor from, Neighbor through, Neighbor to) const
    {
        if (from != to)
            fmt::print(out, "{} {} {} {} {} {}\n", from->vertex, from->value, through->vertex, through->value, to->vertex, to->value);
        else
            fmt::print(out, "{} {} {} --\n",    from->vertex,  from->value, (negate ? "-inf" : "inf"));
    }

    std::ostream&       out;
    bool                negate;
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

    std::string infn, outfn;
    if (  ops >> Present('h', "help", "show help message") ||
        !(ops >> PosOption(infn) >> PosOption(outfn)))
    {
        fmt::print("Usage: {} IN.mt OUT.dgm\n{}", argv[0], ops);
        return 1;
    }

    dlog::add_stream(std::cerr, dlog::severity(log_level))
        << dlog::stamp() << dlog::color_pre() << dlog::level() << dlog::color_post() >> dlog::flush();

    std::ofstream   profile_stream;
    if (profile_path.empty())
        dlog::prof.add_stream(std::cerr);
    else
    {
        std::string profile_fn = fmt::format("{}.prf", profile_path);
        profile_stream.open(profile_fn.c_str());
        dlog::prof.add_stream(profile_stream);
    }

    MergeTree mt; Vertex shape;

    diy::MemoryBuffer bb; bb.read(infn);
    diy::load(bb, shape);
    diy::load(bb, mt);
    fmt::print("Tree loaded: {}\n", mt.size());

    std::ofstream ofs(outfn.c_str());
    r::traverse_persistence(mt, OutputPairs(ofs, mt.negate()));
}
