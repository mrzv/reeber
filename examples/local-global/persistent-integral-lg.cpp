#include <iostream>
#include <string>

#include <boost/foreach.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/range/adaptor/map.hpp>

#include <opts/opts.h>
#include <dlog/log.h>
#include <dlog/stats.h>

#include <diy/decomposition.hpp>
#include <diy/reduce-operations.hpp>
#include <diy/io/block.hpp>

#include "format.h"

#include "merge-tree-block.h"
#include "persistent-integral-block.h"

namespace ba = boost::adaptors;

typedef diy::RegularDecomposer<diy::DiscreteBounds>                 Decomposer;
typedef MergeTreeBlock::Vertex                                      Vertex;
typedef MergeTreeBlock::Value                                       Value;

void compute_average(void *b_, const diy::Master::ProxyWithLink& cp, void*)
{
    const MergeTreeBlock& b     = *static_cast<MergeTreeBlock*>(b_);
    double                value = 0;
    size_t                count = 0;

    BOOST_FOREACH(const Neighbor node, const_cast<const MergeTree&>(b.mt).nodes() | ba::map_values)
    {
        if (b.core.contains(node->vertex))
        {
            value += node->value;
            count += 1;
        }
        BOOST_FOREACH(const MergeTreeNode::ValueVertex& x, node->vertices)
        {
            if (b.core.contains(x.second))
            {
                value += x.first;
                count += 1;
            }
        }
    }

    cp.all_reduce(value, std::plus<double>());
    cp.all_reduce(count, std::plus<size_t>());
}

class TreeTracer
{
    private:
        typedef std::map<MergeTreeNode::Vertex, MinIntegral>        MinIntegralMap;

    public:
                    TreeTracer(const Decomposer& decomposer_, diy::Master& pi_master_, Real m_, Real t_) :
                         decomposer(decomposer_), pi_master(pi_master_), m(m_), t(t_)
        {}

        void        operator()(void *b, const diy::ReduceProxy& rp) const
        {
            MergeTreeBlock &block = *static_cast<MergeTreeBlock*>(b);

            if (rp.in_link().size() == 0)
            {
                MergeTreeBlock::OffsetGrid::GridProxy gp(0, block.global.shape());
                trace(block.mt.find_root(), block, gp, rp);
            }
            else
            {
                MinIntegralMap mi_map;
                for (int i = 0; i < rp.in_link().size(); ++i)
                {
                    int gid = rp.in_link().target(i).gid;
                    assert(gid == i);

                    // ALTERNATIVE: diy::MemoryBiffer& incoming = rp.incoming(gid);
                    while (rp.incoming(gid)) // ALTERNATIVE: while (incoming)
                    {
                        // ALTERNATIVE: diy::load(incoming, mi);
                        MinIntegral mi;
                        rp.dequeue(gid, mi);
                        if (mi_map.find(mi.min_vtx) == mi_map.end())
                            mi_map[mi.min_vtx] = mi;
                        else
                            mi_map[mi.min_vtx].combine(mi);
                    }
                }

                PersistentIntegralBlock *pi_block = new PersistentIntegralBlock(block);
                for (MinIntegralMap::const_iterator it = mi_map.begin(); it != mi_map.end(); ++it)
                    pi_block->add_integral(it->second);
                pi_master.add(rp.gid(), pi_block, new diy::Link);
            }
        }

        void        trace(Neighbor n, const MergeTreeBlock& block, const MergeTreeBlock::OffsetGrid::GridProxy& gp, const diy::ReduceProxy& rp) const
        {
            BOOST_FOREACH(Neighbor child, n->children)
            {
                if (block.mt.cmp(t, child->value))                                   // still too high
                {
                    trace(child, block, gp, rp);
                }
                else                                                                 // crossing the threshold
                {
                    MinIntegral mi = integrate(child, block, gp);
                    if (block.mt.cmp(m, mi.min_val))                                 // the component is not deep enough
                        continue;

                    if (mi.integral == 0)                                            // non-local extrema are redundant
                        continue;

                    int dest_gid = decomposer.point_to_gid(block.global.position(mi.min_vtx));
                    diy::BlockID dest = rp.out_link().target(dest_gid);              // out_link targets are ordered as gids
                    assert(dest.gid == dest_gid);
                    rp.enqueue(dest, mi);
                }
            }
        }

        MinIntegral integrate(Neighbor n, const MergeTreeBlock& block, const MergeTreeBlock::OffsetGrid::GridProxy& gp) const
        {
            MinIntegral mi_res(n);

            // Find contribution from n
            if (block.core.contains(n->vertex))
            {
                mi_res.integral += n->value * block.cell_size[0] * block.cell_size[1] * block.cell_size[2];
                ++mi_res.n_cells;
                mi_res.push_back(MergeTreeNode::ValueVertex(n->value, n->vertex));
            }

            BOOST_FOREACH(const MergeTree::Node::ValueVertex& x, n->vertices)
                if (block.core.contains(x.second) && block.mt.cmp(x.first, t))
                {
                    mi_res.integral += x.first * block.cell_size[0] * block.cell_size[1] * block.cell_size[2];
                    ++mi_res.n_cells;
                    mi_res.push_back(x);
                }

            // Find contribution from children + min
            BOOST_FOREACH(Neighbor child, n->children)
            {
                MinIntegral mi = integrate(child, block, gp);
                if (block.mt.cmp(mi, mi_res))
                {
                    mi_res.min_val = mi.min_val;
                    mi_res.min_vtx = mi.min_vtx;
                }
                mi_res.combine(mi);
            }

            return mi_res;
        }

    private:
        const Decomposer&  decomposer;
        diy::Master&       pi_master;
        Real               m;
        Real               t;
};

bool vv_cmp(const MergeTreeNode::ValueVertex& a, const MergeTreeNode::ValueVertex& b)
{
    return a.second < b.second;
}

class OutputIntegrals {
    public:
                    OutputIntegrals(std::string outfn_, bool verbose_) : outfn(outfn_), verbose(verbose_) {}
       void         operator()(void *b, const diy::Master::ProxyWithLink& cp, void* aux) const
       {
           PersistentIntegralBlock&  block = *static_cast<PersistentIntegralBlock*>(b);

           std::string   dgm_fn = fmt::format("{}-b{}.comp", outfn, block.gid);
           std::ofstream ofs(dgm_fn.c_str());
 
           MergeTreeBlock::OffsetGrid::GridProxy gp(0, block.global.shape());
           BOOST_FOREACH(MinIntegral &mi, block.persistent_integrals)
           {
               Vertex v = block.local.position(mi.min_vtx);
               ofs << v[0] * block.cell_size[0] << " " << v[1] * block.cell_size[1] << " " << v[2] * block.cell_size[2] << " ";
               if (verbose)
                   ofs << v[0] << "x" << v[1] << "x" << v[2] << " (" << mi.min_vtx << ") ";
               ofs <<  mi.integral;
               if (verbose)
                   ofs << " " << mi.n_cells;
              ofs << std::endl;
#ifdef REEBER_PERSISTENT_INTEGRAL_TRACE_VTCS
               std::sort(mi.vertices.begin(), mi.vertices.end(), vv_cmp);
               for (std::vector< MergeTreeNode::ValueVertex >::const_iterator it = mi.vertices.begin(); it != mi.vertices.end(); ++it)
                   ofs << "   " << it->second << " (" << block.local.position(it->second) <<  ")" << std::endl;
#if 0
               // Consistency check for debugging
               for (size_t i =0; i < mi.vertices.size()-1; ++i)
                   if (mi.vertices[i].second == mi.vertices[i+1].second)
                       LOG_SEV(info) << "Duplicate vertex " << mi.vertices[i].second << " in component " << mi.min_vtx;
#endif
#endif
           }
       }

    private:
        std::string outfn;
        bool        verbose;
};

int main(int argc, char** argv)
{
    diy::mpi::environment   env(argc, argv);
    diy::mpi::communicator  world;

    using namespace opts;

    std::string prefix      = "./DIY.XXXXXX";
    int         in_memory   = -1;
    int         threads     = 1;
    int         k           = 2;
    Real        m           = 200;
    Real        t           = 82;

    std::string profile_path;
    std::string log_level = "info";

    Options ops(argc, argv);
    ops
        >> Option('m', "memory",    in_memory,    "maximum blocks to store in memory")
        >> Option('j', "jobs",      threads,      "threads to use during the computation")
        >> Option('k', "k",         k,            "use k-ary swap")
        >> Option('s', "storage",   prefix,       "storage prefix")
        >> Option('p', "profile",   profile_path, "path to keep the execution profile")
        >> Option('l', "log",       log_level,    "log level")
        >> Option('x', "max",       m,            "maximum threshold")
        >> Option('i', "iso",       t,            "isofind threshold")
    ;
    bool absolute = ops >> Present('a', "absolute", "use absolute values for thresholds (instead of multiples of mean");
    bool verbose  = ops >> Present('v', "verbose",  "verbose output: logical coordiantes and number of cells");

    std::string infn, outfn;
    if (  ops >> Present('h', "help", "show help message") ||
        !(ops >> PosOption(infn) >> PosOption(outfn)))
    {
        if (world.rank() == 0)
        {
            fmt::print("Usage: {} IN.lgt OUT.pi\n{}", argv[0], ops);
        }
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

    world.barrier();
    dlog::Timer timer;
    LOG_SEV_IF(world.rank() == 0, info) << "Starting computation";

    diy::FileStorage            storage(prefix);

    diy::Master                 mt_master(world,
                                          threads,
                                          in_memory,
                                          &MergeTreeBlock::create,
                                          &MergeTreeBlock::destroy,
                                          &storage,
                                          &MergeTreeBlock::save,
                                          &MergeTreeBlock::load);

    diy::Master                 pi_master(world,
                                          threads,
                                          in_memory);

    diy::ContiguousAssigner     assigner(world.size(), 0);

    // load the trees
    LOG_SEV_IF(world.rank() == 0, debug) << "Reading blocks from " << infn;
    diy::io::read_blocks(infn, world, assigner, mt_master);
    LOG_SEV_IF(world.rank() ==0 , info) << "Blocks read: " << mt_master.size();

    world.barrier();
    LOG_SEV_IF(world.rank() == 0, info) << "Time to read data:                    " << dlog::clock_to_string(timer.elapsed());
    timer.restart();

    // get the domain bounds from any block that's in memory (they are all the same) and set up a decomposer
    MergeTreeBlock::Box global = static_cast<MergeTreeBlock*>(((const diy::Master&) mt_master).block(mt_master.loaded_block()))->global;
    diy::DiscreteBounds domain;
    for (unsigned i = 0; i < 3; ++i)
    {
        domain.min[i] = global.from()[i];
        domain.max[i] = global.to()[i];
    }
    diy::RegularDecomposer<diy::DiscreteBounds>     decomposer(3, domain, assigner, Decomposer::BoolVector(3, true));

    // Compute average
    if (!absolute)
    {
        mt_master.foreach(&compute_average);
        mt_master.exchange();

        const diy::Master::ProxyWithLink& proxy = mt_master.proxy(mt_master.loaded_block());
        double mean = proxy.get<double>() / proxy.get<size_t>();
        m *= mean;
        t *= mean;

        LOG_SEV_IF(world.rank() == 0, info) << "Average value is " << mean << ". Using isofind threshold of " << t << " and maximum threshold of " << m;
    }

    // Compute and combine persistent integrals
    diy::all_to_all(mt_master, assigner, TreeTracer(decomposer, pi_master, m, t), k);

    world.barrier();
    LOG_SEV_IF(world.rank() == 0, info) << "Time to compute persistent integrals: " << dlog::clock_to_string(timer.elapsed());
    timer.restart();

    // Save persistent integrals to file
    pi_master.foreach(OutputIntegrals(outfn, verbose));

    world.barrier();
    LOG_SEV_IF(world.rank() == 0, info) << "Time to output persistent integrals:  " << dlog::clock_to_string(timer.elapsed());
    timer.restart();

    dlog::prof.flush();     // TODO: this is necessary because the profile file will close before
                            //       the global dlog::prof goes out of scope and flushes the events.
                            //       Need to eventually fix this.
    dlog::stats.flush();
}
