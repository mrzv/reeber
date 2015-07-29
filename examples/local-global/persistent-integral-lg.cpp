#include <iostream>
#include <string>

#include <opts/opts.h>
#include <dlog/log.h>
#include <dlog/stats.h>

#include <diy/decomposition.hpp>
#include <diy/reduce-operations.hpp>
#include <diy/io/block.hpp>

#include "format.h"

#include "merge-tree-block.h"

//#define REEBER_PERSISTENT_INTEGRAL_TRACE_VTCS

typedef diy::RegularDecomposer<diy::DiscreteBounds>                 Decomposer;
typedef MergeTreeBlock::MergeTree                                   MergeTree;
typedef MergeTreeBlock::Vertex                                      Vertex;
typedef MergeTreeBlock::Value                                       Value;
typedef MergeTreeBlock::MergeTree::Neighbor                         Neighbor;
typedef MergeTree::Node                                             MergeTreeNode;

struct MinIntegral
{
                     MinIntegral() : min_vtx(0), min_val(0), integral(0), n_cells(0)  {}
                     MinIntegral(const Neighbor min_node_, Real integral_ = 0, size_t n_cells_ = 0) :
                         min_vtx(min_node_->vertex), min_val(min_node_->value), integral(integral_), n_cells(n_cells_)  {}

    void             combine(const MinIntegral& other)              { integral += other.integral; n_cells += other.n_cells; append(other); }
#ifdef REEBER_PERSISTENT_INTEGRAL_TRACE_VTCS
    void             append(const MinIntegral& other)               { vertices.insert(vertices.end(), other.vertices.begin(), other.vertices.end()); }
    void             push_back(const MergeTreeNode::ValueVertex& v) { vertices.push_back(v); }
#else
    void             append(const MinIntegral& other)               { }
    void             push_back(const MergeTreeNode::ValueVertex& v) { }
#endif
    MergeTreeNode::Vertex   min_vtx;
    MergeTreeNode::Value    min_val;
    Real                    integral;
    size_t                  n_cells;
#ifdef REEBER_PERSISTENT_INTEGRAL_TRACE_VTCS
    std::vector<MergeTreeNode::ValueVertex> vertices;
#endif
};

#ifdef REEBER_PERSISTENT_INTEGRAL_TRACE_VTCS
namespace diy {
    template<>
    struct Serialization<MinIntegral>
    {
        static void      save(diy::BinaryBuffer& bb, const MinIntegral &mi)
        {
            diy::save(bb, mi.min_vtx);
            diy::save(bb, mi.min_val);
            diy::save(bb, mi.integral);
            diy::save(bb, mi.n_cells);
            diy::save(bb, mi.vertices);
        }
        static void      load(diy::BinaryBuffer& bb, MinIntegral &mi)
        {
            diy::load(bb, mi.min_vtx);
            diy::load(bb, mi.min_val);
            diy::load(bb, mi.integral);
            diy::load(bb, mi.n_cells);
            diy::load(bb, mi.vertices);
        }
    };
}
#endif

struct PersistentIntegralBlock
{
    typedef std::vector<MinIntegral>                                MinIntegralVector;
    int                                                             gid;
    MergeTreeBlock::Box                                             local;
    MergeTreeBlock::Box                                             global;
    MinIntegralVector                                               persistent_integrals;

                     PersistentIntegralBlock()                      { }
                     PersistentIntegralBlock(const MergeTreeBlock& mtb) :
                         gid(mtb.gid), local(mtb.local),
                         global(mtb.global), persistent_integrals() { }
    void             addIntegral(const MinIntegral& mi)             { persistent_integrals.push_back(mi); }
    static void*     create()                                       { return new PersistentIntegralBlock; }
    static void      destroy(void*b)  { delete static_cast<PersistentIntegralBlock*>(b); }
    static void      save(const void *b, diy::BinaryBuffer& bb)     { diy::save(bb, *static_cast<const PersistentIntegralBlock*>(b)); }
    static void      load(      void *b, diy::BinaryBuffer& bb)     { diy::load(bb, *static_cast<PersistentIntegralBlock*>(b)); }
};

class TreeTracer
{
    private:
        typedef std::map<MergeTreeNode::Vertex, MinIntegral>        MinIntegralMap;

    public:
                     TreeTracer(const Decomposer& decomposer_, diy::Master& pi_master_, Real m_, Real t_, Real cell_volume_) :
                         decomposer(decomposer_), pi_master(pi_master_), m(m_), t(t_), cell_volume(cell_volume_)
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
                    pi_block->addIntegral(it->second);
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

                    int dest_gid = decomposer.point_to_gid(block.local.position(mi.min_vtx));
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
            if (decomposer.point_to_gid(block.local.position(n->vertex)) == block.gid) // FIXME: More efficient way? block.local also contains ghosts, so using that won't work
            {
                mi_res.integral += n->value * cell_volume;
                ++mi_res.n_cells;
                mi_res.push_back(MergeTreeNode::ValueVertex(n->value, n->vertex));
            }

            BOOST_FOREACH(const MergeTree::Node::ValueVertex& x, n->vertices)
                if (decomposer.point_to_gid(block.local.position(x.second)) == block.gid && block.mt.cmp(x.first, t)) // FIXME: More efficient way? block.local also contains ghosts, so using that won't work
                {
                    mi_res.integral += x.first * cell_volume;
                    ++mi_res.n_cells;
                    mi_res.push_back(x);
                }

            // Find contribution from children + min
            BOOST_FOREACH(Neighbor child, n->children)
            {
                MinIntegral mi = integrate(child, block, gp);
                if (block.mt.cmp(mi.min_val, mi_res.min_val))
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
        Real               cell_volume;
};

bool vv_cmp(const MergeTreeNode::ValueVertex& a, const MergeTreeNode::ValueVertex& b)
{
    return a.second < b.second;
}

class OutputIntegrals {
    public:
                    OutputIntegrals(std::string outfn_) : outfn(outfn_) {}
       void         operator()(void *b, const diy::Master::ProxyWithLink& cp, void* aux) const
       {
           PersistentIntegralBlock&  block = *static_cast<PersistentIntegralBlock*>(b);

           std::string   dgm_fn = fmt::format("{}-b{}.comp", outfn, block.gid);
           std::ofstream ofs(dgm_fn.c_str());
 
           MergeTreeBlock::OffsetGrid::GridProxy gp(0, block.global.shape());
           BOOST_FOREACH(MinIntegral &mi, block.persistent_integrals)
           {
               ofs << block.local.position(mi.min_vtx) << " (" << mi.min_vtx << " " << mi.min_val << ") " << mi.integral << " " << mi.n_cells;
#ifdef REEBER_PERSISTENT_INTEGRAL_TRACE_VTCS
               ofs << " [ ";
               std::sort(mi.vertices.begin(), mi.vertices.end(), vv_cmp);
               for (std::vector<MergeTreeNode::ValueVertex>::const_iterator it = mi.vertices.begin(); it != mi.vertices.end(); ++it)
                   ofs << "(" << block.local.position(it->second) << " ," <<it->second << ",  " << it->first << ") ";
               ofs << "]";
#endif
               ofs << std::endl;
           }
       }

    private:
        std::string outfn;
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
    Real        t           = 80;
    Real        cell_volume = 1;;

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
        >> Option('v', "volume",    cell_volume,  "cell volume")
    ;
    //bool verbose = ops >> Present('v', "verbose", "verbose output");

    std::string infn, outfn;
    if (  ops >> Present('h', "help", "show help message") ||
        !(ops >> PosOption(infn) >> PosOption(outfn)))
    {
        fmt::print("Usage: {} IN.lgt OUT.pi\n{}", argv[0], ops);
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
        std::string profile_fn = fmt::format("{}.prf", profile_path);
        profile_stream.open(profile_fn.c_str());
        dlog::prof.add_stream(profile_stream);
    }

    LOG_SEV(info) << "Starting computation";
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
    LOG_SEV(debug) << "Reading blocks from " << infn;
    diy::io::read_blocks(infn, world, assigner, mt_master);
    LOG_SEV(info) << "Blocks read: " << mt_master.size();

    // get the domain bounds from any block that's in memory (they are all the same) and set up a decomposer
    MergeTreeBlock::Box global = static_cast<MergeTreeBlock*>(((const diy::Master&) mt_master).block(mt_master.loaded_block()))->global;
    diy::DiscreteBounds domain;
    for (unsigned i = 0; i < 3; ++i)
    {
        domain.min[i] = global.from()[i];
        domain.max[i] = global.to()[i];
    }
    diy::RegularDecomposer<diy::DiscreteBounds>     decomposer(3, domain, assigner, Decomposer::BoolVector(3, true));

    // Compute and combine persistent integrals
    diy::all_to_all(mt_master, assigner, TreeTracer(decomposer, pi_master, m, t, cell_volume), k);
    // Save persistent integrals to file
    pi_master.foreach(OutputIntegrals(outfn));
    dlog::prof.flush();
}
