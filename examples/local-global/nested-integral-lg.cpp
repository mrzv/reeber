#include <iostream>
#include <string>

#include <boost/algorithm/string/split.hpp>

#include <opts/opts.h>
#include <dlog/log.h>
#include <dlog/stats.h>

#include <diy/decomposition.hpp>
#include <diy/reduce-operations.hpp>
#include <diy/io/block.hpp>

#include <reeber/format.h>

#include "reader-interfaces.h"

#include "merge-tree-block.h"
#include "persistent-integral-block.h"

typedef diy::RegularDecomposer<diy::DiscreteBounds>                 Decomposer;
typedef MergeTreeBlock::Vertex                                      Vertex;
typedef MergeTreeBlock::Value                                       Value;

bool vv_cmp(const MergeTreeNode::ValueVertex& a, const MergeTreeNode::ValueVertex& b)
{
    return a.second < b.second;
}

class TreeTracer
{
    private:
        typedef std::map<MergeTreeNode::Vertex, MinIntegral>        MinIntegralMap;
        typedef MergeTreeBlock::OffsetGrid                          OffsetGrid;
        typedef std::vector<OffsetGrid*>                            OffsetGridVector;

    public:
                    TreeTracer(const Decomposer& decomposer_, diy::Master& pi_master_,
                               Real m_, Real t_, Real e_,
                               std::vector<std::string> avg_fn_list_,
                               std::string density_fn_,
                               bool density_weighted_):
                         decomposer(decomposer_), pi_master(pi_master_),
                         m(m_), t(t_), e(e_),
                         density_reader(0), density_weighted(density_weighted_)
        {
            diy::mpi::communicator& world = pi_master_.communicator();

            for(const std::string& fn : avg_fn_list_)
                avg_var_readers.push_back(Reader::create(fn, world));

            if (!density_fn_.empty())
                density_reader = Reader::create(density_fn_, world);
        }

                    ~TreeTracer()
        {
            for(Reader* reader : avg_var_readers)
                delete reader;
            delete density_reader;
        }

        void        operator()(void* b, const diy::ReduceProxy& rp) const
        {
            MergeTreeBlock &block = *static_cast<MergeTreeBlock*>(b);

            MinIntegralMap mi_map;

            if (rp.in_link().size() == 0)
            {
                OffsetGridVector add_data;
                for(Reader* reader : avg_var_readers)
                    add_data.push_back(reader->read(block.core));
                MergeTreeBlock::OffsetGrid *density_data = density_reader ? density_reader->read(block.core) : 0;

                reeber::traverse_persistence(block.mt, Integrator(m,t,e,density_weighted, &mi_map, &block, &add_data, density_data));

                for(const MinIntegral& mi : mi_map | reeber::range::map_values)
                {
                    if (block.mt.cmp(m, mi.min_val))
                        continue;

                    if (mi.integral == 0)                                            // this block doesn't contribute anything to this extremum
                        continue;

                    int dest_gid = decomposer.point_to_gid(block.global.position(mi.min_vtx));
                    diy::BlockID dest = rp.out_link().target(dest_gid);              // out_link targets are ordered as gids
                    assert(dest.gid == dest_gid);
                    rp.enqueue(dest, mi);
                }

                for(MergeTreeBlock::OffsetGrid* og : add_data)
                    delete og;
                delete density_data;
            }
            else
            {
                for (int i = 0; i < rp.in_link().size(); ++i)
                {
                    int gid = rp.in_link().target(i).gid;
                    assert(gid == i);

                    while (rp.incoming(gid))
                    {
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

        struct Integrator
        {
            Real                m,t,e;
            bool                density_weighted;
            MinIntegralMap*     mi_map_;
            MergeTreeBlock*     b;
            OffsetGridVector*   add_data_;
            OffsetGrid*         density_data;

                        Integrator(Real m_, Real t_, Real e_, bool density_weighted_,
                                   MinIntegralMap* mi_map, MergeTreeBlock* b_, OffsetGridVector* add_data, OffsetGrid* density_data_):
                            m(m_), t(t_), e(e_), density_weighted(density_weighted_),
                            mi_map_(mi_map), b(b_), add_data_(add_data), density_data(density_data_)    {}

            void        operator()(Neighbor from, Neighbor through, Neighbor to) const
            {
                MinIntegralMap& mi_map = *mi_map_;

                // contribution of the edge from -- from->parent to the integral of from
                mi_map[from->vertex].combine(integrate(from));

                // figure out the contribution of the edge through -- through->parent to the integral of to
                bool record_through = through->children.size() == 2;
                if (!record_through && from != to)    // if we have more than 2 children record through when processing the highest one (lexicographically)
                {
                    Neighbor n = 0;
                    for (size_t i = 0; i < through->children.size(); ++i)
                    {
                        Neighbor nc = through->children[i];
                        if (static_cast<Neighbor>(nc->aux) != to && nc > n)
                            n = nc;
                    }

                    if (n == 0)
                        assert(from == to);
                    else if (static_cast<Neighbor>(n->aux) == from)
                        record_through = true;
                }

                if (record_through && from != to)
                    mi_map[to->vertex].combine(integrate(through));

                // if we are not looking at the global pair, add the entire from-integral (now correct) to the integral of to
                if (from != to && !b->mt.cmp(t, through->value))
                    mi_map[to->vertex].combine(mi_map[from->vertex]);

                // from-through pair is not persistent enough, erase
                if (from->value/through->value < e)     // we should make this more generic somehow
                    mi_map.erase(from->vertex);
                else
                {
                    // fix the integral id (initialized to zero up until now)
                    MinIntegral mi(from);
                    mi.combine(mi_map[from->vertex]);
                    mi_map[from->vertex] = mi;
                }
            }

            void        integrate(MinIntegral& mi, Value val, MergeTree::Vertex vrt) const
            {
                OffsetGridVector& add_data = *add_data_;

                if (b->core.contains(vrt) && !b->mt.cmp(t, val))
                {
                    mi.integral += val * b->cell_size[0] * b->cell_size[1] * b->cell_size[2];
                    ++mi.n_cells;
                    for (size_t i = 0; i < add_data.size(); ++i)
                    {
                        Real new_val = (*add_data[i])(vrt);
                        if (density_data)
                            new_val /= (*density_data)(vrt);
                        if (density_weighted)
                            new_val *= val * b->cell_size[0] * b->cell_size[1] * b->cell_size[2];
                        mi.add_sums[i] += new_val;
                    }
                    mi.push_back(MergeTreeNode::ValueVertex(val, vrt));
                }
            }

            MinIntegral integrate(Neighbor n) const
            {
                MinIntegral mi;

                integrate(mi, n->value, n->vertex);

                for(const MergeTree::Node::ValueVertex& x : n->vertices)
                    integrate(mi, x.first, x.second);

                return mi;
            }
        };

    private:
        const Decomposer&    decomposer;
        diy::Master&         pi_master;
        Real                 m;
        Real                 t;
        Real                 e;
        std::vector<Reader*> avg_var_readers;
        Reader*              density_reader;
        bool                 density_weighted;
};

struct OutputIntegrals
{
                    OutputIntegrals(std::string outfn_, bool density_weighted_, bool verbose_):
                        outfn(outfn_), density_weighted(density_weighted_), verbose(verbose_)   {}

       void         operator()(PersistentIntegralBlock* b, const diy::Master::ProxyWithLink& cp) const
       {
           PersistentIntegralBlock&  block = *b;

           std::string   dgm_fn = fmt::format("{}-b{}.comp", outfn, block.gid);
           std::ofstream ofs(dgm_fn.c_str());
 
           for(MinIntegral &mi : block.persistent_integrals)
           {
               Vertex v = block.global.position(mi.min_vtx);
               ofs << v[0] * block.cell_size[0] << " " << v[1] * block.cell_size[1] << " " << v[2] * block.cell_size[2] << " ";
               if (verbose)
                   ofs << v[0] << "x" << v[1] << "x" << v[2] << " (" << mi.min_vtx << ") ";
               ofs <<  mi.integral;
               if (verbose)
                   ofs << " " << mi.n_cells;
               for(Real sum : mi.add_sums)
                   ofs << " " << sum / (density_weighted ? mi.integral : mi.n_cells);
               ofs << std::endl;
#ifdef REEBER_PERSISTENT_INTEGRAL_TRACE_VTCS
               std::sort(mi.vertices.begin(), mi.vertices.end(), vv_cmp);
               for (std::vector< MergeTreeNode::ValueVertex >::const_iterator it = mi.vertices.begin(); it != mi.vertices.end(); ++it)
                   ofs << "   " << it->second << " (" << block.global.position(it->second) <<  ")" << std::endl;
#if 1
               // Consistency check for debugging
               for (size_t i = 0; i < mi.vertices.size() - 1; ++i)
                   if (mi.vertices[i].second == mi.vertices[i+1].second)
                       LOG_SEV(fatal) << "Duplicate vertex " << mi.vertices[i].second << " in component " << mi.min_vtx;
#endif
#endif
           }
       }

    private:
        std::string outfn;
        bool        density_weighted;
        bool        verbose;
};

int main(int argc, char** argv)
{
    diy::mpi::environment   env(argc, argv);
    diy::mpi::communicator  world;
#ifdef REEBER_USE_BOXLIB_READER
    reeber::io::BoxLib::environment boxlib_env(argc, argv, world);
#endif

    using namespace opts;

    std::string prefix      = "./DIY.XXXXXX";
    int         in_memory   = -1;
    int         threads     = 1;
    int         k           = 2;
    Real        m           = 200;
    Real        t           = 82;
    Real        e           = m - t;

    std::string profile_path;
    std::string log_level = "info";
    std::string avg_fn_str = "";
    std::string density_fn = "";

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
        >> Option('e', "epsilon",   e,            "persistence threshold")
        >> Option('f', "mean",      avg_fn_str,   "list of additionals files/variables to average separated by ','")
        >> Option('q', "quotient",  density_fn,   "divide by density in file")
    ;
    bool absolute         = ops >> Present('a', "absolute", "use absolute values for thresholds (instead of multiples of mean)");
    bool verbose          = ops >> Present('v', "verbose",  "verbose output: logical coordiantes and number of cells");
    bool density_weighted = ops >> Present('w', "weight",   "compute density-weighted averages");
    bool split            = ops >> Present(     "split",    "use split IO");

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

    std::vector<std::string> avg_fn_list;
    if (!avg_fn_str.empty())
        boost::split(avg_fn_list, avg_fn_str, std::bind1st(std::equal_to<char>(), ','));

    dlog::add_stream(std::cerr, dlog::severity(log_level))
        << dlog::stamp() << dlog::aux_reporter(world.rank()) << dlog::color_pre() << dlog::level() << dlog::color_post() >> dlog::flush();

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
    if (!split)
        diy::io::read_blocks(infn, world, assigner, mt_master);
    else
        diy::io::split::read_blocks(infn, world, assigner, mt_master);
    LOG_SEV_IF(world.rank() == 0, info) << "Blocks read: " << mt_master.size();

    world.barrier();
    LOG_SEV_IF(world.rank() == 0, info) << "Time to read data:                    " << dlog::clock_to_string(timer.elapsed());
    timer.restart();

    // get the domain bounds from any block that's in memory (they are all the same) and set up a decomposer
    MergeTreeBlock::Box global = static_cast<MergeTreeBlock*>(((const diy::Master&) mt_master).block(mt_master.loaded_block()))->global;
    diy::DiscreteBounds domain {3};
    for (unsigned i = 0; i < 3; ++i)
    {
        domain.min[i] = global.from()[i];
        domain.max[i] = global.to()[i];
    }
    diy::RegularDecomposer<diy::DiscreteBounds>     decomposer(3, domain, assigner.nblocks(), Decomposer::BoolVector(3, true));

    // Compute average
    if (!absolute)
    {
        mt_master.foreach(&MergeTreeBlock::compute_average);
        mt_master.exchange();

        const diy::Master::ProxyWithLink& proxy = mt_master.proxy(mt_master.loaded_block());
        double mean = proxy.get<double>() / proxy.get<size_t>();
        m *= mean;
        t *= mean;

        LOG_SEV_IF(world.rank() == 0, info) << "Average value is " << mean << ". Using isofind threshold of " << t << " and maximum threshold of " << m;
    }

    // Compute and combine persistent integrals
    diy::all_to_all(mt_master, assigner, TreeTracer(decomposer, pi_master, m, t, e, avg_fn_list, density_fn, density_weighted), k);

    world.barrier();
    LOG_SEV_IF(world.rank() == 0, info) << "Time to compute persistent integrals: " << dlog::clock_to_string(timer.elapsed());
    timer.restart();

    // Save persistent integrals to file
    pi_master.foreach(OutputIntegrals(outfn, density_weighted, verbose));

    world.barrier();
    LOG_SEV_IF(world.rank() == 0, info) << "Time to output persistent integrals:  " << dlog::clock_to_string(timer.elapsed());
    timer.restart();

    dlog::prof.flush();     // TODO: this is necessary because the profile file will close before
                            //       the global dlog::prof goes out of scope and flushes the events.
                            //       Need to eventually fix this.
    dlog::stats.flush();
}
