//#define SEND_COMPONENTS

#define DO_DETAILED_TIMING

#include "reeber-real.h"

// to print nice backtrace on segfault signal
#include <signal.h>
#include <execinfo.h>
#include <cxxabi.h>

#include <diy/master.hpp>
#include <diy/io/block.hpp>
#include <diy/io/shared.hpp>
#include <diy/decomposition.hpp>

#include <dlog/stats.h>
#include <dlog/log.h>
#include <opts/opts.h>

#include <AMReX.H>
#include <error.h>

#include "fab-block.h"
#include "fab-tmt-block.h"
#include "reader-interfaces.h"

#include "../../local-global/output-persistence.h"

#include "read-npy.h"

#ifdef AMR_MT_SEND_COMPONENTS
#include "amr-merge-tree-send-componentwise.h"
#else
#include "amr-merge-tree-send-simple.h"
#endif

#include "amr-plot-reader.h"


// block-independent types
using AMRLink = diy::AMRLink;

using Bounds = diy::DiscreteBounds;
using AmrVertexId = r::AmrVertexId;
using AmrEdge = reeber::AmrEdge;

#define DIM 3

using FabBlockR = FabBlock<Real, DIM>;

using Block = FabTmtBlock<Real, DIM>;
using Vertex = Block::Vertex;
using Component = Block::Component;
using VertexNeighborMap = Block::TripletMergeTree::VertexNeighborMap;
using AmrTripletMergeTree = Block::TripletMergeTree;
using MaskedBox = Block::MaskedBox;
using GidVector = Block::GidVector;
using GidSet = Block::GidContainer;

using Neighbor = AmrTripletMergeTree::Neighbor;

using Diagram = Block::Diagram;

constexpr bool abort_on_segfault = true;

struct IsAmrVertexLocal
{
    bool operator()(const Block& b, const Neighbor& from) const
    {
        return from->vertex.gid == b.gid;
    }
};

template<class Real, class LocalFunctor>
struct ComponentDiagramsFunctor
{

    ComponentDiagramsFunctor(Block* b, const LocalFunctor& lf)
            :
            block_(b),
            negate_(b->get_merge_tree().negate()),
            ignore_zero_persistence_(true),
            test_local(lf)
    { }

    void operator()(Neighbor from, Neighbor through, Neighbor to) const
    {
        if (!test_local(*block_, from))
            return;

        AmrVertexId current_vertex = from->vertex;

        Real birth_time = from->value;
        Real death_time = through->value;

        if (ignore_zero_persistence_ and birth_time == death_time)
            return;

        AmrVertexId root = block_->final_vertex_to_deepest_[current_vertex];
        block_->local_diagrams_[root].emplace_back(birth_time, death_time);
    }

    Block* block_;
    const bool negate_;
    const bool ignore_zero_persistence_;
    LocalFunctor test_local;
};

using OutputPairsR = OutputPairs<Block, IsAmrVertexLocal>;

inline bool file_exists(const std::string& s)
{
    std::ifstream ifs(s);
    return ifs.good();

}

inline bool ends_with(const std::string& s, const std::string& suffix)
{
    if (suffix.size() > s.size())
        return false;
    return std::equal(suffix.rbegin(), suffix.rend(), s.rbegin());
}

void read_from_file(std::string infn,
        diy::mpi::communicator& world,
        diy::Master& master_reader,
        diy::ContiguousAssigner& assigner,
        diy::MemoryBuffer& header,
        diy::DiscreteBounds& domain,
        bool split,
        int nblocks)
{
    if (not file_exists(infn))
        throw std::runtime_error("Cannot read file " + infn);

    if (ends_with(infn, ".npy"))
    {
        read_from_npy_file<DIM>(infn, world, nblocks, master_reader, assigner, header, domain);
    } else
    {
        if (split)
            diy::io::split::read_blocks(infn, world, assigner, master_reader, header, FabBlockR::load);
        else
            diy::io::read_blocks(infn, world, assigner, master_reader, header, FabBlockR::load);
        diy::load(header, domain);
    }
}

void catch_sig(int signum)
{
    LOG_SEV_IF(true, fatal) << "caught signal " << signum;
    dlog::flush();
    //    << ",  local group = " << pro {}, local rank = {}", signum, active_puppet, proc_map->group(), proc_map->local_rank());

    // print backtrace
    void* callstack[128];
    int frames = backtrace(callstack, 128);
    char** strs = backtrace_symbols(callstack, frames);

    size_t funcnamesize = 256;
    char* funcname = (char*) malloc(funcnamesize);

    // iterate over the returned symbol lines. skip the first, it is the
    // address of this function.
    for(int i = 1; i < frames; i++)
    {
        char* begin_name = 0, * begin_offset = 0, * end_offset = 0;

        // find parentheses and +address offset surrounding the mangled name:
        // ./module(function+0x15c) [0x8048a6d]
        for(char* p = strs[i]; *p; ++p)
        {
            if (*p == '(')
                begin_name = p;
            else if (*p == '+')
                begin_offset = p;
            else if (*p == ')' && begin_offset)
            {
                end_offset = p;
                break;
            }
        }

        if (begin_name && begin_offset && end_offset && begin_name < begin_offset)
        {
            *begin_name++ = '\0';
            *begin_offset++ = '\0';
            *end_offset = '\0';

            // mangled name is now in [begin_name, begin_offset) and caller
            // offset in [begin_offset, end_offset). now apply __cxa_demangle():

            int status;
            char* ret = abi::__cxa_demangle(begin_name, funcname, &funcnamesize, &status);
            if (status == 0)
            {
                funcname = ret; // use possibly realloc()-ed string
                LOG_SEV_IF(true, fatal) << "  " << strs[i] << " : " << funcname << "+" << begin_offset;
                dlog::flush();
            } else
            {
                // demangling failed. Output function name as a C function with no arguments.
                //                logger->critical("  {} : {}()+{}", strs[i], begin_name, begin_offset);
                LOG_SEV_IF(true, fatal)  << "  " << strs[i] << " : " << funcname << "+" << begin_offset;
                dlog::flush();
            }
        } else
        {
            // couldn't parse the line? print the whole line.
            LOG_SEV_IF(true, fatal)  << "  " << strs[i];
            dlog::flush();
            //            logger->critical("  {}", strs[i]);
        }
    }

    free(funcname);
    free(strs);

    signal(signum, SIG_DFL);    // restore the default signal
    if (abort_on_segfault)
        MPI_Abort(MPI_COMM_WORLD, 1);
}

int main(int argc, char** argv)
{
    bool debug = false;
    signal(SIGSEGV, catch_sig);
    signal(SIGABRT, catch_sig);

    diy::mpi::environment env(argc, argv);
    diy::mpi::communicator world;

//    if (argc < 2)
//    {
//        fmt::print(std::cerr, "Usage: {} IN.amr rho\n", argv[0]);
//        return 1;
//    }

    int nblocks = world.size();
    std::string prefix = "./DIY.XXXXXX";
    int in_memory = -1;
    int threads = 1;
    std::string profile_path;
    std::string log_level = "info";

    // threshold
    Real rho = 81.66;
    Real absolute_rho;
    int min_cells = 10;
    int n_runs = 1;

    std::string fields_to_read;

    using namespace opts;

    //-----------
    opts::Options ops(argc, argv);
    ops
            >> Option('b', "blocks", nblocks, "number of blocks to use")
            >> Option('m', "memory", in_memory, "maximum blocks to store in memory")
            >> Option('j', "jobs", threads, "threads to use during the computation")
            >> Option('s', "storage", prefix, "storage prefix")
            >> Option('i', "rho", rho, "iso threshold")
            >> Option('x', "mincells", min_cells, "minimal number of cells to output halo")
            >> Option('f', "fields", fields_to_read, "comma-separated list of fields to read")
            >> Option('r', "runs", n_runs, "number of runs")
            >> Option('p', "profile", profile_path, "path to keep the execution profile")
            >> Option('l', "log", log_level, "log level");

    bool absolute =
            ops >> Present('a', "absolute", "use absolute values for thresholds (instead of multiples of mean)");
    bool read_plotfile = ops >> Present("plotfile", "read AMR plotfiles");
    bool negate = ops >> opts::Present('n', "negate", "sweep superlevel sets");
    // ignored for now, wrap is always assumed
    bool wrap = ops >> opts::Present('w', "wrap", "wrap");
    bool split = ops >> opts::Present("split", "use split IO");

    bool print_stats = ops >> opts::Present("stats", "print statistics");
    std::string input_filename, output_filename, output_diagrams_filename, output_integral_filename;

    std::vector<std::string> all_var_names = split_by_delim(fields_to_read,
            ',');  //{"particle_mass_density", "density", "xmom", "ymom", "zmom"};

    int n_mt_vars = all_var_names.size();

    if (all_var_names.empty() or
            ops >> Present('h', "help", "show help message") or
            not(ops >> PosOption(input_filename)) or
            not(ops >> PosOption(output_filename)))
    {
        if (world.rank() == 0)
        {
            fmt::print("Usage: {} INPUT-PLTFILE OUTPUT.mt [OUT_DIAGRAMS] [OUT_INTEGRAL]\n", argv[0]);
            fmt::print("Compute local-global tree from AMR data\n");
            fmt::print("{}", ops);
        }
        return 1;
    }

    //-----------

    bool write_diag = (ops >> PosOption(output_diagrams_filename));
    if (output_diagrams_filename == "none")
        write_diag = false;

    bool write_integral = (ops >> PosOption(output_integral_filename));
    if (output_integral_filename == "none")
        write_integral = false;

    diy::FileStorage storage(prefix);

    diy::Master master_reader(world, 1, in_memory, &FabBlockR::create, &FabBlockR::destroy);
    diy::ContiguousAssigner assigner(world.size(), nblocks);
    diy::MemoryBuffer header;
    diy::DiscreteBounds domain(DIM);

    dlog::add_stream(std::cerr, dlog::severity(log_level))
            << dlog::stamp() << dlog::aux_reporter(world.rank()) << dlog::color_pre() << dlog::level()
            << dlog::color_post() >> dlog::flush();

    world.barrier();
    dlog::Timer timer;
    dlog::Timer timer_all;
    LOG_SEV_IF(world.rank() == 0, info) << "Starting computation, input_filename = " << input_filename << ", nblocks = "
                                                                                     << nblocks
                                                                                     << ", rho = " << rho
                                                                                     << ", summing first " << n_mt_vars
                                                                                     << " of " << container_to_string(
            all_var_names);
    dlog::flush();
    world.barrier();

#ifdef DO_DETAILED_TIMING
    // detailed timings
    using DurationType = decltype(timer.elapsed());

    dlog::Timer timer_send;
    dlog::Timer timer_receieve;
    dlog::Timer timer_tmt_exchange;

    DurationType time_to_construct_blocks;
    DurationType time_to_init_blocks;
    DurationType time_to_get_average;
    DurationType tmt_send_time = timer_send.elapsed();
    DurationType tmt_receive_time = timer_receieve.elapsed();
    DurationType tmt_exchange_1_time = timer_receieve.elapsed();
    DurationType tmt_exchange_2_time = timer_receieve.elapsed();
    DurationType time_to_delete_low_edges;

    DurationType min_time_to_construct_blocks;
    DurationType min_time_to_init_blocks;
    DurationType min_time_to_get_average;
    DurationType min_tmt_send_time;
    DurationType min_tmt_receive_time;
    DurationType min_tmt_exchange_1_time;
    DurationType min_tmt_exchange_2_time;
    DurationType min_time_to_delete_low_edges;

    DurationType max_time_to_construct_blocks;
    DurationType max_time_to_init_blocks;
    DurationType max_time_to_get_average;
    DurationType max_tmt_send_time;
    DurationType max_tmt_receive_time;
    DurationType max_tmt_exchange_1_time;
    DurationType max_tmt_exchange_2_time;
    DurationType max_time_to_delete_low_edges;

#endif

    if (read_plotfile)
    {
        read_amr_plotfile(input_filename, all_var_names, n_mt_vars, world, nblocks, master_reader, header, domain);
    } else
    {
        read_from_file(input_filename, world, master_reader, assigner, header, domain, split, nblocks);
    }

    world.barrier();

    auto time_to_read_data = timer.elapsed();
    LOG_SEV_IF(world.rank() == 0, info) << "Data read, local size = " << master_reader.size();
    LOG_SEV_IF(world.rank() == 0, info) << "Time to read data:       " << dlog::clock_to_string(timer.elapsed());
    dlog::flush();
    for(int n_run = 0; n_run < n_runs; ++n_run)
    {

        diy::Master master(world, threads, in_memory, &Block::create, &Block::destroy, &storage, &Block::save,
                &Block::load);
        timer.restart();
        timer_all.restart();

        // copy FabBlocks to FabTmtBlocks
        // in FabTmtConstructor mask will be set and local trees will be computed
        // FabBlock can be safely discarded afterwards

        master_reader.foreach(
                [&master, domain, rho, negate, absolute](FabBlockR* b, const diy::Master::ProxyWithLink& cp) {
                    auto* l = static_cast<AMRLink*>(cp.link());
                    AMRLink* new_link = new AMRLink(*l);

                    // prepare neighbor box info to save in MaskedBox
                    // TODO: vector refinement
                    int local_ref = l->refinement()[0];
                    int local_lev = l->level();

                    master.add(cp.gid(),
                            new Block(b->fab, local_ref, local_lev, domain, l->bounds(), l->core(), cp.gid(),
                                    new_link, rho, negate, absolute),
                            new_link);

                });

        auto time_for_local_computation = timer.elapsed();

#ifdef DO_DETAILED_TIMING
        time_to_construct_blocks = timer.elapsed();
#endif

        Real mean = std::numeric_limits<Real>::min();

        if (debug)
        {
            master.foreach([debug](Block* b, const diy::Master::ProxyWithLink& cp) {
                auto* l = static_cast<diy::AMRLink*>(cp.link());
                {
                    fmt::print(
                            "master, FabTmtBlocks: gid = {}: level = {}, shape = {}, core = {} - {}, bounds = {} - {}, sum_  = {}, n_active = {}, local = {}\n",
                            cp.gid(), l->level(), b->fab_.shape(),
                            l->core().min, l->core().max,
                            l->bounds().min, l->bounds().max,
                            b->sum_,
                            b->n_active_,
                            b->local_);
                }
            });
        }

        if (absolute)
        {
            absolute_rho = rho;
            LOG_SEV_IF(world.rank() == 0, info) << "Time to compute local trees and components:  "
                    << dlog::clock_to_string(timer.elapsed());
        } else
        {
            LOG_SEV_IF(world.rank() == 0, info) << "Time to construct FabTmtBlocks: "
                    << dlog::clock_to_string(timer.elapsed());
            dlog::flush();
            timer.restart();

            master.foreach([debug](Block* b, const diy::Master::ProxyWithLink& cp) {
                cp.collectives()->clear();
                cp.all_reduce(b->sum_ * b->scaling_factor(), std::plus<Real>());
                cp.all_reduce(static_cast<Real>(b->n_unmasked_) * b->scaling_factor(), std::plus<Real>());
                if (debug)
                    fmt::print("BEFORE EXCHANgE gid = {}, sum = {}, n_unmasked = {}\n", b->gid, b->sum_,
                            b->n_unmasked_);
            });

            master.exchange();

            const diy::Master::ProxyWithLink& proxy = master.proxy(master.loaded_block());

            mean = proxy.get<Real>() / proxy.get<Real>();
            absolute_rho = rho * mean;                                            // now rho contains absolute threshold

            LOG_SEV_IF(world.rank() == 0, info) << "Average = " << mean << ", rho = " << rho
                                                                << ", absolute_rho =  " << absolute_rho
                                                                << ", time to compute average: "
                                                                << dlog::clock_to_string(timer.elapsed());
#ifdef DO_DETAILED_TIMING
            time_to_get_average = timer.elapsed();
#endif
            time_for_local_computation += timer.elapsed();

            dlog::flush();

            if (mean < 0 or std::isnan(mean) or std::isinf(mean) or mean > 1e+40)
            {
                LOG_SEV_IF(world.rank() == 0, error) << "Bad average = " << mean << ", do not proceed";
                if (read_plotfile)
                    amrex::Finalize();
                return 1;
            }

            timer.restart();
            master.foreach([absolute_rho](Block* b, const diy::Master::ProxyWithLink& cp) {
                AMRLink* l = static_cast<AMRLink*>(cp.link());
                b->init(absolute_rho, l);
                cp.collectives()->clear();
            });

#ifdef DO_DETAILED_TIMING
            time_to_init_blocks = timer.elapsed();
#endif

            LOG_SEV_IF(world.rank() == 0, info) <<
            "Time to initialize FabTmtBlocks (low vertices, local trees, components, outgoing edges): "
                    << timer.elapsed();
            time_for_local_computation += timer.elapsed();
        }
        dlog::flush();
        timer.restart();

        int global_n_undone = 1;

        master.foreach(&send_edges_to_neighbors<DIM>);
        master.exchange();
        master.foreach(&delete_low_edges<DIM>);

        LOG_SEV_IF(world.rank() == 0, info)  << "edges symmetrized, time elapsed " << timer.elapsed();
        auto time_for_communication = timer.elapsed();

#ifdef DO_DETAILED_TIMING
        time_to_delete_low_edges = timer.elapsed();
#endif

        dlog::flush();
        timer.restart();

        if (debug)
        {

            master.foreach([](Block* b, const diy::Master::ProxyWithLink& cp) {
                auto* l = static_cast<AMRLink*>(cp.link());

                std::set<diy::BlockID> receivers;
                for(int i = 0; i < l->size(); ++i)
                {
                    if (l->target(i).gid != b->gid)
                    {
                        receivers.insert(l->target(i));
                    }
                }

                for(const diy::BlockID& receiver : receivers)
                {
                    int receiver_gid = receiver.gid;
                    cp.enqueue(receiver, (int) (b->new_receivers_.count(receiver_gid)));
                }
            });

            master.exchange();

            master.foreach([](Block* b, const diy::Master::ProxyWithLink& cp) {
                auto* l = static_cast<AMRLink*>(cp.link());
                std::set<diy::BlockID> senders;
                for(int i = 0; i < l->size(); ++i)
                {
                    if (l->target(i).gid != b->gid)
                    {
                        senders.insert(l->target(i));
                    }
                }

                for(const diy::BlockID& sender : senders)
                {
                    int t;
                    cp.dequeue(sender, t);

                    if (t != (int) (b->new_receivers_.count(sender.gid)))
                    {
                        throw std::runtime_error("Asymmetry");
                    }

                }
            });

            LOG_SEV_IF(world.rank() == 0, info)  << "Symmetry checked in "
                    << dlog::clock_to_string(timer.elapsed());
            time_for_communication += timer.elapsed();
            dlog::flush();
            timer.restart();
        }

#ifdef DO_DETAILED_TIMING
        tmt_exchange_2_time = 0;
        tmt_exchange_1_time = 0;
        tmt_receive_time = 0;
        tmt_send_time = 0;
#endif

        int rounds = 0;
        while(global_n_undone)
        {
            rounds++;
#ifdef DO_DETAILED_TIMING
            timer_send.restart();
#endif
            master.foreach(&amr_tmt_send<Real, DIM>);

#ifdef DO_DETAILED_TIMING
            tmt_send_time += timer_send.elapsed();
            timer_tmt_exchange.restart();
#endif

            master.exchange();

#ifdef DO_DETAILED_TIMING
            tmt_exchange_1_time += timer_tmt_exchange.elapsed();
            timer_receieve.restart();
#endif
            master.foreach(&amr_tmt_receive<Real, DIM>);

#ifdef DO_DETAILED_TIMING
            tmt_receive_time += timer_receieve.elapsed();
#endif

            LOG_SEV_IF(world.rank() == 0, info) << "MASTER round " << rounds << ", get OK";
            dlog::flush();
#ifdef DO_DETAILED_TIMING
            timer_tmt_exchange.restart();
#endif
            master.exchange();

#ifdef DO_DETAILED_TIMING
            tmt_exchange_2_time += timer_tmt_exchange.elapsed();
#endif
            // to compute total number of undone blocks
            global_n_undone = master.proxy(master.loaded_block()).read<int>();
            LOG_SEV_IF(world.rank() == 0, info) << "MASTER round " << rounds << ", global_n_undone = "
                                                                   << global_n_undone;

            if (print_stats)
            {
                int local_n_undone = 0;
                master.foreach(
                        [&local_n_undone](Block* b, const diy::Master::ProxyWithLink& cp) {
                            local_n_undone += (b->done_ != 1);
                        });

                LOG_SEV(info) << "STAT MASTER round " << rounds << ", rank = " << world.rank() << ", local_n_undone = "
                                                      << local_n_undone;
            }

            dlog::flush();
        }
        world.barrier();

        auto time_for_run = timer_all.elapsed();
        //    fmt::print("world.rank = {}, time for exchange = {}\n", world.rank(), dlog::clock_to_string(timer.elapsed()));

        LOG_SEV_IF(world.rank() == 0, info) << "Time for exchange:  " << dlog::clock_to_string(timer.elapsed());


        LOG_SEV_IF(world.rank() == 0, info) << "Time for run:  " << time_for_run;
        time_for_communication += timer.elapsed();
        dlog::flush();
        timer.restart();

        size_t local_n_active = 0;
        size_t local_n_components = 0;
        size_t local_n_blocks = 0;

        {
            master.foreach([&local_n_active, &local_n_blocks, &local_n_components](Block* b,
                    const diy::Master::ProxyWithLink& cp) {
//            fmt::print("PRINT ADF gid = {}, mt.size = {}\n", b->gid, b->current_merge_tree_.size());
                auto sum_n_vertices_pair = b->get_local_stats();
                int n_low = b->n_low_;
                int n_active = b->n_active_;
                cp.collectives()->clear();
                cp.all_reduce(sum_n_vertices_pair.first, std::plus<Real>());
                cp.all_reduce(sum_n_vertices_pair.second, std::plus<size_t>());
                cp.all_reduce(n_low, std::plus<size_t>());
                cp.all_reduce(n_active, std::plus<size_t>());

                local_n_active += n_active;
                local_n_blocks += 1;
                local_n_components += b->original_deepest_.size();

            });

            master.exchange();

            const diy::Master::ProxyWithLink& proxy = master.proxy(master.loaded_block());

            Real total_value = proxy.get<Real>();
            size_t total_vertices = proxy.get<size_t>();
            size_t total_n_low = proxy.get<size_t>();
            size_t total_n_active = proxy.get<size_t>();

            LOG_SEV_IF(world.rank() == 0, info) << "Total value = " << total_value << ", total # vertices = "
                                                                    << total_vertices
                                                                    << ", mean = " << mean << ", total_n_low = "
                                                                    << total_n_low << ", total_n_active = "
                                                                    << total_n_active << ", master.size = "
                                                                    << master.size();
            dlog::flush();
        }

        // save the result
        if (output_filename != "none")
        {
            if (!split)
                diy::io::write_blocks(output_filename, world, master);
            else
                diy::io::split::write_blocks(output_filename, world, master);
        }

        LOG_SEV_IF(world.rank() == 0, info) << "Time to write tree:  " << dlog::clock_to_string(timer.elapsed());
        auto time_for_output = timer.elapsed();
        dlog::flush();
        timer.restart();

        bool verbose = false;

        if (write_diag)
        {
            bool ignore_zero_persistence = true;
            OutputPairsR::ExtraInfo extra(output_diagrams_filename, verbose, world);
            IsAmrVertexLocal test_local;
            master.foreach(
                    [&extra, &test_local, ignore_zero_persistence, absolute_rho](Block* b,
                            const diy::Master::ProxyWithLink& cp) {
                        output_persistence(b, cp, extra, test_local, absolute_rho, ignore_zero_persistence);
                    });
        }

        LOG_SEV_IF(world.rank() == 0, info) << "Time to write diagrams:  " << dlog::clock_to_string(timer.elapsed());
        time_for_output += timer.elapsed();
        dlog::flush();
        timer.restart();

        if (write_integral)
        {
            master.foreach([](Block* b, const diy::Master::ProxyWithLink& cp) {

                bool debug = false;

                if (debug) fmt::print("PI gid = {}, started final_cc\n", b->gid);
                b->compute_final_connected_components();
                if (debug) fmt::print("PI gid = {}, final_cc done\n", b->gid);

                b->compute_local_integral();
                if (debug) fmt::print("PI gid = {}, compute_local_integral done\n", b->gid);

                AMRLink* l = static_cast<AMRLink*>(cp.link());

                for(const auto& vertex_value_pair : b->local_integral_)
                {
                    int receiver_gid = vertex_value_pair.first.gid;
                    if (receiver_gid == b->gid)
                        continue;
                    auto receiver = l->target(l->find(receiver_gid));
                    cp.enqueue(receiver, vertex_value_pair);
                }
                if (debug) fmt::print("PI gid = {}, enqueing done\n", b->gid);
            });

            master.exchange();

            diy::io::SharedOutFile integral_file(output_integral_filename, world);

            master.foreach([&integral_file](Block* b, const diy::Master::ProxyWithLink& cp) {
                AMRLink* l = static_cast<AMRLink*>(cp.link());
                bool debug = false;
                for(auto bid : l->neighbors())
                {
                    while(cp.incoming(bid.gid))
                    {
                        Block::LocalIntegral::value_type x;
                        cp.dequeue(bid, x);
                        assert(x.first.gid == b->gid);
                        b->local_integral_[x.first] += x.second;
                    }
                }

                if (debug) fmt::print("PI gid = {}, dequeing done\n", b->gid);

                for(const auto& root_value_pair : b->local_integral_)
                {
                    AmrVertexId root = root_value_pair.first;
                    if (root.gid != b->gid)
                        continue;

//                integral_file << fmt::format("{} {} {}\n", root, b->local_.global_position(root), root_value_pair.second);
                    integral_file << fmt::format("{} {}\n", b->local_.global_position(root), root_value_pair.second);
                }
                if (debug) fmt::print("PI gid = {}, writing to file done\n", b->gid);
            });

            world.barrier();
            LOG_SEV_IF(world.rank() == 0, info) << "Time to compute and write integral:  "
                    << dlog::clock_to_string(timer.elapsed());
            time_for_output += timer.elapsed();
            dlog::flush();
            timer.restart();
        }

        master.foreach([](Block* b, const diy::Master::ProxyWithLink& cp) {
            auto sum_n_vertices_pair = b->get_local_stats();
            int n_low = b->n_low_;
            int n_active = b->n_active_;
            cp.collectives()->clear();
            cp.all_reduce(sum_n_vertices_pair.first, std::plus<Real>());
            cp.all_reduce(sum_n_vertices_pair.second, std::plus<size_t>());
            cp.all_reduce(n_low, std::plus<size_t>());
            cp.all_reduce(n_active, std::plus<size_t>());

        });

        master.exchange();

        const diy::Master::ProxyWithLink& proxy = master.proxy(master.loaded_block());

        Real total_value = proxy.get<Real>();
        size_t total_vertices = proxy.get<size_t>();
        size_t total_n_low = proxy.get<size_t>();
        size_t total_n_active = proxy.get<size_t>();

        LOG_SEV_IF(world.rank() == 0, info) << "Total value = " << total_value << ", total # vertices = "
                                                                << total_vertices
                                                                << ", mean = " << mean << ", total_n_low = "
                                                                << total_n_low << ", total_n_active = "
                                                                << total_n_active;
        dlog::flush();

        std::string final_timings = fmt::format("run: {} read: {} local: {} exchange: {} output: {} total: {}\n",
                n_run, time_to_read_data, time_for_local_computation, time_for_communication, time_for_output, time_for_run);
        LOG_SEV_IF(world.rank() == 0, info) << final_timings;
        dlog::flush();
#ifdef DO_DETAILED_TIMING
        std::string final_detailed_timings = fmt::format(
                "run: {} construct_blocks = {} init_blocks = {} average = {} delete_low_edges = {} tmt_send = {} tmt_receive = {} tmt_exchange_1 = {} tmt_exchange_2 = {}\n",
                n_run, time_to_construct_blocks, time_to_init_blocks, time_to_get_average, time_to_delete_low_edges,
                tmt_send_time, tmt_receive_time, tmt_exchange_1_time, tmt_exchange_2_time);
        LOG_SEV_IF(world.rank() == 0, info) << final_detailed_timings;
#endif

#if 0
        // output diagrams componentwise
        if (write_halo_diagrams)
        {
            timer.restart();
            master.foreach([](Block* b, const diy::Master::ProxyWithLink& cp) {
                cp.collectives()->clear();
                IsAmrVertexLocal test_local;
                reeber::traverse_persistence(b->get_merge_tree(),
                                             ComponentDiagramsFunctor<Real, IsAmrVertexLocal>(b, test_local));
                auto* l = static_cast<AMRLink*>(cp.link());
                for(auto& root_diagram_pair : b->local_diagrams_)
                {
                    auto root_gid = root_diagram_pair.first.gid;
                    if (root_gid == b->gid)
                        continue;
                    diy::BlockID root_block = l->target(l->find(root_gid));
                    cp.enqueue(root_block, root_diagram_pair.first);
                    cp.enqueue(root_block, root_diagram_pair.second);

                }
            });

            master.exchange();

            master.foreach([](Block* b, const diy::Master::ProxyWithLink& cp) {
                auto* l = static_cast<AMRLink*>(cp.link());
                for(const diy::BlockID& sender : link_unique(l, b->gid))
                {
                    while (cp.incoming(sender.gid))
                    {
                        AmrVertexId component_root;
                        Diagram received_diagram;
                        cp.dequeue(sender, component_root);
                        cp.dequeue(sender, received_diagram);
                        b->local_diagrams_[component_root].insert(b->local_diagrams_[component_root].end(),
                                                                  received_diagram.begin(), received_diagram.end());
                    }
                }
            });

            std::string halos_dgm_fname = fmt::format("{0}-halo-diagrams.txt", output_diagrams_filename);

            diy::io::SharedOutFile halo_diagrams_file(halos_dgm_fname, world);

            master.foreach([&halo_diagrams_file](Block* b, const diy::Master::ProxyWithLink& cp) {
                for(const auto& root_diagram_pair : b->local_diagrams_)
                {
                    AmrVertexId root = root_diagram_pair.first;
                    if (root.gid != b->gid)
                        continue;
                    if (b->local_integral_.count(root) == 0)
                        continue;
                    Real integral_value = b->local_integral_[root];
                    for(const auto dgm_point : root_diagram_pair.second)
                    {
                        std::string s = fmt::format("{}\t{}\t{}\t{}\t{}", integral_value, b->gid, b->local_.global_position(root), dgm_point.first, dgm_point.second);
                        halo_diagrams_file << s << "\n";
                    }
                }
                halo_diagrams_file.flush();
            });

            LOG_SEV_IF(world.rank() == 0, info) << "Time to write individual halo diagrams: "
                    << dlog::clock_to_string(timer.elapsed());
            dlog::flush();
            timer.restart();
        }
#endif
//        // cleanup
//        for(int i = master.size() - 1; i >= 0; --i)
//        {
//            master.destroy(i);
//        }

#ifdef DO_DETAILED_TIMING
        if (n_run == n_runs - 1 or n_run == 0)
        {

            master.foreach(
                    [time_to_construct_blocks, time_to_init_blocks, time_to_get_average, tmt_send_time, tmt_exchange_1_time, tmt_receive_time, tmt_exchange_2_time](
                            Block* b, const diy::Master::ProxyWithLink& cp) {
                        cp.collectives()->clear();

                        cp.all_reduce(time_to_construct_blocks, diy::mpi::maximum<DurationType>());
                        cp.all_reduce(time_to_init_blocks, diy::mpi::maximum<DurationType>());
                        cp.all_reduce(time_to_get_average, diy::mpi::maximum<DurationType>());
                        cp.all_reduce(tmt_send_time, diy::mpi::maximum<DurationType>());
                        cp.all_reduce(tmt_exchange_1_time, diy::mpi::maximum<DurationType>());
                        cp.all_reduce(tmt_receive_time, diy::mpi::maximum<DurationType>());
                        cp.all_reduce(tmt_exchange_2_time, diy::mpi::maximum<DurationType>());

                        cp.all_reduce(time_to_construct_blocks, diy::mpi::minimum<DurationType>());
                        cp.all_reduce(time_to_init_blocks, diy::mpi::minimum<DurationType>());
                        cp.all_reduce(time_to_get_average, diy::mpi::minimum<DurationType>());
                        cp.all_reduce(tmt_send_time, diy::mpi::minimum<DurationType>());
                        cp.all_reduce(tmt_exchange_1_time, diy::mpi::minimum<DurationType>());
                        cp.all_reduce(tmt_receive_time, diy::mpi::minimum<DurationType>());
                        cp.all_reduce(tmt_exchange_2_time, diy::mpi::minimum<DurationType>());

                    });

            master.exchange();

            const diy::Master::ProxyWithLink& proxy = master.proxy(master.loaded_block());

            max_time_to_construct_blocks = proxy.get<DurationType>();
            max_time_to_init_blocks = proxy.get<DurationType>();
            max_time_to_get_average = proxy.get<DurationType>();
            max_tmt_send_time = proxy.get<DurationType>();
            max_tmt_exchange_1_time = proxy.get<DurationType>();
            max_tmt_receive_time = proxy.get<DurationType>();
            max_tmt_exchange_2_time = proxy.get<DurationType>();

            min_time_to_construct_blocks = proxy.get<DurationType>();
            min_time_to_init_blocks = proxy.get<DurationType>();
            min_time_to_get_average = proxy.get<DurationType>();
            min_tmt_send_time = proxy.get<DurationType>();
            min_tmt_exchange_1_time = proxy.get<DurationType>();
            min_tmt_receive_time = proxy.get<DurationType>();
            min_tmt_exchange_2_time = proxy.get<DurationType>();

            LOG_SEV_IF(world.rank() == 0, info) << "max_time_to_construct_blocks = " << max_time_to_construct_blocks;
            LOG_SEV_IF(world.rank() == 0, info) << "min_time_to_construct_blocks = " << min_time_to_construct_blocks;
            LOG_SEV_IF(world.rank() == 0, info) << "---------------------------------------";

            LOG_SEV_IF(world.rank() == 0, info) << "max_time_to_init_blocks = " << max_time_to_init_blocks;
            LOG_SEV_IF(world.rank() == 0, info) << "min_time_to_init_blocks = " << min_time_to_init_blocks;
            LOG_SEV_IF(world.rank() == 0, info) << "---------------------------------------";

            LOG_SEV_IF(world.rank() == 0, info) << "max_time_to_get_average = " << max_time_to_get_average;
            LOG_SEV_IF(world.rank() == 0, info) << "min_time_to_get_average = " << min_time_to_get_average;
            LOG_SEV_IF(world.rank() == 0, info) << "---------------------------------------";

            LOG_SEV_IF(world.rank() == 0, info) << "max_tmt_send_time = " << max_tmt_send_time;
            LOG_SEV_IF(world.rank() == 0, info) << "min_tmt_send_time = " << min_tmt_send_time;
            LOG_SEV_IF(world.rank() == 0, info) << "---------------------------------------";

            LOG_SEV_IF(world.rank() == 0, info) << "max_tmt_exchange_1_time = " << max_tmt_exchange_1_time;
            LOG_SEV_IF(world.rank() == 0, info) << "min_tmt_exchange_1_time = " << min_tmt_exchange_1_time;
            LOG_SEV_IF(world.rank() == 0, info) << "---------------------------------------";

            LOG_SEV_IF(world.rank() == 0, info) << "max_tmt_receive_time = " << max_tmt_receive_time;
            LOG_SEV_IF(world.rank() == 0, info) << "min_tmt_receive_time = " << min_tmt_receive_time;
            LOG_SEV_IF(world.rank() == 0, info) << "---------------------------------------";

            LOG_SEV_IF(world.rank() == 0, info) << "max_tmt_exchange_2_time = " << max_tmt_exchange_2_time;
            LOG_SEV_IF(world.rank() == 0, info) << "min_tmt_exchange_2_time = " << min_tmt_exchange_2_time;
            LOG_SEV_IF(world.rank() == 0, info) << "---------------------------------------";

            if (tmt_exchange_2_time == max_tmt_exchange_2_time or
                    tmt_exchange_2_time == min_tmt_exchange_1_time or
                    tmt_receive_time == max_tmt_receive_time or
                    tmt_receive_time == min_tmt_receive_time)
            {
                LOG_SEV(info) << "Rank = " << world.rank() << ", tmt_exchange_2_time = " << tmt_exchange_2_time
                                           << ", receive_time = "
                                           << tmt_receive_time << ", local_blocks = " << local_n_blocks
                                           << ", local_n_active = " << local_n_active
                                           << ", local_n_components = " << local_n_components;
            }

            if (tmt_exchange_2_time == min_tmt_exchange_1_time or
                    tmt_receive_time == max_tmt_receive_time)
            {
                master.foreach([](Block* b, const diy::Master::ProxyWithLink& cp) {

                    if (b->receive_trees_and_gids_time > 0)
                        LOG_SEV(info) << "MAX RECEIVE TIME details, gid = " << b->gid
                                                                            << ", time_to_receive_trees_and_gids = "
                                                                            << b->receive_trees_and_gids_time;
                    if (b->rl_loop_time > 0)
                        LOG_SEV(info) << "MAX RECEIVE TIME details, gid = " << b->gid << ", rl_loop_time = "
                                                                            << b->rl_loop_time;
                    if (b->repair_time > 0)
                        LOG_SEV(info) << "MAX RECEIVE TIME details, gid = " << b->gid << ", repair_time = "
                                                                            << b->repair_time;
                    if (b->whole_merge_tree_time > 0)
                        LOG_SEV(info) << "MAX RECEIVE TIME details, gid = " << b->gid << ", time_to_merge_trees = "
                                                                            << b->whole_merge_tree_time;

                    if (b->merge_call_time > 0)
                        LOG_SEV(info) << "MAX RECEIVE TIME details, gid = " << b->gid << ", merge_call_time = "
                                                                            << b->merge_call_time;
                    if (b->union_find_time > 0)
                        LOG_SEV(info) << "MAX RECEIVE TIME details, gid = " << b->gid << ", union_find_time = "
                                                                            << b->union_find_time;
                    if (b->sparsify_time > 0)
                        LOG_SEV(info) << "MAX RECEIVE TIME details, gid = " << b->gid << ", sparsify_time = "
                                                                            << b->sparsify_time;
                    if (b->expand_link_time > 0)
                        LOG_SEV(info) << "MAX RECEIVE TIME details, gid = " << b->gid << ", expand_link_time = "
                                                                            << b->expand_link_time;
                    if (b->is_done_time > 0)
                        LOG_SEV(info) << "MAX RECEIVE TIME details, gid = " << b->gid << ", is_done_time = "
                                                                            << b->is_done_time;
                });
            }
        }
#endif
        world.barrier();
    }

    if (read_plotfile)
    {
        amrex::Finalize();
    }

    return 0;
}
