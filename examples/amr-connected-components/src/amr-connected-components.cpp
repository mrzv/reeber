//#define SEND_COMPONENTS

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

#include "../../amr-merge-tree/include/fab-block.h"
#include "fab-cc-block.h"
#include "reader-interfaces.h"

#include "../../local-global/output-persistence.h"

#include "../../amr-merge-tree/include/read-npy.h"

#include "amr-connected-components-simple.h"


// block-independent types
using AMRLink = diy::AMRLink;

using Bounds = diy::DiscreteBounds;
using AmrVertexId = r::AmrVertexId;
using AmrEdge = reeber::AmrEdge;

#define DIM 3

using FabBlockR = FabBlock<Real, DIM>;

using Block = FabComponentBlock<Real, DIM>;
using Vertex = Block::Vertex;
using Component = Block::Component;
using MaskedBox = Block::MaskedBox;
using GidVector = Block::GidVector;
using GidContainer = Block::GidContainer;

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


int main(int argc, char** argv)
{
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
    Real theta = 90.0;

    using namespace opts;

    opts::Options ops(argc, argv);
    ops
            >> Option('b', "blocks", nblocks, "number of blocks to use")
            >> Option('m', "memory", in_memory, "maximum blocks to store in memory")
            >> Option('j', "jobs", threads, "threads to use during the computation")
            >> Option('s', "storage", prefix, "storage prefix")
            >> Option('i', "rho", rho, "iso threshold")
            >> Option('x', "theta", theta, "integral threshold")
            >> Option('p', "profile", profile_path, "path to keep the execution profile")
            >> Option('l', "log", log_level, "log level");

    bool absolute =
            ops >> Present('a', "absolute", "use absolute values for thresholds (instead of multiples of mean)");
    bool negate = ops >> opts::Present('n', "negate", "sweep superlevel sets");
    // ignored for now, wrap is always assumed
    bool wrap = ops >> opts::Present('w', "wrap", "wrap");
    bool split = ops >> opts::Present("split", "use split IO");
    bool write_halo_diagrams = ops >> opts::Present("write-halo-diagrams", "Write diagrams of each halo");

    std::string input_filename, output_integral_filename;

    if (ops >> Present('h', "help", "show help message") or
        not(ops >> PosOption(input_filename))
        or not(ops >> PosOption(output_integral_filename)))
    {
        if (world.rank() == 0)
        {
            fmt::print("Usage: {} INPUT.AMR OUTPUT-INTEGRAL \n", argv[0]);
            fmt::print("Compute persistent integral from AMR data\n");
            fmt::print("{}", ops);
        }
        return 1;
    }

    bool write_integral = not (output_integral_filename == "none");

    if (write_integral)
    {
        if ((negate and theta < rho) or (not negate and theta > rho))
            throw std::runtime_error("Bad integral threshold");
    }

    diy::FileStorage storage(prefix);

    diy::Master master_reader(world, 1, in_memory, &FabBlockR::create, &FabBlockR::destroy);
    diy::Master master(world, threads, in_memory, &Block::create, &Block::destroy, &storage, &Block::save,
                       &Block::load);
    diy::ContiguousAssigner assigner(world.size(), nblocks);
    diy::MemoryBuffer header;
    diy::DiscreteBounds domain;

    dlog::add_stream(std::cerr, dlog::severity(log_level))
            << dlog::stamp() << dlog::aux_reporter(world.rank()) << dlog::color_pre() << dlog::level()
            << dlog::color_post() >> dlog::flush();

    world.barrier();
    dlog::Timer timer;
    LOG_SEV_IF(world.rank() == 0, info) << "Starting computation, input_filename = " << input_filename << ", nblocks = "
                                                                                     << nblocks
                                                                                     << ", rho = " << rho;
    dlog::flush();
    world.barrier();

    read_from_file(input_filename, world, master_reader, assigner, header, domain, split, nblocks);

    world.barrier();

    auto time_to_read_data = timer.elapsed();
    LOG_SEV_IF(world.rank() == 0, info) << "Data read, local size = " << master.size();
    LOG_SEV_IF(world.rank() == 0, info) << "Time to read data:       " << dlog::clock_to_string(timer.elapsed());
    dlog::flush();
    timer.restart();

    world.barrier();

    // copy FabBlocks to FabComponentBlocks
    // in FabTmtConstructor mask will be set and local trees will be computed
    // FabBlock can be safely discarded afterwards

    master_reader.foreach(
            [&master, domain, rho, negate, absolute](FabBlockR* b, const diy::Master::ProxyWithLink& cp) {
                auto* l = static_cast<AMRLink*>(cp.link());
                AMRLink* new_link = new AMRLink(*l);

                // prepare neighbor box info to save in MaskedBox
                int local_ref = l->refinement();
                int local_lev = l->level();

                master.add(cp.gid(),
                           new Block(b->fab, local_ref, local_lev, domain, l->bounds(), l->core(), cp.gid(),
                                     new_link, rho, negate, absolute),
                           new_link);

            });


    auto time_for_local_computation = timer.elapsed();

    Real mean = std::numeric_limits<Real>::min();

    if (absolute)
    {
        LOG_SEV_IF(world.rank() == 0, info) << "Time to compute local trees and components:  "
                << dlog::clock_to_string(timer.elapsed());
        dlog::flush();
        timer.restart();
    } else
    {
        LOG_SEV_IF(world.rank() == 0, info) << "Time to construct FabComponentBlocks: "
                << dlog::clock_to_string(timer.elapsed());
        dlog::flush();
        timer.restart();

        master.foreach([](Block* b, const diy::Master::ProxyWithLink& cp) {
            cp.collectives()->clear();
            cp.all_reduce(b->sum_, std::plus<Real>());
            cp.all_reduce(static_cast<Real>(b->n_unmasked_) * b->scaling_factor(), std::plus<Real>());
        });

        master.exchange();

        const diy::Master::ProxyWithLink& proxy = master.proxy(master.loaded_block());

        mean = proxy.get<Real>() / proxy.get<Real>();
        rho *= mean;                                            // now rho contains absolute threshold
        theta *= mean;

        world.barrier();
        LOG_SEV_IF(world.rank() == 0, info) << "Average = " << mean << ", rho = " << rho
                                                            << ", time to compute average: "
                                                            << dlog::clock_to_string(timer.elapsed());

        time_for_local_computation += timer.elapsed();
        dlog::flush();
        timer.restart();

        master.foreach([rho](Block* b, const diy::Master::ProxyWithLink& cp) {
            AMRLink* l = static_cast<AMRLink*>(cp.link());
            b->init(rho, l);
            cp.collectives()->clear();
        });


        world.barrier();
        LOG_SEV_IF(world.rank() == 0, info) << "Time to initialize FabComponentBlocks (low vertices, local trees, components, outgoing edges): " << timer.elapsed();
        time_for_local_computation += timer.elapsed();
        dlog::flush();
        timer.restart();
    }

    int global_n_undone = 1;

    master.foreach(&send_edges_to_neighbors<Block>);
    master.exchange();
    master.foreach(&delete_low_edges<Block>);

    world.barrier();
    LOG_SEV_IF(world.rank() == 0, info)  << "edges symmetrized, time elapsed " << timer.elapsed();
    auto time_for_communication = timer.elapsed();
    dlog::flush();
    timer.restart();

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

    int rounds = 0;
    while (global_n_undone)
    {
        rounds++;

        master.foreach(&amr_tmt_send<Real, DIM>);
        master.exchange();
        master.foreach(&amr_tmt_receive<Real, DIM>);

        LOG_SEV_IF(world.rank() == 0, info) << "MASTER round " << rounds << ", get OK";
        dlog::flush();
        master.exchange();
        // to compute total number of undone blocks
        global_n_undone = master.proxy(master.loaded_block()).read<int>();
        LOG_SEV_IF(world.rank() == 0, info) << "MASTER round " << rounds << ", global_n_undone = " << global_n_undone;
        dlog::flush();
    }

    world.barrier();

    //    fmt::print("world.rank = {}, time for exchange = {}\n", world.rank(), dlog::clock_to_string(timer.elapsed()));

    LOG_SEV_IF(world.rank() == 0, info) << "Time for exchange:  " << dlog::clock_to_string(timer.elapsed());
    time_for_communication += timer.elapsed();
    dlog::flush();
    timer.restart();

    //    fmt::print("----------------------------------------\n");
    //    master.foreach([](Block* b, const diy::Master::ProxyWithLink& cp) {
    //        int my_nodes = 0;
    //        int not_my_nodes = 0;
    //        const auto& const_tree = b->mt;
    //        for(const auto& x : const_tree.nodes()) {
    //            if (x.second->vertex != x.first)
    //                fmt::print("Bad node vertex = {}, node vertex = {}\n", x.first, x.second->vertex);
    //            if (b->gid == x.second->vertex.gid) {
    //                my_nodes++;
    //            } else {
    //                not_my_nodes++;
    //            }
    //        }
    //        fmt::print("finished gid = {}, mt.size = {}, my_nodes = {}, not_my_nodes = {}\n", b->gid, b->mt.size(), my_nodes, not_my_nodes);
    //    });
    //    fmt::print("----------------------------------------\n");

    // save the result
//    if (output_filename != "none")
//    {
//        if (!split)
//            diy::io::write_blocks(output_filename, world, master);
//        else
//            diy::io::split::write_blocks(output_filename, world, master);
//    }

//    world.barrier();
//    LOG_SEV_IF(world.rank() == 0, info) << "Time to write tree:  " << dlog::clock_to_string(timer.elapsed());
    auto time_for_output = timer.elapsed();
//    dlog::flush();
//    timer.restart();

//    bool verbose = false;
//
//    world.barrier();
//    LOG_SEV_IF(world.rank() == 0, info) << "Time to write diagrams:  " << dlog::clock_to_string(timer.elapsed());
//    time_for_output += timer.elapsed();
//    dlog::flush();
//    timer.restart();

    if (write_integral)
    {
        master.foreach([rho, theta](Block* b, const diy::Master::ProxyWithLink& cp) {

            b->compute_local_integral(rho, theta);

            AMRLink* l = static_cast<AMRLink*>(cp.link());

            for(const auto& vertex_value_pair : b->local_integral_)
            {
                int receiver_gid = vertex_value_pair.first.gid;
                if (receiver_gid == b->gid)
                    continue;
                auto receiver = l->target(l->find(receiver_gid));
                cp.enqueue(receiver, vertex_value_pair);
            }
        });

        master.exchange();

        diy::io::SharedOutFile integral_file(output_integral_filename, world);

        master.foreach([&integral_file](Block* b, const diy::Master::ProxyWithLink& cp) {
            AMRLink* l = static_cast<AMRLink*>(cp.link());
            for(auto bid : l->neighbors())
            {
                while (cp.incoming(bid.gid))
                {
                    Block::LocalIntegral::value_type x;
                    cp.dequeue(bid, x);
                    assert(x.first.gid == b->gid);
                    b->local_integral_[x.first] += x.second;
                }
            }

            for(const auto& root_value_pair : b->local_integral_)
            {
                AmrVertexId root = root_value_pair.first;
                if (root.gid != b->gid)
                    continue;

                integral_file << fmt::format("{} {}\n", b->local_.global_position(root), root_value_pair.second);
            }
        });

        world.barrier();
        LOG_SEV_IF(world.rank() == 0, info) << "Time to compute and write integral:  "
                << dlog::clock_to_string(timer.elapsed());
        time_for_output += timer.elapsed();
        dlog::flush();
        timer.restart();
    }

//    master.foreach([](Block* b, const diy::Master::ProxyWithLink& cp) {
//        auto sum_n_vertices_pair = b->get_local_stats();
//        cp.collectives()->clear();
//        cp.all_reduce(sum_n_vertices_pair.first, std::plus<Real>());
//        cp.all_reduce(sum_n_vertices_pair.second, std::plus<size_t>());
//
//    });
//
//    master.exchange();
//
//    world.barrier();
//
//    const diy::Master::ProxyWithLink& proxy = master.proxy(master.loaded_block());
//
//    Real total_value = proxy.get<Real>();
//    size_t total_vertices = proxy.get<size_t>();
//
//    LOG_SEV_IF(world.rank() == 0, info) << "Total value = " << total_value << ", total # vertices = " << total_vertices
//                                                            << ", mean = " << mean;
//    dlog::flush();
//
//    std::string final_timings = fmt::format("read: {} local: {} exchange: {} output: {}\n", time_to_read_data,
//                                            time_for_local_computation, time_for_communication, time_for_output);
//    LOG_SEV_IF(world.rank() == 0, info) << final_timings;
//    dlog::flush();


    return 0;
}
