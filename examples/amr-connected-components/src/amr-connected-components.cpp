//#define SEND_COMPONENTS

#include "reeber-real.h"

// to print nice backtrace on segfault signal
#include <signal.h>
#include <execinfo.h>
#include <cxxabi.h>

#include <AMReX.H>

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

#include "amr-plot-reader.h"
#include "amr-connected-components-complex.h"


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
using GidSet = Block::GidSet;

using TripletMergeTree = Block::TripletMergeTree;
using Neighbor = TripletMergeTree::Neighbor;

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

    ComponentDiagramsFunctor(Block* b, const LocalFunctor& lf) :
            block_(b),
            negate_(b->get_merge_tree().negate()),
            ignore_zero_persistence_(true),
            test_local(lf)
    {}

    void operator()(Neighbor from, Neighbor through, Neighbor to) const
    {
        if (!test_local(*block_, from))
        {
            return;
        }

        AmrVertexId current_vertex = from->vertex;

        Real birth_time = from->value;
        Real death_time = through->value;

        if (ignore_zero_persistence_ and birth_time == death_time)
        {
            return;
        }

        AmrVertexId root = block_->vertex_to_deepest_[current_vertex];
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
    {
        return false;
    }
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
    {
        throw std::runtime_error("Cannot read file " + infn);
    }

    if (ends_with(infn, ".npy"))
    {
        read_from_npy_file<DIM>(infn, world, nblocks, master_reader, assigner, header, domain);
    } else
    {
        if (split)
        {
            diy::io::split::read_blocks(infn, world, assigner, master_reader, header, FabBlockR::load);
        } else
        {
            diy::io::read_blocks(infn, world, assigner, master_reader, header, FabBlockR::load);
        }
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
    int min_cells = 10;

    using namespace opts;

    opts::Options ops(argc, argv);
    ops
            >> Option('b', "blocks", nblocks, "number of blocks to use")
            >> Option('m', "memory", in_memory, "maximum blocks to store in memory")
            >> Option('j', "jobs", threads, "threads to use during the computation")
            >> Option('s', "storage", prefix, "storage prefix")
            >> Option('i', "rho", rho, "iso threshold")
            >> Option('x', "mincells", min_cells, "minimal number of cells to output halo")
            >> Option('p', "profile", profile_path, "path to keep the execution profile")
            >> Option('l', "log", log_level, "log level");

    bool absolute =
            ops >> Present('a', "absolute", "use absolute values for thresholds (instead of multiples of mean)");
    bool read_plotfile = ops >> Present("plotfile", "read AMR plotfiles");
    bool negate = ops >> opts::Present('n', "negate", "sweep superlevel sets");
    // ignored for now, wrap is always assumed
    bool wrap = ops >> opts::Present('w', "wrap", "wrap");
    bool split = ops >> opts::Present("split", "use split IO");

    std::string input_filename, output_filename, output_diagrams_filename, output_integral_filename;


    std::vector<std::string> all_var_names { "density", "particle_mass_density", "xmom", "ymom", "zmom" };
    int n_mt_vars = 2;  // use sum of density + particle_mass_density

    if (ops >> Present('h', "help", "show help message") or
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

    bool write_diag = (ops >> PosOption(output_diagrams_filename));
    if (output_diagrams_filename == "none")
    {
        write_diag = false;
    }

    bool write_integral = (ops >> PosOption(output_integral_filename));
    if (output_integral_filename == "none")
    {
        write_integral = false;
    }

    diy::FileStorage storage(prefix);

    diy::Master master_reader(world, 1, in_memory, &FabBlockR::create, &FabBlockR::destroy);
    diy::Master master(world, threads, in_memory, &Block::create, &Block::destroy, &storage, &Block::save,
                       &Block::load);
    diy::ContiguousAssigner assigner(world.size(), nblocks);
    diy::MemoryBuffer header;
    diy::DiscreteBounds domain(DIM);

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

    read_plotfile = true;

    if (read_plotfile)
    {
        read_amr_plotfile(input_filename, all_var_names, n_mt_vars, world, nblocks, master_reader, header, domain);
    } else
    {
        read_from_file(input_filename, world, master_reader, assigner, header, domain, split, nblocks);
    }

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
                // TODO: refinment vector
                int local_ref = l->refinement()[0];
                int local_lev = l->level();

                master.add(cp.gid(),
                           new Block(b->fab, b->extra_names_, b->extra_fabs_,  local_ref, local_lev, domain, l->bounds(), l->core(), cp.gid(),
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

        world.barrier();
        LOG_SEV_IF(world.rank() == 0, info) << "Average = " << mean << ", rho = " << rho
                                                            << ", time to compute average: "
                                                            << dlog::clock_to_string(timer.elapsed());

        if (mean < 0 or std::isnan(mean) or std::isinf(mean) or mean > 1e+40)
        {
            LOG_SEV_IF(world.rank() == 0, error) << "Bad average = " << mean << ", do not proceed";
            if (read_plotfile)
                amrex::Finalize();
            return 1;
        }




        time_for_local_computation += timer.elapsed();
        dlog::flush();
        timer.restart();


        master.foreach([rho](Block* b, const diy::Master::ProxyWithLink& cp) {
            AMRLink* l = static_cast<AMRLink*>(cp.link());
            b->init(rho, l);
            cp.collectives()->clear();
        });


        world.barrier();
        LOG_SEV_IF(world.rank() == 0, info) <<
        "Time to initialize FabComponentBlocks (low vertices, local trees, components, outgoing edges): "
                << timer.elapsed();
        time_for_local_computation += timer.elapsed();
        dlog::flush();
        timer.restart();


    }

    int global_n_undone = 1;

    master.foreach(&send_edges_to_neighbors_cc<Real, DIM>);
    master.exchange();
    master.foreach(&delete_low_edges_cc<Real, DIM>);

    world.barrier();
    LOG_SEV_IF(world.rank() == 0, info)  << "edges symmetrized, time elapsed " << timer.elapsed();
    auto time_for_communication = timer.elapsed();
    dlog::flush();
    timer.restart();


    // debug: check symmetry
    master.foreach([](Block* b, const diy::Master::ProxyWithLink& cp) {
        auto* l = static_cast<AMRLink*>(cp.link());

        for(const diy::BlockID& receiver : link_unique(l, b->gid))
        {
            cp.enqueue(receiver, b->components_);
        }
    });

    master.exchange();

    master.foreach([](Block* b, const diy::Master::ProxyWithLink& cp) {
        auto* l = static_cast<AMRLink*>(cp.link());
        for(const diy::BlockID& sender : link_unique(l, b->gid))
        {
            std::vector<Component> received_components;
            cp.dequeue(sender, received_components);
            b->check_symmetry(sender.gid, received_components);
        }
    });

    LOG_SEV_IF(world.rank() == 0, info)  << "Symmetry checked in " << dlog::clock_to_string(timer.elapsed());
    time_for_communication += timer.elapsed();
    dlog::flush();
    timer.restart();
    // end symmetry checking

    int rounds = 0;
    while(global_n_undone)
    {
        rounds++;

        master.foreach(&amr_cc_send<Real, DIM>);
        master.exchange();
        master.foreach(&amr_cc_receive<Real, DIM>);

        world.barrier();
        LOG_SEV_IF(world.rank() == 0, info) << "MASTER round " << rounds << ", get OK";
        dlog::flush();
        master.exchange();
        //LOG_SEV_IF(world.rank() == 0, info) << "MASTER round " << rounds << ", collectives exchange OK";
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

    if (true)
    {
        master.foreach([](Block* b, const diy::Master::ProxyWithLink& cp) {
            size_t n_low = b->n_low_;
            size_t n_active = b->n_active_;
            size_t n_masked = b->n_masked_;
            cp.collectives()->clear();
            cp.all_reduce(n_low, std::plus<size_t>());
            cp.all_reduce(n_active, std::plus<size_t>());
            cp.all_reduce(n_masked, std::plus<size_t>());
        });

        master.exchange();

        world.barrier();

        const diy::Master::ProxyWithLink& proxy = master.proxy(master.loaded_block());

        size_t total_n_low = proxy.get<size_t>();
        size_t total_n_active = proxy.get<size_t>();
        size_t total_n_masked = proxy.get<size_t>();

        LOG_SEV_IF(world.rank() == 0, info) << "Total_n_low = " << total_n_low << ", total_n_active = "
                                                                << total_n_active << ", total_n_masked = " <<  total_n_masked;
        dlog::flush();
        world.barrier();
        timer.restart();
    }

#if 0
    auto time_for_output = timer.elapsed();
#else
    // save the result
    if (output_filename != "none")
    {
        if (!split)
        {
            diy::io::write_blocks(output_filename, world, master);
        } else
        {
            diy::io::split::write_blocks(output_filename, world, master);
        }
    }

    world.barrier();
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
                [&extra, &test_local, ignore_zero_persistence, rho](Block* b, const diy::Master::ProxyWithLink& cp) {
                    b->compute_final_connected_components();
                    output_persistence(b, cp, extra, test_local, rho, ignore_zero_persistence);
                });
    }

    world.barrier();
    LOG_SEV_IF(world.rank() == 0, info) << "Time to write diagrams:  " << dlog::clock_to_string(timer.elapsed());
    time_for_output += timer.elapsed();
    dlog::flush();
    timer.restart();

    if (write_integral)
    {
        master.foreach([](Block* b, const diy::Master::ProxyWithLink& cp)
        {
            b->compute_final_connected_components();
            b->compute_local_integral();

        });

        LOG_SEV_IF(world.rank() == 0, info) << "Local integrals computed";
        dlog::flush();
        world.barrier();

        master.foreach([output_integral_filename, domain, min_cells](Block* b, const diy::Master::ProxyWithLink& cp) {

            std::string integral_local_fname = fmt::format("{}-b{}.comp", output_integral_filename, b->gid);
            std::ofstream ofs(integral_local_fname);

            diy::Point<int, 3> domain_shape;
            for(int i = 0 ; i < 3; ++i)
            {
                domain_shape[i] = domain.max[i] - domain.min[i] + 1;
            }

            diy::GridRef<void*, 3> domain_box(nullptr, domain_shape, /* c_order = */ false);

            // local integral already stores number of vertices (set in init)
            // so we add it here just to print it
            b->extra_names_.insert(b->extra_names_.begin(), std::string("n_vertices"));

            for(const auto& root_values_pair : b->local_integral_)
            {
                AmrVertexId root = root_values_pair.first;
                if (root.gid != b->gid)
                    continue;

                auto& values = root_values_pair.second;

                if (values.count("n_vertices") == 0) {
                    fmt::print("ERROR HERE, no n_vertices, gid = {}\n", b->gid);
                }

                if (values.count("xmom") == 0) {
                    fmt::print("ERROR HERE, no xmom gid = {}\n", b->gid);
                }

                if (values.count("ymom") == 0) {
                    fmt::print("ERROR HERE, no ymom gid = {}\n", b->gid);
                }

                if (values.count("zmom") == 0) {
                    fmt::print("ERROR HERE, no zmom gid = {}\n", b->gid);
                }

                if (values.count("density") == 0) {
                    fmt::print("ERROR HERE, no density gid = {}\n", b->gid);
                }

                if (values.count("particle_mass_density") == 0) {
                    fmt::print("ERROR HERE, no particle_mass_density gid = {}\n", b->gid);
                }


                Real n_vertices = values.at("n_vertices");

                if (n_vertices < min_cells)
                    continue;

                Real vx = values.at("xmom") / n_vertices;
                Real vy = values.at("ymom") / n_vertices;
                Real vz = values.at("zmom") / n_vertices;

                Real m_gas = values.at("density");
                Real m_particles = values.at("particle_mass_density");
                Real m_total = m_gas + m_particles;

                fmt::print(ofs, "{} {} {} {} {} {} {} {} {}\n",
                        domain_box.index(b->local_.global_position(root)), // TODO: fix for non-flat AMR
                        n_vertices,
                        b->local_.global_position(root),
                        vx, vy, vz,
                        m_gas, m_particles, m_total);
            }
            ofs.close();
        });

        world.barrier();
        LOG_SEV_IF(world.rank() == 0, info) << "Time to compute and write integral:  "
                << dlog::clock_to_string(timer.elapsed());
        time_for_output += timer.elapsed();
        dlog::flush();
        timer.restart();
    }
#endif
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
    dlog::flush();

    world.barrier();

    std::string final_timings = fmt::format("read: {} local: {} exchange: {} output: {}\n", time_to_read_data,
                                            time_for_local_computation, time_for_communication, time_for_output);
    LOG_SEV_IF(world.rank() == 0, info) << final_timings;

    dlog::flush();

    if (false)
    {
        master.foreach([](Block* b, const diy::Master::ProxyWithLink& cp) {
            size_t n_low = b->n_low_;
            size_t n_active = b->n_active_;
            size_t n_masked = b->n_masked_;
            cp.collectives()->clear();
            cp.all_reduce(n_low, std::plus<size_t>());
            cp.all_reduce(n_active, std::plus<size_t>());
            cp.all_reduce(n_masked, std::plus<size_t>());
        });

        master.exchange();

        world.barrier();

        const diy::Master::ProxyWithLink& proxy = master.proxy(master.loaded_block());

        size_t total_n_low = proxy.get<size_t>();
        size_t total_n_active = proxy.get<size_t>();
        size_t total_n_masked = proxy.get<size_t>();

        LOG_SEV_IF(world.rank() == 0, info) << "Total_n_low = " << total_n_low << ", total_n_active = "
                                                                << total_n_active << ", total_n_masked = " <<  total_n_masked;
        dlog::flush();
        timer.restart();
    }

    if (read_plotfile)
    {
        amrex::Finalize();
    }

    return 0;
}
