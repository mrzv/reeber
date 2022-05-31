#include <string>
#include <iostream>

#include "AMReX_ParmParse.H"
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_PlotFileUtil.H>

#include <diy/master.hpp>
#include <diy/assigner.hpp>
#include <diy/link.hpp>
#include <diy/resolve.hpp>
#include <diy/io/block.hpp> // for saving blocks in DIY format

#include <reeber/format.h>
#include "reeber/amr_helper.h"
#include "fab-block.h"

static constexpr unsigned DIY_DIM = 3;

using namespace amrex;
// by this we also use amrex_real, not reeber_real


using Block = FabBlock<Real, DIY_DIM>;

diy::AMRLink::Bounds bounds(const amrex::Box& box)
{
    diy::AMRLink::Bounds bounds(3);
    for(int i = 0; i < 3; ++i)
    {
        bounds.min[i] = box.loVect()[i];
        bounds.max[i] = box.hiVect()[i];
    }
    return bounds;
}


void read_amr_plotfile(std::string infile,
        std::vector<std::string> all_var_names,
        int n_mt_vars,
        diy::mpi::communicator& world,
        int nblocks,
        diy::Master& master_reader,
        diy::MemoryBuffer& header,
        Real& cell_volume,
        diy::DiscreteBounds& domain_diy)
{
    amrex::Initialize(world);

    bool debug = false;

    PlotFileData plotfile(infile);
    const auto& pf_comp_names = plotfile.varNames();
    const int n_levels = plotfile.finestLevel() + 1;
    const int finest_level = plotfile.finestLevel();

    cell_volume = 1.0;
    for(size_t i = 0; i < DIY_DIM; ++i)
    {
        cell_volume *= plotfile.cellSize(0)[i];
    }

    // TODO: fix wrap
    int periodic = 0;
    std::array<bool, DIY_DIM> is_periodic;
    for(size_t i = 0; i < DIY_DIM; ++i)
        is_periodic[i] = 1; // periodic & (1 << i);

    const Box& domain = plotfile.probDomain(0);

    for(size_t i = 0; i < DIY_DIM; ++i)
    {
        domain_diy.min[i] = domain.loVect()[i];
        domain_diy.max[i] = domain.hiVect()[i];
    }

    std::vector<int> gid_offsets = {0};
    std::vector<int> refinements = {1};
    nblocks = 0;

    // iterate over all levels to collect refinements and box array sizes (=gid offsets)
    for(int level = 0; level < n_levels; ++level)
    {
        const MultiFab& mf = plotfile.get(level, all_var_names[0]);

        BoxArray ba = mf.boxArray();
        nblocks += ba.size();
        gid_offsets.push_back(nblocks);

        refinements.push_back(refinements.back() * plotfile.refRatio(level));
    }

    Real* fab_ptr_copy{nullptr};

    Real total_sum = 0;
    Real total_sum_wo = 0;
    long int n_nans = 0;
    long int n_infs = 0;
    long int n_negs = 0;
    long int n_wo = 0;

    Real total_sum_1 = 0;
    long int n_nans_1 = 0;
    long int n_infs_1 = 0;
    long int n_negs_1 = 0;
    long int n_wo_1 = 0;
    long int n_additions = 0;

    long long int total_fab_vertices = 0;
    long long int total_a_vertices = 0;

    std::map<int, Real*> gid_to_fab;

    std::map<int, std::vector<Real*>> gid_to_extra_pointers;

    for(size_t var_idx = 0; var_idx < all_var_names.size(); ++var_idx)
    {
        for(int level = 0; level < n_levels; ++level)
        {
            const MultiFab& mf = plotfile.get(level, all_var_names[var_idx]);
            const BoxArray ba = mf.boxArray();
            BoxArray ba_finer;
            if (level < finest_level)
            {
                // TODO: this will load actual data; we only boxes from finer level
                const MultiFab& mf_finer = plotfile.get(level + 1, all_var_names[var_idx]);
                ba_finer = mf_finer.boxArray();
            }

            // false is for no tiling in MFIter; we want boxes exactly as they are in plotfile
            for(MFIter mfi(mf, false); mfi.isValid(); ++mfi)
            {
                const FArrayBox& my_fab = mf[mfi];
                Block::Shape valid_shape, a_shape;
                const Box& valid_box = mfi.validbox();
                for(size_t i = 0; i < DIY_DIM; ++i)
                    valid_shape[i] = valid_box.bigEnd()[i] - valid_box.smallEnd()[i] + 1;

                total_fab_vertices += valid_shape[0] * valid_shape[1] * valid_shape[2];

                // This is the Box on which the FArrayBox is defined.
                // Note that "abox" includes ghost cells (if there are any),
                // and is thus larger than or equal to "box".
                Box abox = my_fab.box();
                for(size_t i = 0; i < DIY_DIM; ++i)
                    a_shape[i] = abox.bigEnd()[i] - abox.smallEnd()[i] + 1;
                total_a_vertices += a_shape[0] * a_shape[1] * a_shape[2];

                if (a_shape != valid_shape)
                {
                    throw std::runtime_error("ghosts in plotfile - not expected");
                }

                int gid = gid_offsets[level] + mfi.index();
                if (debug) fmt::print( "amr-plot-reader: ALL SHAPES rank = {}, gid = {}; fab shape = ({}, {}, {}), a_shape = ({}, {}, {})\n", world.rank(), gid, valid_shape[0], valid_shape[1], valid_shape[2], a_shape[0], a_shape[1], a_shape[2]);
                if (debug) fmt::print( "amr-plot-reader: ALL BOXES rank = {}, gid = {}; FABBOX smallEnd = ({}, {}, {}), bigEnd = ({}, {}, {}),  ABOX  smallEnd = ({}, {}, {}), bigEnd = ({}, {}, {}),\n", world.rank(), gid,
                        valid_box.smallEnd()[0], valid_box.smallEnd()[1], valid_box.smallEnd()[2], valid_box.bigEnd()[0], valid_box.bigEnd()[1], valid_box.bigEnd()[2],
                        abox.smallEnd()[0], abox.smallEnd()[1], abox.smallEnd()[2], abox.bigEnd()[0], abox.bigEnd()[1], abox.bigEnd()[2]);


                std::vector<std::pair<int, Box>> isects;
                diy::AMRLink* link = new diy::AMRLink(3, level, refinements[level], bounds(valid_box), bounds(abox));
                // init fab
                // TODO: c_order?
                if (debug) fmt::print("n_comps = {}, var_idx = {}, smallestPtr = {}, nextPtr = {}, diff = {}, max_ptr = {}\n", my_fab.nComp(), var_idx,
                        (void*) my_fab.dataPtr(0), (void*) my_fab.dataPtr(1), (my_fab.dataPtr(1) - my_fab.dataPtr(0)), (void*) my_fab.dataPtr(my_fab.nComp() - 1));

                Real* fab_ptr = const_cast<Real*>(my_fab.dataPtr(0));
                long long int fab_size = a_shape[0] * a_shape[1] * a_shape[2];
                if (var_idx == 0)
                {
                    fab_ptr_copy = new Real[fab_size];
                    gid_to_fab[gid] = fab_ptr_copy;
                    memcpy(fab_ptr_copy, fab_ptr, sizeof(Real) * fab_size);

                    // allocate memory for all fields that we store in FabBlock
                    // actual copying for next fields will happen later
                    std::vector<Real*> extra_pointers;
                    for(size_t i = 0; i < all_var_names.size(); ++i)
                    {
                        Real* extra_ptr_copy = new Real[fab_size];
                        if (debug) fmt::print("amr-plot-reader: gid = {}, size of Real = {}, allocated memory for {}\n", gid, sizeof(Real), all_var_names[var_idx]);
                        extra_pointers.push_back(extra_ptr_copy);
                    }

                    gid_to_extra_pointers[gid] = extra_pointers;
                    memcpy(extra_pointers[0], fab_ptr, sizeof(Real) * fab_size);

                    for(int i = 0; i < fab_size; ++i)
                    {
                        total_sum += fab_ptr[i];
                        n_nans += std::isnan(fab_ptr[i]);
                        n_infs += std::isinf(fab_ptr[i]);
                        n_negs += (fab_ptr[i] < 0);
                        if (not std::isnan(fab_ptr[i]) and not std::isinf(fab_ptr[i]))
                        {
                            total_sum_wo += fab_ptr[i];
                            n_wo += 1;
                        }
                    }
                    if (debug) { fmt::print("FIELD 0 rank = {}, gid = {}, sum = {}, fabs_size = {}, avg_in_fab = {}, n_nans = {}, n_infs = {}, n_negs = {}, n_wo = {}, avg_wo = {}\n", world.rank(), gid, total_sum, fab_size, total_sum / fab_size, n_nans, n_infs, n_negs, n_wo, total_sum_wo / n_wo); }

                    master_reader.add(gid, new Block(fab_ptr_copy, all_var_names, extra_pointers, a_shape), link);

                    // record wrap
                    for(int dir_x : {-1, 0, 1})
                    {
                        if (!is_periodic[0] && dir_x) continue;
                        if (dir_x < 0 && valid_box.loVect()[0] != domain.loVect()[0]) continue;
                        if (dir_x > 0 && valid_box.hiVect()[0] != domain.hiVect()[0]) continue;

                        for(int dir_y : {-1, 0, 1})
                        {
                            if (!is_periodic[1] && dir_y) continue;
                            if (dir_y < 0 && valid_box.loVect()[1] != domain.loVect()[1]) continue;
                            if (dir_y > 0 && valid_box.hiVect()[1] != domain.hiVect()[1]) continue;
                            for(int dir_z : {-1, 0, 1})
                            {
                                if (dir_x == 0 && dir_y == 0 && dir_z == 0)
                                    continue;

                                if (!is_periodic[2] && dir_z) continue;
                                if (dir_z < 0 && valid_box.loVect()[2] != domain.loVect()[2]) continue;
                                if (dir_z > 0 && valid_box.hiVect()[2] != domain.hiVect()[2]) continue;

                                link->add_wrap(diy::Direction{dir_x, dir_y, dir_z});
                            }
                        }
                    }
                    // record neighbors
                    for(int nbr_lev = std::max(0, level - 1); nbr_lev <= std::min(finest_level, level + 1); ++nbr_lev)
                    {
                        // gotta do this yoga to work around AMReX's static variables
                        const Box& nbr_lev_domain = plotfile.probDomain(nbr_lev);
                        Periodicity periodicity(IntVect(AMREX_D_DECL(nbr_lev_domain.length(0) * is_periodic[0],
                                nbr_lev_domain.length(1) * is_periodic[1],
                                nbr_lev_domain.length(2) * is_periodic[2])));

                        const std::vector<IntVect>& pshifts = periodicity.shiftIntVect();
                        // TODO: here we always assume ghosts, get this information somehow
                        int ng = 0;
                        const BoxArray& ba = plotfile.boxArray(nbr_lev);

//                        fmt::print("gid = {}, level = {}, refRatio = {}\n", gid, level, plotfile.refRatio(level));

                        // TODO!
//                        int ratio = mesh.RefRatio().at(std::min(lev, nbr_lev));
                        int ratio = plotfile.refRatio(level);

                        if (ratio == 0 and level != finest_level)
                            throw std::runtime_error("ration!");

                        if (ratio == 0)
                            ratio = plotfile.refRatio(level - 1);

                        Box gbx = valid_box;
                        if (nbr_lev < level)
                            gbx.coarsen(ratio);
                        else if (nbr_lev > level)
                            gbx.refine(ratio);
                        gbx.grow(1);

                        for(const auto& piv : pshifts)
                        {
                            ba.intersections(gbx + piv, isects);
                            for(const auto& is : isects)
                            {
                                // is.first is the index of neighbor box
                                // ba[is.first] is the neighbor box
                                int nbr_gid = gid_offsets[nbr_lev] + is.first;
                                const Box& nbr_box = ba[is.first];
                                Box nbr_ghost_box = grow(nbr_box, ng);
                                if (debug) { fmt::print("nbr_box = {}, nbr_lev = {}, my box = {}, my level = {}, ratio = {}\n", nbr_box, nbr_lev, valid_box, level, ratio); }
                                link->add_neighbor(diy::BlockID{nbr_gid,
                                                                -1});        // we don't know the proc, but we'll figure it out later through DynamicAssigner
                                // ghosts not expected, hence nbr_box in the 2 last parameter
                                link->add_bounds(nbr_lev, refinements[nbr_lev], bounds(nbr_box), bounds(nbr_box));
                            }
                        }
                    }
                } else
                {
                    Real* block_extra_ptr = gid_to_extra_pointers.at(gid).at(var_idx);

                    Real* block_fab_ptr = gid_to_fab.at(gid);
                    bool add_to_fab = (static_cast<decltype(n_mt_vars)>(var_idx) < n_mt_vars);
                    if (debug) fmt::print("Adding next field, block_fab_ptr = {}, fab_ptr = {}, gid = {}, fab_size = {}\n", (void*) block_fab_ptr, (void*) fab_ptr, gid, fab_size);
                    for(int i = 0; i < fab_size; ++i)
                    {
                        total_sum_1 += fab_ptr[i];
                        n_nans_1 += std::isnan(fab_ptr[i]);
                        n_infs_1 += std::isinf(fab_ptr[i]);
                        n_negs_1 += (fab_ptr[i] < 0);
                        if (not std::isnan(fab_ptr[i]) and not std::isinf(fab_ptr[i]))
                        {
                            n_wo_1 += 1;
                        }

                        if (add_to_fab)
                        {
                            block_fab_ptr[i] += fab_ptr[i];
                            n_additions += 1;
                        }

                        block_extra_ptr[i] = fab_ptr[i];
                    }
                    if (debug) fmt::print( "Added next field, block_fab_ptr = {}, fab_ptr = {}, gid = {}, n_nans_1 = {}, n_negs_1 = {}, n_infs_1 = {}, totao_sum_1 = {}\n", (void*) block_fab_ptr, (void*) fab_ptr, gid, n_nans_1, n_negs_1, n_infs_1, total_sum_1);
                }
            } // loop over tiles
        } // loop over all_var_names
    } // loop over levels




    if (debug) fmt::print("ADDED n_additions = {}\n", n_additions);
    if (debug) fmt::print("SHAPE_TOT  total_a_vertices = {}\n", total_a_vertices);
    // fill dynamic assigner and fix links
    diy::DynamicAssigner assigner(master_reader.communicator(), master_reader.communicator().size(), nblocks);
    diy::fix_links(master_reader, assigner);

#ifdef REEBER_COMPUTE_GAS_VELOCITIES

    // find index of gas density
    auto d_iter = std::find(all_var_names.begin(), all_var_names.end(), "density");
    if (d_iter == all_var_names.end())
    {
        return;
    }
    auto d_idx = d_iter - all_var_names.begin();

    // divide momenta by density to get velocities
    master_reader.foreach([d_idx](Block* b, const diy::Master::ProxyWithLink& cp)
    {
        for(size_t i = 0; i < b->extra_names_.size(); ++i)
        {
            if (b->extra_names_[i] == "xmom" or b->extra_names_[i] == "ymom" or b->extra_names_[i] == "zmom")
            {
                Real* momentum_ptr = b->extra_fabs_[i].data();
                Real* density_ptr = b->extra_fabs_[d_idx].data();
                for(size_t j = 0; j < b->extra_fabs_[i].size(); ++j)
                {
                    if (density_ptr[j] > 0)
                    {
                        momentum_ptr[j] /= density_ptr[j];
                    }
                }
                std::string vel_name = "gas_v_";
                vel_name.append(b->extra_names_[i], 0, 1);
                b->extra_names_[i] = vel_name;
                amrex::Print() << "Velocities computed, names = " << container_to_string(b->extra_names_) << std::endl;
            }
        }
    });

#endif



    if (debug)
    {
        master_reader.foreach([debug](Block* b, const diy::Master::ProxyWithLink& cp)
        {
            auto* l = static_cast<diy::AMRLink*>(cp.link());
            auto receivers = link_unique(l, cp.gid());
            fmt::print("{}: level = {}, shape = {}, core = {} - {}, bounds = {} - {}, neighbors = {}\n", cp.gid(), l->level(), b->fab.shape(), l->core().min, l->core().max, l->bounds().min, l->bounds().max, container_to_string(receivers));
        });
    }
}
