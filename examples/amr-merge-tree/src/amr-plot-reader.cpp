#include <string>
#include <iostream>

#include "AMReX_ParmParse.H"
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_DataServices.H>


#include <iostream>
#include <stdexcept>
#include <vector>
#include <array>

//#include <Nyx.H>

#include <diy/master.hpp>
#include <diy/assigner.hpp>
#include <diy/link.hpp>
#include <diy/resolve.hpp>
#include <diy/fmt/format.h>
#include <diy/io/block.hpp> // for saving blocks in DIY format

#include "reeber/amr_helper.h"
#include "fab-block.h"

static constexpr unsigned DIY_DIM = 3;

using namespace amrex;

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

struct AMReXMeshHierarchy
{
/*
  A minimal class to describe the AMR hierarchy for analysis routines
 */
public:
    AMReXMeshHierarchy()
    {}

    void define(const AmrData& ad)
    {
        finestLevel = ad.FinestLevel();
        int nlevs = finestLevel + 1;
        ba.resize(nlevs);
        probSize = ad.ProbSize();
        probDomain = ad.ProbDomain();
        refRatio = ad.RefRatio();
        for(int lev = 0; lev < nlevs; ++lev)
        {
            ba[lev] = &ad.boxArray(lev);
        }
    }

    int FinestLevel() const
    { return finestLevel; }

    const BoxArray& boxArray(int level) const
    { return *ba[level]; }

    const Vector<int>& RefRatio() const
    { return refRatio; }

    const Vector<Real>& ProbSize() const
    { return probSize; }

    const Vector<Box>& ProbDomain() const
    { return probDomain; }

protected:
    int finestLevel;
    std::vector<const BoxArray*> ba;
    Vector<int> refRatio;
    Vector<Real> probSize;
    Vector<Box> probDomain;
};

struct AMReXDataHierarchy
{
/*
  Data on a AMReXMeshHierarchy, currently pointing to MultiFabs of
  named variables managed by an AmrData object.
*/
public:
    AMReXDataHierarchy(AmrData& ad, const Vector<std::string>& varNames)
    {
        mesh.define(ad);
        const Vector<std::string>& plotVarNames = ad.PlotVarNames();
        int nComp = varNames.size();
        int nlevs = mesh.FinestLevel() + 1;
        for(int i = 0; i < nComp; ++i)
        {
            int idx = -1;
            for(int j = 0; j < plotVarNames.size() && idx < 0; ++j)
            {
                if (plotVarNames[j] == varNames[i])
                { idx = j; }
            }
            if (ParallelDescriptor::IOProcessor() && idx < 0)
            {
                Abort("Cannot find variable=" + varNames[i] + " in pltfile");
            }
            std::vector<MultiFab*> mfs(nlevs);
            for(int lev = 0; lev < nlevs; ++lev)
            {
                mfs[lev] = &ad.GetGrids(lev, idx); // Note: This lazily triggers a MultiFab read in the AmrData
            }
            // Arnur: save idx in varMap
            varMap[varNames[i]] = std::make_pair(idx, mfs);
        }
    }

    MultiFab& GetGrids(int level, const std::string& name)
    {
        if (varMap.find(name) == varMap.end())
        {
            Abort("Unknown component requested");
        }
        return *(varMap[name].second[level]);
    }

    int GetIndex(const std::string& name) const
    {
        auto iter = varMap.find(name);
        if (iter != varMap.end())
            return iter->second.first;
        else
            return -1;
    }

    const AMReXMeshHierarchy& Mesh() const
    { return mesh; }

protected:
    AMReXMeshHierarchy mesh;
    std::map<std::string, std::pair<int, std::vector<MultiFab*>>> varMap;
};


void print_3v(const int* v)
{
    std::cout << "(" << v[0] << ", " << v[1] << ", " << v[2] << ")";
}

void print_box(const Box& b)
{
    std::cout << "Box(lo";
    print_3v(b.loVect());
    std::cout << ", hi";
    print_3v(b.hiVect());
    std::cout << ")" << std::endl;
}


void read_amr_plotfile(std::string infile,
                       std::string varName,
                       diy::mpi::communicator& world,
                       int nblocks,
                       diy::Master& master_reader,
//                       diy::ContiguousAssigner& assigner,
                       diy::MemoryBuffer& header,
                       diy::DiscreteBounds& domain_diy)
{
    amrex::Initialize(world);

    // Create the AmrData object from a pltfile on disk

    DataServices::SetBatchMode();
    Amrvis::FileType fileType(Amrvis::NEWPLT);
    DataServices* pdataServices =  new DataServices(infile, fileType);
    DataServices& dataServices =  *pdataServices;
    if (!dataServices.AmrDataOk())
    {
        DataServices::Dispatch(DataServices::ExitRequest, NULL);
    }
    AmrData& amrData = dataServices.AmrDataRef();

    amrex::Vector<string> varNames;
    varNames.push_back(varName);

    // Make a data struct for just the variables needed
    AMReXDataHierarchy* pdata = new AMReXDataHierarchy(amrData, varNames);
    AMReXDataHierarchy& data = *pdata;
    const AMReXMeshHierarchy& mesh = data.Mesh();

    int periodic = 0;
    std::array<bool, 3> is_periodic;
    for(int i = 0; i < 3; ++i)
        is_periodic[i] = 1; // periodic & (1 << i);

    const Box& domain = mesh.ProbDomain()[0];

//    std::cout << "printing ProbDomain" << std::endl;
//
//    for(const Box& d : mesh.ProbDomain())
//    {
//        print_box(d);
//    }
//
//    std::cout << "finisned printing ProbDomain" << std::endl;

    // Compute the volume integrals
    const int finestLevel = mesh.FinestLevel();
    const int nLev = finestLevel + 1;

    std::vector<int> gid_offsets = {0};
    std::vector<int> refinements = {1};
    for(int lev = 0; lev <= finestLevel; lev++)
    {

        if (world.rank() == 0)
        { fmt::print("Processing lev = {}", lev); }

        const MultiFab& mf = data.GetGrids(lev, varName);
        const BoxArray& ba = mf.boxArray();

        fmt::print("rank = {}, nblocks = {}, ba.size() = {}\n", world.rank(), nblocks, ba.size());

        nblocks += ba.size();
        gid_offsets.push_back(nblocks);

        auto refinement = mesh.RefRatio();

        if (world.rank() == 0)
        { fmt::print("refinement = {}", refinement[0]); }

        if (refinement.size() != 1)
            throw std::runtime_error("Unexpected uneven refinement");
        refinements.push_back(refinements.back() * refinement[0]);
    }

    for(int lev = 0; lev < nLev; ++lev)
    {
        const BoxArray ba = mesh.boxArray(lev);

        // Make boxes that are projection of finer ones (if exist)
        const BoxArray baf = lev < finestLevel
                             ? BoxArray(mesh.boxArray(lev + 1)).coarsen(mesh.RefRatio()[lev])
                             : BoxArray();

        // For each component listed...
        int componentIndex = data.GetIndex(varName);

        // Make copy of original data because we will modify here
        const MultiFab& mf = data.GetGrids(lev, varName);

        for(MFIter mfi(mf, false); mfi.isValid(); ++mfi)
        {
            const FArrayBox& myFab = mf[mfi];

            const Box& box = mfi.tilebox();
            Block::Shape shape;
            for(int i = 0; i < DIY_DIM; ++i)
                shape[i] = box.bigEnd()[i] - box.smallEnd()[i] + 1;

            // This is the Box on which the FArrayBox is defined.
            // Note that "abox" includes ghost cells (if there are any),
            // and is thus larger than or equal to "box".
            Box abox = myFab.box();

            int gid = gid_offsets[lev] + mfi.index();

            std::vector<std::pair<int, Box>> isects;

            diy::AMRLink* link = new diy::AMRLink(3, lev, refinements[lev], bounds(box), bounds(abox));
            // init fab
            // TODO: c_order!
            Real* fab_ptr = const_cast<Real*>(myFab.dataPtr(componentIndex));

            int fab_size = shape[0] * shape[1] * shape[2];
            Real total_sum = 0;
            for(int i = 0 ; i < fab_size; ++i)
            {
//                if (gid == 0) { fmt::print("rank = {}, gid = {}, value = {}, i  = {}, size = {}\n", world.rank(), gid, fab_ptr[i], i, sizeof(fab_ptr[i])); }
                total_sum += fab_ptr[i];
            }
            fmt::print("rank = {}, gid = {}, fab_ptr = {}, sum = {}, fabs_soize = {}\n",
                    world.rank(), gid, (void*)fab_ptr, total_sum, fab_size);

//            Real* fab_ptr = myFab.dataPtr(componentIndex);
            master_reader.add(gid, new Block(fab_ptr, shape), link);
            fmt::print("rank = {}, ADDED\n", world.rank());

            fmt::print(
                    "rank = {}, gid = {},  smallEnd = ({}, {}, {}), bigEnd = ({}, {}, {}), mfi.index - {}\n",
                    world.rank(), gid,  box.smallEnd()[0], box.smallEnd()[1], box.smallEnd()[2],
                    box.bigEnd()[0], box.bigEnd()[1], box.bigEnd()[2], mfi.index());


            // record wrap
            for(int dir_x : {-1, 0, 1})
            {
                if (!is_periodic[0] && dir_x) continue;
                if (dir_x < 0 && box.loVect()[0] != domain.loVect()[0]) continue;
                if (dir_x > 0 && box.hiVect()[0] != domain.hiVect()[0]) continue;

                for(int dir_y : {-1, 0, 1})
                {
                    if (!is_periodic[1] && dir_y) continue;
                    if (dir_y < 0 && box.loVect()[1] != domain.loVect()[1]) continue;
                    if (dir_y > 0 && box.hiVect()[1] != domain.hiVect()[1]) continue;

                    for(int dir_z : {-1, 0, 1})
                    {
                        if (dir_x == 0 && dir_y == 0 && dir_z == 0)
                            continue;

                        if (!is_periodic[2] && dir_z) continue;
                        if (dir_z < 0 && box.loVect()[2] != domain.loVect()[2]) continue;
                        if (dir_z > 0 && box.hiVect()[2] != domain.hiVect()[2]) continue;

                        link->add_wrap(diy::Direction{dir_x, dir_y, dir_z});
                    }
                }
            }

            // record neighbors
            for(int nbr_lev = std::max(0, lev - 1); nbr_lev <= std::min(finestLevel, lev + 1); ++nbr_lev)
            {

//                std::cout << "in nbr_lev loop, nbr_lev = " << nbr_lev << std::endl;
                // gotta do this yoga to work around AMReX's static variables
                const Box& nbr_lev_domain = mesh.ProbDomain().at(nbr_lev);
                Periodicity periodicity(IntVect(AMREX_D_DECL(nbr_lev_domain.length(0) * is_periodic[0],
                                                             nbr_lev_domain.length(1) * is_periodic[1],
                                                             nbr_lev_domain.length(2) * is_periodic[2])));


                const std::vector<IntVect>& pshifts = periodicity.shiftIntVect();

                // TODO: here we always assume ghosts, get this information somehow
                int ng = 1;
                const BoxArray& ba = mesh.boxArray(nbr_lev);
                // TODO: need to divide?
                int ratio = mesh.RefRatio().at(nbr_lev);

                Box gbx = box;

                if (nbr_lev < lev)
                    gbx.coarsen(ratio);
                else if (nbr_lev > lev)
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
                        fmt::print("{}: gid = {}, adding neighbor gid = {}\n", world.rank(), gid, nbr_gid);
                        const Box& nbr_box = ba[is.first];
                        Box nbr_ghost_box = grow(nbr_box, ng);

                        link->add_neighbor(diy::BlockID{nbr_gid,
                                                        -1});        // we don't know the proc, but we'll figure it out later through DynamicAssigner
                        link->add_bounds(nbr_lev, refinements[nbr_lev], bounds(nbr_box), bounds(nbr_ghost_box));
                    }
                }
            }
        }
    }

    fmt::print("{}: started fixing links\n", world.rank());
    // fill dynamic assigner and fix links
    diy::DynamicAssigner assigner(master_reader.communicator(), master_reader.communicator().size(), nblocks);
    diy::fix_links(master_reader, assigner);

    master_reader.foreach([](Block* b, const diy::Master::ProxyWithLink& cp) {
        auto* l = static_cast<diy::AMRLink*>(cp.link());

        auto receivers = link_unique(l, cp.gid());

        fmt::print("{}: level = {}, shape = {}, core = {} - {}, bounds = {} - {}, neighbors = {}\n",
                   cp.gid(), l->level(), b->fab.shape(),
                   l->core().min, l->core().max,
                   l->bounds().min, l->bounds().max,
                   container_to_string(receivers));
    });
}
