#ifndef REEBER_IO_BOXLIB_H
#define REEBER_IO_BOXLIB_H

// STL
#include <cassert>
#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>

// DIY
#include <diy/types.hpp>
#include <diy/mpi.hpp>
#include <diy/serialization.hpp>

// Logging
#include <dlog/log.h>

// Reeber
#include <reeber/grid.h>

// BoxLib
#include "ParallelDescriptor.H"
#include "DataServices.H"
#include "Geometry.H"
#include "Utility.H"

namespace reeber
{

namespace io
{

namespace BoxLib
{
    // Sets up BoxLib and puts DataServices into batch mode
    struct environment
    {
        environment()                                                    { int argc = 0; char** argv; ::BoxLib::Initialize(argc, argv, false); DataServices::SetBatchMode(); }
        environment(int argc, char* argv[], diy::mpi::communicator comm) { ::BoxLib::Initialize(argc, argv, false, comm); DataServices::SetBatchMode(); }
        ~environment()                                                   { ::BoxLib::Finalize(false); }
    };

    class reader
    {
        public:
                  typedef       std::vector<int>                                         Shape;
                  typedef       std::vector<int>                                         Vertex;
                  typedef       std::vector<std::string>                                 VarNameList;
        public:
                  // Construct reader, get domain and variable names in file.
                  reader(std::string            filename,
                         diy::mpi::communicator communicator):
                      communicator_(communicator),
                      dataServices_(filename.c_str(), Amrvis::NEWPLT)
                  {
                      if(!dataServices_.AmrDataOk())
                          throw std::runtime_error("Cannot open BoxLib file");

                      // Get domain extents for all levels
                      const int finest_level = dataServices_.AmrDataRef().FinestLevel();
                      shape_.resize(finest_level + 1);
                      from_.resize(finest_level + 1);
                      to_.resize(finest_level + 1);
                      for (int curr_level = 0; curr_level <= finest_level; ++curr_level)
                      {
                          const Box& domain_box = dataServices_.AmrDataRef().ProbDomain()[curr_level];
                          for (int d = 0; d < BL_SPACEDIM; ++d)
                          {
                              from_[curr_level].push_back(domain_box.smallEnd()[d]);
                              to_[curr_level].push_back(domain_box.bigEnd()[d]);
                              shape_[curr_level].push_back(to_[curr_level][d] - from_[curr_level][d] + 1);
                          }
                      }

                      // Get variable names
                      const Array<string>& plotVarNames = dataServices_.AmrDataRef().PlotVarNames();
                      for (int i = 0; i < plotVarNames.size(); ++i)
                          varnames_.push_back(plotVarNames[i]);
                  }

                  const int          finestLevel()         const                          { return dataServices_.AmrDataRef().FinestLevel(); }
                  const Shape&       shape(int level = 0)  const                          { return shape_[level]; }
                  const Vertex&      from(int level = 0)   const                          { return from_[level]; }
                  const Vertex&      to(int level = 0)     const                          { return to_[level]; }
                  const VarNameList& varnames()            const                          { return varnames_; }

                  void               read(const diy::DiscreteBounds& bounds,    // Region to read
                                          Real* buffer,                         // Buffer where data will be copied to
                                          std::string varname = std::string(),  // Name of variable to use; if empty string, read first variable in file
                                          int finestFillLevel = 0) const        // Finest level to reead/level to flatten to
                  {
                      // If varname is not set, use first variable defined in file
                      if (varname.empty())
                      {
                          // No varname speciied, use first in file
                          varname = varnames_[0];
                          LOG_SEV(info) << "No variable name specified. Using " << varname << " (first variable in file)";
                      }

                      // BoxLib expects a list of boxes to read on all ranks.
                      // Gather the bounds for the read requests on all processors.
                      std::vector< std::vector<char> > buffer_vector;
                      diy::MemoryBuffer                bb;
                      diy::save(bb, bounds);
                      diy::mpi::all_gather(communicator_, bb.buffer, buffer_vector);

                      // Create the list of all boxes to read as well as the distribution mapping
                      // specifying on which rank a box is reead.
                      BoxList    partition_boxes;
                      Array<int> proc_for_box(communicator_.size() + 1); // + 1 is for historic reasons
                      for (int i = 0; i < buffer_vector.size(); ++i)
                      {
                          diy::DiscreteBounds bnds;
                          diy::MemoryBuffer bb; bb.buffer.swap(buffer_vector[i]);
                          diy::load(bb, bnds);
                          IntVect from(bnds.min), to(bnds.max);
                          partition_boxes.push_back(Box(from, to));
                          proc_for_box[i] = i;
                      }
                      DistributionMapping dm(proc_for_box);

                      // Create the multi FAB and read data into it
                      MultiFab mf;
                      BoxArray mf_boxes(partition_boxes);
                      mf.define(mf_boxes, 1, 0, dm, Fab_allocate);
                      dataServices_.AmrDataRef().FillVar(mf, finestFillLevel, varname);

                      // Copy the data from BoxLib to Reeber2
                      const FArrayBox& my_fab = mf[ParallelDescriptor::MyProc()];
                      const Box& my_partition_box = mf_boxes[ParallelDescriptor::MyProc()];
                      typedef GridRef< Real, BL_SPACEDIM > RealGridRef;
                      typedef RealGridRef::Vertex GridVertex;
                      RealGridRef grid(buffer, GridVertex(bounds.max) - GridVertex(bounds.min) + GridVertex::one());
                      for (IntVect iv = my_partition_box.smallEnd(); iv <= my_partition_box.bigEnd(); my_partition_box.next(iv))
                      {
                          GridVertex pos((iv - my_partition_box.smallEnd()).getVect());
                          grid(pos) = my_fab(iv);
                      }
                  }

        private:
                  diy::mpi::communicator     communicator_;
                  mutable DataServices       dataServices_; // Hack, declare "mutable" since BoxLib does not declare const funrctions for reading
                  std::vector < Shape >      shape_;
                  std::vector < Vertex >     from_;
                  std::vector < Vertex >     to_;
                  std::vector< std::string > varnames_;
    };
}
}
}
#endif
