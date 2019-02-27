#ifndef REEBER_AMR_PLOT_READER_H
#define REEBER_AMR_PLOT_READER_H

#include <string>

#include <diy/master.hpp>
#include <diy/assigner.hpp>

void read_amr_plotfile(std::string infile,
                       std::string varName,
                       diy::mpi::communicator& world,
                       int nblocks,
                       diy::Master& master_reader,
//                       diy::ContiguousAssigner& assigner,
                       diy::MemoryBuffer& header,
                       diy::DiscreteBounds& domain);

#endif //REEBER_AMR_PLOT_READER_H
