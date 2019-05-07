#ifndef REEBER_AMR_PLOT_READER_H
#define REEBER_AMR_PLOT_READER_H

#include <string>

#include <diy/master.hpp>
#include <diy/assigner.hpp>

void read_amr_plotfile(std::string infile,
                       std::vector<std::string> all_var_names, // all fields that will be read from plotfile
                       int n_mt_vars,                          // sum of first n_mt_vars in all_var_names will be stored in fab of FabBlock,
                                                               // for each variable listed in all_var_names FabBlock will have an extra GridRef
                       diy::mpi::communicator& world,
                       int nblocks,
                       diy::Master& master_reader,
                       diy::MemoryBuffer& header,
                       diy::DiscreteBounds& domain);

#endif //REEBER_AMR_PLOT_READER_H
