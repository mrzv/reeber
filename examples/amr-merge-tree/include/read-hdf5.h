#pragma once

#include <string>

#include <diy/master.hpp>
#include <diy/assigner.hpp>
#include <diy/types.hpp>
#include <diy/serialization.hpp>
#include <diy/link.hpp>
#include <reeber/format.h>

void read_from_hdf5_file(std::string infn,
                         std::vector<std::string> all_var_names, // HDF5 only: all fields that will be read from plotfile
                         int n_mt_vars,                          // sum of first n_mt_vars in all_var_names will be stored in fab of FabBlock,
                                                                 // for each variable listed in all_var_names FabBlock will have an extra GridRef
                         diy::mpi::communicator& world,
                         int nblocks,
                         diy::Master& master_reader,
                         diy::ContiguousAssigner& assigner,
                         diy::MemoryBuffer& header,
                         diy::DiscreteBounds& domain);
