#ifndef REEBER_READER_INTERFACES_H
#define REEBER_READER_INTERFACES_H

#include <boost/algorithm/string/predicate.hpp>

#ifdef REEBER_USE_BOXLIB_READER
#include <algorithm>
#include <reeber/io/boxlib.h>
namespace r = reeber;
#endif

#include "real.h"

struct Reader
{
    typedef                 std::vector<int>                            Shape;
    virtual const Shape&    shape() const                               =0;
    virtual void            read(const diy::DiscreteBounds& bounds,     // Region to read
                                 Real* buffer,                          // Buffer where data will be copied to
                                 bool collective = true) const          =0;
    virtual                 ~Reader()                                   {}
};

struct NumPyReader: public Reader
{
                            NumPyReader(std::string             infn,
                                        diy::mpi::communicator  world):
                                in(world, infn, diy::mpi::io::file::rdonly),
                                numpy_reader(in)                    { numpy_reader.read_header(); }

    virtual const Shape&    shape() const                           { return numpy_reader.shape(); }
    virtual void            read(const diy::DiscreteBounds& bounds,
                                 Real* buffer,
                                 bool collective = true) const      { numpy_reader.read(bounds, buffer, collective); }

    diy::mpi::io::file      in;
    diy::io::NumPy          numpy_reader;
};

#ifdef REEBER_USE_BOXLIB_READER
struct BoxLibReader: public Reader
{
                            BoxLibReader(std::string             infn,
                                        diy::mpi::communicator  world):
                                boxlib_reader(infn, world)          {}


    virtual const Shape&    shape() const                           { return boxlib_reader.shape(); }
    virtual void            read(const diy::DiscreteBounds& bounds,
                                 Real* buffer,
                                 bool collective = true) const      { boxlib_reader.read(bounds, buffer, collective); }

    r::io::BoxLib::Reader   boxlib_reader;
};
#endif

#endif
