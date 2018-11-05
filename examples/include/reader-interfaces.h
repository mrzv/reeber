#ifndef REEBER_READER_INTERFACES_H
#define REEBER_READER_INTERFACES_H

#include <stdexcept>
#include <boost/algorithm/string/predicate.hpp>

#include <diy/io/numpy.hpp>

#include <reeber/box.h>

#ifdef REEBER_USE_BOXLIB_READER
#include <algorithm>
#include <reeber/io/boxlib.h>
#endif

namespace r = reeber;

#include "reeber-real.h"

struct Reader
{
    typedef                 std::vector<int>                            Shape;
    typedef                 std::vector<Real>                           Size;
    typedef                 r::OffsetGrid<Real, 3>                      OffsetGrid;
    virtual const Shape&    shape() const                               =0;
    virtual const Size&     cell_size() const                           =0;
    virtual void            read(const diy::DiscreteBounds& bounds,     // Region to read
                                 Real* buffer,                          // Buffer where data will be copied to
                                 bool collective = true) const          =0;
    virtual OffsetGrid*     read(const r::Box<3> &core);
    virtual                 ~Reader()                                   {}
    static Reader*          create(std::string, diy::mpi::communicator);
};

inline Reader::OffsetGrid* Reader::read(const r::Box<3> &bounds)
{
    diy::DiscreteBounds read_bounds{3};
    for (int d=0; d<3; ++d)
    {
        read_bounds.min[d] = bounds.from()[d];
        read_bounds.max[d] = bounds.to()[d];
    }
    OffsetGrid *og = new OffsetGrid(shape(), &read_bounds.min[0], &read_bounds.max[0]);
    read(read_bounds, og->data(), true);
    return og;
}

struct NumPyReader: public Reader
{
                            NumPyReader(std::string             infn,
                                        diy::mpi::communicator  world):
                                in(world, infn, diy::mpi::io::file::rdonly),
                                numpy_reader(in), dx(3, 1.0)
    {
        unsigned word_size = numpy_reader.read_header();
        if (word_size != sizeof(Real))
            throw std::runtime_error("Data type does not match");
    }

    virtual const Shape&    shape() const                           { return numpy_reader.shape(); }
    virtual const Size&     cell_size() const                       { return dx; }

    virtual void            read(const diy::DiscreteBounds& bounds,
                                 Real* buffer,
                                 bool collective = true) const      { numpy_reader.read(bounds, buffer, collective); }

    diy::mpi::io::file      in;
    diy::io::NumPy          numpy_reader;
    Size                    dx;
};

#ifdef REEBER_USE_BOXLIB_READER
struct BoxLibReader: public Reader
{
                            BoxLibReader(std::string            infn,
                                         diy::mpi::communicator world):
                                boxlib_reader(infn, world)          {}


    virtual const Shape&    shape() const                           { return boxlib_reader.shape(); }
    virtual const Size&     cell_size() const                       { return boxlib_reader.cell_size(); }
    virtual void            read(const diy::DiscreteBounds& bounds,
                                 Real* buffer,
                                 bool collective = true) const      { boxlib_reader.read(bounds, buffer, collective); }

    r::io::BoxLib::Reader   boxlib_reader;
};

struct BoxLibInSituCopier: public Reader
{
                               BoxLibInSituCopier(const MultiFab&        simulation_data,
                                                  const Geometry&        geometry,
                                                  int                    component,
                                                  diy::mpi::communicator world):
                                   boxlib_copier(simulation_data, geometry,
                                                 component, world)     {}


    virtual const Shape&       shape() const                           { return boxlib_copier.shape(); }
    virtual const Size&        cell_size() const                       { return boxlib_copier.cell_size(); }
    virtual void               read(const diy::DiscreteBounds& bounds,
                                     Real* buffer,
                                     bool collective = true) const     { boxlib_copier.read(bounds, buffer, collective); }

    r::io::BoxLib::InSituCopier boxlib_copier;
};

#endif

inline Reader* Reader::create(std::string infn, diy::mpi::communicator world)
{
    Reader* reader_ptr;
#ifdef REEBER_USE_BOXLIB_READER
    if (boost::algorithm::ends_with(infn, ".npy"))
        reader_ptr = new NumPyReader(infn, world);
    else
        reader_ptr = new BoxLibReader(infn, world);
#else
    reader_ptr      = new NumPyReader(infn, world);
#endif
    return reader_ptr;
}

#endif
