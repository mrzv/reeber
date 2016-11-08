#ifndef REEBER_POINT_H
#define REEBER_POINT_H

#include <iostream>
#include <vector>
#include <string>
#include <sstream>

#include <diy/point.hpp>

namespace reeber
{

template<class Coordinate, unsigned D>
using Point = ::diy::Point<Coordinate, D>;

} // reeber

namespace opts
{
    template<class T>
    struct Traits;

    template<class C, unsigned D>
    struct Traits< reeber::Point<C,D> >
    {
        static
        std::string     type_string()           { return "POINT"; }
    };
}


#endif // POINT_H
