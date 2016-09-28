#pragma once

#include <diy/serialization.hpp>

namespace reeber
{

#ifdef REEBER_USE_TBB
template<class U>
struct Serialization< vector<U> >
{
    typedef             vector<U>          Vector;

    static void         save(::diy::BinaryBuffer& bb, const Vector& v)
    {
      size_t s = v.size();
      diy::save(bb, s);
      for (auto& x : v)
          diy::save(bb, x);
    }

    static void         load(::diy::BinaryBuffer& bb, Vector& v)
    {
      size_t s;
      diy::load(bb, s);
      v.resize(s);
      for (size_t i = 0; i < s; ++i)
        diy::load(bb, v[i]);
    }
};
#endif

}

namespace diy
{

#ifdef REEBER_USE_TBB
template<class T>
struct Serialization<::reeber::vector<T>>
{
    typedef     ::reeber::vector<T>               Vector;

    static void save(BinaryBuffer& bb, const Vector& v)     { ::reeber::Serialization<Vector>::save(bb, v); }
    static void load(BinaryBuffer& bb, Vector& v)           { ::reeber::Serialization<Vector>::load(bb, v); }
};
#endif

}
