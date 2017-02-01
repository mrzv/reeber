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

template<class K, class V, class H, class E, class A>
struct Serialization< map<K,V,H,E,A> >
{
    typedef             map<K,V,H,E,A>   Map;

    static void         save(::diy::BinaryBuffer& bb, const Map& m)
    {
      size_t s = m.size();
      diy::save(bb, s);
      for (auto& x : m)
        diy::save(bb, x);
    }

    static void         load(::diy::BinaryBuffer& bb, Map& m)
    {
      size_t s;
      diy::load(bb, s);
      for (size_t i = 0; i < s; ++i)
      {
        std::pair<K,V> p;
        diy::load(bb, p);
        m.emplace(std::move(p));
      }
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
    typedef     ::reeber::vector<T>             Vector;

    static void save(BinaryBuffer& bb, const Vector& v)     { ::reeber::Serialization<Vector>::save(bb, v); }
    static void load(BinaryBuffer& bb, Vector& v)           { ::reeber::Serialization<Vector>::load(bb, v); }
};

template<class K, class V, class H, class E, class A>
struct Serialization<::reeber::map<K,V,H,E,A>>
{
    typedef     ::reeber::map<K,V,H,E,A>        Map;

    static void save(BinaryBuffer& bb, const Map& m)        { ::reeber::Serialization<Map>::save(bb, m); }
    static void load(BinaryBuffer& bb, Map& m)              { ::reeber::Serialization<Map>::load(bb, m); }
};
#endif

}
