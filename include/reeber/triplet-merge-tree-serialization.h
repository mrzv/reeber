#ifndef REEBER_TRIPLET_MERGE_TREE_SERIALIZATION_H
#define REEBER_TRIPLET_MERGE_TREE_SERIALIZATION_H

#include <boost/range/adaptor/map.hpp>
#include <boost/foreach.hpp>
#include <diy/serialization.hpp>
#include "parallel-tbb.h"
#include "triplet-merge-tree.h"

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

template<class Vertex, class Value>
struct Serialization< TripletMergeTree<Vertex, Value> >
{
    typedef     ::reeber::TripletMergeTree<Vertex, Value>          TripletMergeTree;
    typedef     typename TripletMergeTree::Neighbor                Neighbor;

    static void save(::diy::BinaryBuffer& bb, const TripletMergeTree& mt, bool save_vertices = true)
    {
        namespace ba = boost::adaptors;
        diy::save(bb, save_vertices);
        diy::save(bb, mt.negate_);
        size_t sz = mt.nodes_.size();
        diy::save(bb, sz);
        BOOST_FOREACH(Neighbor n, mt.nodes_ | ba::map_values)
        {
            diy::save(bb, n->vertex);
            diy::save(bb, n->value);
            Neighbor s, v;
            std::tie(s, v) = n->parent();
            diy::save(bb, s->vertex);
            diy::save(bb, v->vertex);
            if (save_vertices)
                diy::save(bb, n->vertices);
        }
    }

    static void load(::diy::BinaryBuffer& bb, TripletMergeTree& mt)
    {
        namespace ba = boost::adaptors;
        bool load_vertices;
        diy::load(bb, load_vertices);
        diy::load(bb, mt.negate_);
        size_t sz;
        diy::load(bb, sz);
        for (size_t i = 0; i < sz; ++i)
        {
            Vertex u, s, v; Value val;
            diy::load(bb, u);
            diy::load(bb, val);
            diy::load(bb, s);
            diy::load(bb, v);

            Neighbor n_u, n_s, n_v;
            n_u = mt.add_or_update(u,val);
            n_s = mt.find_or_add(s,0);
            n_v = mt.find_or_add(v,0);

            mt.link(n_u, n_s, n_v);

            if (load_vertices)
                diy::load(bb, n_u->vertices);
        }
    }
};

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

template<class Vertex, class Value>
struct Serialization< ::reeber::TripletMergeTree<Vertex, Value> >
{
    typedef     ::reeber::TripletMergeTree<Vertex, Value>              TripletMergeTree;

    static void save(BinaryBuffer& bb, const TripletMergeTree& mt)     { ::reeber::Serialization<TripletMergeTree>::save(bb, mt, true); }
    static void load(BinaryBuffer& bb, TripletMergeTree& mt)           { ::reeber::Serialization<TripletMergeTree>::load(bb, mt); }
};

}

#endif
