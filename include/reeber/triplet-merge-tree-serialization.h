#ifndef REEBER_TRIPLET_MERGE_TREE_SERIALIZATION_H
#define REEBER_TRIPLET_MERGE_TREE_SERIALIZATION_H

#include <boost/range/adaptor/map.hpp>
#include <diy/serialization.hpp>
#include "triplet-merge-tree.h"

namespace reeber
{

template<class Vertex, class Value>
struct Serialization< TripletMergeTree<Vertex, Value> >
{
    typedef     ::reeber::TripletMergeTree<Vertex, Value>          TripletMergeTree;
    typedef     typename TripletMergeTree::Neighbor                Neighbor;

    static void save(::diy::BinaryBuffer& bb, const TripletMergeTree& mt)
    {
        namespace ba = boost::adaptors;
        diy::save(bb, mt.negate_);
        diy::save(bb, mt.nodes_.size());
        BOOST_FOREACH(Neighbor n, mt.nodes_ | ba::map_values)
        {
            diy::save(bb, n->vertex);
            diy::save(bb, n->value);
            Neighbor s, v;
            std::tie(s, v) = n->parent;
            diy::save(bb, s->vertex);
            diy::save(bb, v->vertex);
        }
    }

    static void load(::diy::BinaryBuffer& bb, TripletMergeTree& mt)
    {
        namespace ba = boost::adaptors;
        diy::load(bb, mt.negate_);
        size_t s;
        diy::load(bb, s);
        for (size_t i = 0; i < s; ++i)
        {
            Vertex u, s, v; Value val;
            diy::load(bb, u);
            diy::load(bb, val);
            diy::load(bb, s);
            diy::load(bb, v);

            Neighbor n_u, n_s, n_v;
            auto it = mt.nodes().find(u);
            if (it != mt.nodes().end())
            {
                n_u = it->second;
                n_u->value = val;
            }
            else n_u = mt.add(u, val);

            it = mt.nodes().find(s);
            if (it != mt.nodes().end()) n_s = it->second;
            else n_s = mt.add(s, 0);

            it = mt.nodes().find(v);
            if (it != mt.nodes().end()) n_v = it->second;
            else n_v = mt.add(v, 0);

            mt.link(n_u, n_s, n_v);
        }
    }
};

}


namespace diy
{

template<class Vertex, class Value>
struct Serialization< ::reeber::TripletMergeTree<Vertex, Value> >
{
    typedef     ::reeber::TripletMergeTree<Vertex, Value>              TripletMergeTree;

    static void save(BinaryBuffer& bb, const TripletMergeTree& mt)     { ::reeber::Serialization<TripletMergeTree>::save(bb, mt); }
    static void load(BinaryBuffer& bb, TripletMergeTree& mt)           { ::reeber::Serialization<TripletMergeTree>::load(bb, mt); }
};

}

#endif
