#ifndef REEBER_PATH_MERGE_TREE_SERIALIZATION_H
#define REEBER_PATH_MERGE_TREE_SERIALIZATION_H

#include <diy/serialization.hpp>
#include "parallel-tbb.h"
#include "parallel-tbb-serialization.h"
#include "path-merge-tree.h"
#include "range/map.h"

namespace reeber
{

template<class Vertex, class Value>
struct Serialization< PathMergeTree<Vertex, Value> >
{
    typedef     ::reeber::PathMergeTree<Vertex, Value>          PathMergeTree;
    typedef     typename PathMergeTree::Neighbor                Neighbor;

    static void save(::diy::BinaryBuffer& bb, const PathMergeTree& mt)
    {
        diy::save(bb, mt.negate_);
        size_t sz = mt.nodes_.size();
        diy::save(bb, sz);
        for(Neighbor n : mt.nodes_ | range::map_values)
        {
            diy::save(bb, n->vertex);
            diy::save(bb, n->value);
            Neighbor parent = n->parent;
            diy::save(bb, parent->vertex);
        }
    }

    static void load(::diy::BinaryBuffer& bb, PathMergeTree& mt)
    {
        diy::load(bb, mt.negate_);
        size_t sz;
        diy::load(bb, sz);
        for (size_t i = 0; i < sz; ++i)
        {
            Vertex u, v; Value val;
            diy::load(bb, u);
            diy::load(bb, val);
            diy::load(bb, v);

            Neighbor n_u, n_v;
            n_u = mt.add_or_update(u,val);
            n_v = mt.find_or_add(v,0);
            n_u->parent = n_v;
        }
    }
};

}


namespace diy
{

template<class Vertex, class Value>
struct Serialization< ::reeber::PathMergeTree<Vertex, Value> >
{
    typedef     ::reeber::PathMergeTree<Vertex, Value>              PathMergeTree;

    static void save(BinaryBuffer& bb, const PathMergeTree& mt)     { ::reeber::Serialization<PathMergeTree>::save(bb, mt); }
    static void load(BinaryBuffer& bb, PathMergeTree& mt)           { ::reeber::Serialization<PathMergeTree>::load(bb, mt); }
};

}

#endif

