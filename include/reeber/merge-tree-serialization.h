#ifndef REEBER_MERGE_TREE_SERIALIZATION_H
#define REEBER_MERGE_TREE_SERIALIZATION_H

#include <boost/range/adaptor/map.hpp>
#include <diy/serialization.hpp>
#include "merge-tree.h"

namespace reeber
{

template<class Vertex, class Value>
struct Serialization< MergeTree<Vertex, Value> >
{
    typedef     ::reeber::MergeTree<Vertex, Value>          MergeTree;
    typedef     typename MergeTree::Neighbor                Neighbor;

    static void save(::diy::BinaryBuffer& bb, const MergeTree& mt, bool save_vertices = true)
    {
        namespace ba = boost::adaptors;
        diy::save(bb, save_vertices);
        diy::save(bb, mt.negate_);
        diy::save(bb, mt.nodes_.size());
        BOOST_FOREACH(Neighbor n, mt.nodes_ | ba::map_values)
        {
            diy::save(bb, n->vertex);
            diy::save(bb, n->value);
            if (n->parent)
                diy::save(bb, n->parent->vertex);
            else
                diy::save(bb, n->vertex);       // self-loop to indicate the root

            if (save_vertices)
                diy::save(bb, n->vertices);
        }
    }

    static void load(::diy::BinaryBuffer& bb, MergeTree& mt)
    {
        namespace ba = boost::adaptors;
        bool load_vertices;
        diy::load(bb, load_vertices);
        diy::load(bb, mt.negate_);
        size_t s;
        diy::load(bb, s);
        for (size_t i = 0; i < s; ++i)
        {
            Vertex vert, parent; Value val;
            diy::load(bb, vert);
            diy::load(bb, val);
            diy::load(bb, parent);

            Neighbor n = mt.add_or_update(vert, val);

            if (vert != parent)
            {
                n->parent = mt.find_or_add(parent, 0);
                n->parent->children.push_back(n);
            }

            if (load_vertices)
                diy::load(bb, n->vertices);
        }
    }
};

}


namespace diy
{

template<class Vertex, class Value>
struct Serialization< ::reeber::MergeTree<Vertex, Value> >
{
    typedef     ::reeber::MergeTree<Vertex, Value>              MergeTree;

    static void save(BinaryBuffer& bb, const MergeTree& mt)     { ::reeber::Serialization<MergeTree>::save(bb, mt, true); }
    static void load(BinaryBuffer& bb, MergeTree& mt)           { ::reeber::Serialization<MergeTree>::load(bb, mt); }
};

}

#endif
