#ifndef REEBER_TREE_TOPOLOGY_H
#define REEBER_TREE_TOPOLOGY_H

#include <stdexcept>
#include <boost/foreach.hpp>
#include <boost/range/adaptor/map.hpp>

namespace reeber
{

namespace ba = boost::adaptors;

template<class MergeTree_>
class TreeUnionTopology
{
    public:
        typedef             MergeTree_                      MergeTree;
        typedef             typename MergeTree::Vertex      Vertex;
        typedef             typename MergeTree::Neighbor    Neighbor;
        typedef             std::set<Vertex>                Link;
        typedef             std::set<Vertex>                Vertices;

    public:
                            TreeUnionTopology(const std::vector<MergeTree>& trees):
                                trees_(trees)               {}

        Vertices            vertices() const;
        Link                link(const Vertex& v) const;

        // nonesense; only used for reserve
        size_t              size() const                    { return 0; }

    private:
        const std::vector<MergeTree>&       trees_;
};

template<class MergeTree_>
class TreeUnionFunction
{
    public:
        typedef             MergeTree_                      MergeTree;
        typedef             typename MergeTree::Vertex      Vertex;
        typedef             typename MergeTree::Value       Value;

    public:
                            TreeUnionFunction(const std::vector<MergeTree>& trees):
                                trees_(trees)               {}

        Value               operator()(const Vertex& v) const
        {
            BOOST_FOREACH(const MergeTree& mt, trees_)
                if (mt.contains(v))
                    return mt[v]->value;
            throw std::domain_error("Value not found in TreeUnionFunction");
        }

    private:
        const std::vector<MergeTree>&       trees_;
};

}

// TODO: the current implementation of vertices() is rather inefficient without return value optimization
template<class MT>
typename reeber::TreeUnionTopology<MT>::Vertices
reeber::TreeUnionTopology<MT>::
vertices() const
{
    Vertices vertices;
    BOOST_FOREACH(const MergeTree& mt, trees_)
        BOOST_FOREACH(const Vertex& v, mt.nodes() | ba::map_keys)
            vertices.insert(v);

    return vertices;
}

template<class MT>
typename reeber::TreeUnionTopology<MT>::Link
reeber::TreeUnionTopology<MT>::
link(const Vertex& v) const
{
    Link link;
    BOOST_FOREACH(const MergeTree& mt, trees_)
    {
        if (mt.contains(v))
        {
            Neighbor n = mt[v];
            if (n->parent)
                link.insert(n->parent->vertex);
            BOOST_FOREACH(Neighbor c, n->children)
                link.insert(c->vertex);
        }
    }
    return link;
}

#endif
