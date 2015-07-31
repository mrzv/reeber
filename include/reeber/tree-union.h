#ifndef REEBER_TREE_TOPOLOGY_H
#define REEBER_TREE_TOPOLOGY_H

#include <stdexcept>
#include <boost/foreach.hpp>
#include <boost/range/adaptor/map.hpp>

namespace reeber
{

namespace ba = boost::adaptors;

template<class MergeTree_, class Edges_>
class TreeUnionTopology
{
    public:
        typedef             MergeTree_                      MergeTree;
        typedef             Edges_                          Edges;
        typedef             typename MergeTree::Vertex      Vertex;
        typedef             typename MergeTree::Neighbor    Neighbor;
        typedef             std::set<Vertex>                Link;
        typedef             std::set<Vertex>                Vertices;
        typedef             std::vector<MergeTree>          Trees;

        class VertexUnionIterator;
        typedef             boost::iterator_range<VertexUnionIterator>      VertexUnionRange;

    public:
                            TreeUnionTopology(const std::vector<MergeTree>& trees, const Edges& edges):
                                trees_(trees), edges_(edges)    {}

        VertexUnionRange    vertices() const;
        Link                link(const Vertex& v) const;

        // nonesense; only used for reserve
        size_t              size() const                    { return 0; }

    private:
        const Trees&                        trees_;
        const Edges&                        edges_;
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

template<class MT, class E>
typename reeber::TreeUnionTopology<MT,E>::VertexUnionRange
reeber::TreeUnionTopology<MT,E>::
vertices() const
{
    return VertexUnionRange(VertexUnionIterator(trees_), VertexUnionIterator(trees_,true));
}

template<class MT, class E>
typename reeber::TreeUnionTopology<MT,E>::Link
reeber::TreeUnionTopology<MT,E>::
link(const Vertex& v) const
{
    Link link;
    BOOST_FOREACH(const MergeTree& mt, trees_)
    {
        if (mt.contains(v))
        {
            Neighbor n = mt[v];
            //if (n->parent)
            //    link.insert(n->parent->vertex);
            BOOST_FOREACH(Neighbor c, n->children)
                link.insert(c->vertex);
        }
    }
    BOOST_FOREACH(const Vertex& u, edges_(v))
        link.insert(u);
    return link;
}

template<class MT, class E>
class reeber::TreeUnionTopology<MT,E>::VertexUnionIterator:
    public boost::iterator_facade<VertexUnionIterator,
                                  Vertex,
                                  boost::forward_traversal_tag,
                                  Vertex,
                                  std::ptrdiff_t>
{
    typedef     boost::iterator_facade<VertexUnionIterator,
                                       Vertex,
                                       boost::forward_traversal_tag,
                                       Vertex,
                                       std::ptrdiff_t>              Parent;


    public:
        typedef     typename Parent::value_type                     value_type;
        typedef     typename Parent::difference_type                difference_type;
        typedef     typename Parent::reference                      reference;

        typedef     typename MergeTree::VertexNeighborMap::const_iterator       NodeIterator;

                    VertexUnionIterator(const Trees& trees, bool end = false):
                        trees_(trees)
        {
            for (unsigned i = 0; i < trees_.size(); ++i)
            {
                if (end)
                    iters_.push_back(trees[i].nodes().end());
                else
                    iters_.push_back(trees[i].nodes().begin());
            }
        }

    private:
        void        increment()
        {
            Vertex v = dereference();
            for (unsigned i = 0; i < iters_.size(); ++i)
                if (iters_[i] != trees_[i].nodes().end() && iters_[i]->first == v)
                    ++iters_[i];
        }

        bool        equal(const VertexUnionIterator& other) const
        {
            for (unsigned i = 0; i < iters_.size(); ++i)
                if (iters_[i] != other.iters_[i])
                    return false;
            return true;
        }

        reference   dereference() const
        {
            NodeIterator min;
            unsigned i = 0;
            for (; i < iters_.size(); ++i)
                if (iters_[i] != trees_[i].nodes().end())
                {
                    min = iters_[i];
                    break;
                }
            for (; i < iters_.size(); ++i)
                if (iters_[i] != trees_[i].nodes().end() && iters_[i]->first < min->first)
                    min = iters_[i];

            return min->first;
        }

        friend class ::boost::iterator_core_access;

    private:
        const Trees&                        trees_;
        std::vector<NodeIterator>           iters_;
};

#endif
