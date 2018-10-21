#pragma once

#include <set>

#include "reeber-real.h"
#include "fab-tmt-block.h"
#include "reeber/amr-vertex.h"

using AmrEdgeContainer = reeber::AmrEdgeContainer;
using AmrEdge = reeber::AmrEdge;
using AmrVertexId = r::AmrVertexId;

std::set<diy::BlockID> link_unique(diy::AMRLink *amr_link, int gid)
{
    std::set<diy::BlockID> result;
    for (int i = 0; i < amr_link->size(); ++i)
    {
        if (amr_link->target(i).gid != gid)
        {
            result.insert(amr_link->target(i));
        }
    }
    return result;
}

/**
 *
 * send outgoing edges computed in b to all its neighbors
 * used to symmetrize the edge set in the beginning of the algorithm
 *
 * @param b FabTmtBlock
 * block that will send its edges
 *
 * @param cp diy::Master::ProxyWithLink
 * communication proxy
 *
 */

template<unsigned D>
void send_edges_to_neighbors(FabTmtBlock<Real, D> *b, const diy::Master::ProxyWithLink& cp)
{
    bool debug = false;

    if (debug) fmt::print("Called send_edges_to_neighbors for block = {}\n", b->gid);

    auto *l = static_cast<diy::AMRLink *>(cp.link());

    for (const diy::BlockID& receiver : link_unique(l, b->gid))
    {
        int receiver_gid = receiver.gid;
        cp.enqueue(receiver, b->gid_to_outgoing_edges_[receiver_gid]);
//        debug = receiver_gid == 0 and (b->gid == 3 or b->gid == 12 || b->gid == 15);
        if (debug)
            fmt::print("In send_edges_to_neighbors for block = {}, sending to receiver= {}, cardinality = {}\n",
                       b->gid, receiver_gid, b->gid_to_outgoing_edges_[receiver_gid].size());

    }
}

/**
 *
 * @tparam D
 * @param b
 * @param cp
 */
template<unsigned D>
void delete_low_edges(FabTmtBlock<Real, D> *b, const diy::Master::ProxyWithLink& cp)
{
    auto *l = static_cast<diy::AMRLink *>(cp.link());

    for (const diy::BlockID& sender : link_unique(l, b->gid))
    {
        AmrEdgeContainer edges_from_neighbor;
        cp.dequeue(sender, edges_from_neighbor);
        b->delete_low_edges(sender.gid, edges_from_neighbor);
    }

    b->adjust_outgoing_edges();
    b->sparsify_prune_original_tree();
}

template<class Link>
bool link_contains_gid(Link *link, int gid)
{
    for (int i = 0; i < link->size(); ++i)
        if (link->target(i).gid == gid)
            return true;
    return false;
}


