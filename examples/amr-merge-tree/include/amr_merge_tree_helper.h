#pragma once

#include <set>

#include "fab-tmt-block.h"
#include "reeber/amr-vertex.h"

using AmrEdgeContainer = reeber::AmrEdgeContainer;
using AmrEdge = reeber::AmrEdge;
using AmrVertexId = r::AmrVertexId;

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
void send_edges_to_neighbors(FabTmtBlock<Real, D>* b, const diy::Master::ProxyWithLink& cp)
{
    bool debug = false;

    if (debug) fmt::print("Called send_edges_to_neighbors for block = {}\n", b->gid);

    auto* l = static_cast<diy::AMRLink*>(cp.link());

    std::set<diy::BlockID> receivers;
    for(int i = 0; i < l->size(); ++i)
    {
        if (l->target(i).gid != b->gid)
        {
            receivers.insert(l->target(i));
        }
    }

    for(const diy::BlockID& receiver : receivers)
    {
        int receiver_gid = receiver.gid;
        if (b->new_receivers_.count(receiver_gid))
        {
            if (debug)
                fmt::print("In send_edges_to_neighbors for block = {}, sending to receiver= {}, cardinality = {}\n",
                           b->gid, receiver_gid, b->get_all_outgoing_edges().size());
            cp.enqueue(receiver, b->get_all_outgoing_edges());
        } else
        {
            if (debug)
                fmt::print("In send_edges_to_neighbors for block = {}, sending to receiver= {} empty container\n",
                           b->gid, receiver_gid);
            cp.enqueue(receiver, r::AmrEdgeContainer());
        }
    }
}

/**
 *
 * @tparam D
 * @param b
 * @param cp
 */
template<unsigned D>
void delete_low_edges(FabTmtBlock<Real, D>* b, const diy::Master::ProxyWithLink& cp)
{
    bool debug = true;

    if (debug) fmt::print("Called delete_low_edges for block = {}\n", b->gid);

    auto* l = static_cast<diy::AMRLink*>(cp.link());

    std::set<diy::BlockID> senders;
    for(int i = 0; i < l->size(); ++i)
    {
        if (l->target(i).gid != b->gid)
        {
            senders.insert(l->target(i));
        }
    }

    for(const diy::BlockID& sender : senders)
    {
        AmrEdgeContainer edges_from_neighbor;
        //if (debug) fmt::print("In delete_low_edges for block = {}, dequeing from sender = {}\n", b->gid, sender.gid);

        cp.dequeue(sender, edges_from_neighbor);

        //if (debug) fmt::print("In delete_low_edges for block = {}, dequed {} edges from sender = {}\n", b->gid, edges_from_neighbor.size(), sender.gid);

        b->delete_low_edges(sender.gid, edges_from_neighbor);

        //if (debug) fmt::print("In delete_low_edges for block = {}, from sender = {}, b->delete_low_edges OK\n", b->gid, sender.gid);
    }

    b->adjust_outgoing_edges();
    if (debug) fmt::print("Exit delete_low_edges for block = {}\n", b->gid);
}

