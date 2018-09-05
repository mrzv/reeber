

template<class Real, unsigned D>
void FabTmtBlock<Real, D>::set_mask(const diy::Point<int, D>& v_bounds,
                                    diy::AMRLink* l,
                                    const Real& rho,
                                    bool is_absolute_threshold)
{
    using Position = diy::Point<int, D>;

    int debug_gid = local_.gid();

    bool debug = false;

    bool is_ghost = local_.is_ghost(v_bounds);
    bool is_low = false;
    if (is_absolute_threshold)
        is_low = not is_ghost and (negate_ ? fab_(v_bounds) < rho : fab_(v_bounds) > rho);

    r::AmrVertexId v_idx;
    if (not is_ghost) v_idx = local_.get_vertex_from_global_position(local_.global_position_from_local(v_bounds));

    if (debug) {
        fmt::print("gid = {}, in set_mask, v_bounds = {}, v_idx = {}, value = {}, is_ghost = {}\n", debug_gid, v_bounds, v_idx,
                   fab_(v_bounds), is_ghost);
    }


    // initialization, actual mask to be set later
    if (is_ghost) {
        if (debug) { fmt::print("in set_mask, gid = {}, v_bounds = {}, GHOST detected\n", debug_gid, v_bounds); }
        local_.set_mask(v_bounds, MaskedBox::GHOST);
    } else {
        local_.set_mask(v_bounds, MaskedBox::ACTIVE);
    }

    const int v_ref = local_.refinement();
    const int v_level = local_.level();

    Position v_glob = wrap_point(v_bounds + local_.bounds_from(), domain_, v_ref);

    if (debug) {
        fmt::print("in set_mask, gid = {}, unwrapped v_glob = {}, wrapped = {}, v_idx = {}\n", local_.gid(),
                   v_bounds + local_.bounds_from(), v_glob, v_idx);
    }


    bool mask_set { false };

    // loop over neighbouring blocks one level above
    for (int i = 0; i < l->size(); ++i)
        if (l->level(i) == v_level + 1 && neighbor_contains<D>(i, l, v_glob, v_ref)) {
            if (debug) {
                fmt::print("Is masked ABOVE {} , is_ghost = {}, gid = {}, v_idx = {}\n", local_.gid(), is_ghost,
                           l->target(i).gid, v_idx);
            }
            mask_set = true;
            local_.set_mask(v_bounds, l->target(i).gid);
            break;
        }

    // real cells can only be masked by blocks above, only ghost cells need to look at same level and below
    if (not mask_set and is_ghost) {
        // loop over neighbours on the same level
        for (int i = 0; i < l->size(); ++i) {
            if (l->level(i) == v_level && neighbor_contains<D>(i, l, v_glob, v_ref)) {
                if (debug) {
                    fmt::print("Is masked SAME {} , is_ghost = {}, gid = {}, v_idx = {}\n", l->target(i).gid, is_ghost,
                               local_.gid(), v_idx);
                }
                mask_set = true;
                local_.set_mask(v_bounds, l->target(i).gid);
                break;
            }
        }
    }

    // loop over neighbours one level below. This is not actual masking,
    // but we save this info in ghost cells to get outgoing edges.
    if (not mask_set and is_ghost) {
        for (int i = 0; i < l->size(); ++i) {
            if (l->level(i) == v_level - 1 && neighbor_contains<D>(i, l, v_glob, v_ref)) {
                if (debug) {
                    fmt::print("Is masked SAME {} , is_ghost = {}, gid = {}, v_idx = {}\n", l->target(i).gid, is_ghost,
                               local_.gid(), v_idx);
                }
                mask_set = true;
                local_.set_mask(v_bounds, l->target(i).gid);
                break;
            }
        }
    }

    if (not mask_set and is_low) {
        local_.set_mask(v_bounds, MaskedBox::LOW);
    }

    if (not is_absolute_threshold and not mask_set and not is_ghost)
    {
        // we need to store local sum and local number of unmasked vertices in a block
        // and use this later to mark low vertices
        n_unmasked_++;
        sum_ += fab_(v_bounds);
    }

    if (debug) {
        fmt::print("in set_mask, final mask = {}, gid = {}, v_bounds = {},  v_idx = {}\n",
                   local_.pretty_mask_value(v_bounds), local_.gid(), v_bounds, v_idx);
    }
}

template<class Real, unsigned D>
void FabTmtBlock<Real, D>::set_low(const diy::Point<int, D>& v_bounds,
                                   const Real& absolute_threshold)
{
    if (local_.mask(v_bounds) != MaskedBox::ACTIVE)
        return;
    bool is_low = negate_ ? fab_(v_bounds) < absolute_threshold : fab_(v_bounds) > absolute_threshold;
    if (is_low)
        local_.set_mask(v_bounds, MaskedBox::LOW);
}

template<unsigned D>
r::AmrEdgeContainer
get_vertex_edges(const diy::Point<int, D>& v_glob, const reeber::MaskedBox<D>& local, diy::AMRLink* l,
                 const diy::DiscreteBounds& domain)
{
    using Position = diy::Point<int, D>;

    r::AmrVertexId v_glob_idx = local.get_vertex_from_global_position(v_glob);

    bool debug = false;

    if (debug) { fmt::print("get_vertex_edges called for {}, vertex {}, box = {}\n", v_glob, v_glob_idx, local); }

    r::AmrEdgeContainer result;

    for (const Position& neighb_v_glob : local.outer_edge_link(v_glob)) {

        Position neighb_v_bounds = neighb_v_glob - local.bounds_from();


        if (debug) {
            fmt::print("In get_vertex_edges, v_glob = {}, gid = {}, neighb_v_bounds = {}, neighb_v_glob = {}\n", v_glob,
                       local.gid(), neighb_v_bounds, neighb_v_glob);
        }

        Position wrapped_neighb_vert_glob = wrap_point(neighb_v_glob, domain, local.refinement());

        int gid = local.mask(neighb_v_bounds);

        size_t link_idx = 0;
        bool link_idx_found = false;
        for (; link_idx < (size_t) l->size(); ++link_idx) {
            if (l->target(link_idx).gid == gid) {
                link_idx_found = true;
                break;
            }
        }

        //assert(link_idx_found);
        if (not link_idx_found)
            throw std::runtime_error("gid not found in link");

        auto nb_level = l->level(link_idx);
        Position nb_from = project_point<D>(l->bounds(link_idx).min);
        Position nb_to = project_point<D>(l->bounds(link_idx).max);
        int nb_refinement = l->refinement(link_idx);

        if (debug) {
            fmt::print(
                    "In get_vertex_edges, v_glob = {}, gid = {}, masked by gid= {}, glob_coord = {}, nb_from = {}, nb_to = {}\n",
                    v_glob,
                    local.gid(), gid, wrapped_neighb_vert_glob, nb_from, nb_to);
        }

        //assert(abs(nb_level - local.level()) <= 1);

        if (nb_level <= local.level()) {
            // neighbour is on the same level or below me, neighbouring vertex corresponds
            // to the unique vertex in a neighbouring block
            size_t neighb_vertex_idx = get_vertex_id(wrapped_neighb_vert_glob, local.refinement(), link_idx, l,
                                                     local.mask_grid().c_order());

            result.emplace_back(v_glob_idx, reeber::AmrVertexId { gid, neighb_vertex_idx });

            if (debug) {
                fmt::print("In get_vertex_edges, v_glob = {}, gid = {}, Added edge to idx = {}\n", v_glob, local.gid(),
                           neighb_vertex_idx);
            }

        } else if (nb_level > local.level()) {
            Position masking_box_from, masking_box_to;
            std::tie(masking_box_from, masking_box_to) = refine_vertex(neighb_v_glob, local.refinement(),
                                                                       nb_refinement);
            if (debug)
                fmt::print("In get_vertex_edges, v_glob = {}, gid = {}, masking_box_from = {}, masking_box_to = {}\n",
                           v_glob, local.gid(), masking_box_from, masking_box_to);

            reeber::Box<D> masking_box(masking_box_from, masking_box_to);
            reeber::Box<D> covering_box(masking_box_from - Position::one(), masking_box_to + Position::one());

            for (const Position& masking_position : masking_box.positions()) {
                for (Position covering_position_glob : covering_box.position_link(masking_position)) {
                    Position covering_position_coarsened = wrap_point(
                            coarsen_point(covering_position_glob, nb_refinement,
                                          local.refinement()),
                            domain, local.refinement());
                    if (debug) {
                        fmt::print("In get_vertex_edges, v_glob = {}, gid = {}, cov_vert = {}, cov_vert_loc = {}\n",
                                   v_glob, local.gid(), covering_position_glob,
                                   covering_position_coarsened);
                    }

                    if (covering_position_coarsened == v_glob) {

                        Position masking_position_global = wrap_point(masking_position, domain, nb_refinement);
                        r::AmrVertexId masking_vertex_idx { gid, get_vertex_id(masking_position_global, nb_refinement,
                                                                               link_idx, l,
                                                                               local.mask_grid().c_order()) };

                        if (debug) {
                            fmt::print(
                                    "IN get_vertex_edges, v_glob = {}, gid = {}, masking_position = {}, nb_from = {}, nb_to = {}, adding edge {}\n",
                                    v_glob, local.gid(), masking_position, nb_from, nb_to,
                                    r::AmrEdge { v_glob_idx, masking_vertex_idx });
                        }

                        result.emplace_back(v_glob_idx, masking_vertex_idx);
                    }
                }
            }
        }
    } // loop over outer edge link
    return result;
}


template<class Real, unsigned D>
void FabTmtBlock<Real, D>::compute_outgoing_edges(diy::AMRLink* l, VertexEdgesMap& vertex_to_outgoing_edges)
{
    for (const Vertex& v_glob : local_.active_global_positions()) {
        AmrEdgeContainer out_edges = get_vertex_edges(v_glob, local_, l, domain());
        if (not out_edges.empty()) {
            vertex_to_outgoing_edges[local_.get_vertex_from_global_position(v_glob)] = out_edges;

            std::copy(out_edges.begin(), out_edges.end(), std::back_inserter(initial_edges_));

            for (const AmrEdge& e : out_edges) {
                assert(std::get<0>(e).gid == gid);
                gid_to_outgoing_edges_[std::get<1>(e).gid].push_back(e);
            }

        }
    }
}

template<class Real, unsigned D>
void FabTmtBlock<Real, D>::adjust_original_gids(int sender_gid, FabTmtBlock::GidVector& edges_from_sender)
{
    auto iter = std::find(original_link_gids_.begin(), original_link_gids_.end(), sender_gid);
    bool i_talk_to_sender = iter != original_link_gids_.end();
    bool sender_talks_to_me = std::find(edges_from_sender.begin(), edges_from_sender.end(), gid) != edges_from_sender.end();
    if (i_talk_to_sender and not sender_talks_to_me)
        original_link_gids_.erase(iter);
}

// delete edges from this block that end in a low vertex of a neighbor block.
// edges_from_gid: outogoing eddes we receive from a neigbor block
// subtract them from the edges we keep in gid_to_outgoing_edges_
template<class Real, unsigned D>
void FabTmtBlock<Real, D>::delete_low_edges(int sender_gid, FabTmtBlock::AmrEdgeContainer& edges_from_sender)
{
    bool debug = false;

    auto iter = gid_to_outgoing_edges_.find(sender_gid);
    if (iter == gid_to_outgoing_edges_.end()) {
        // all edges from neighbor must end in low vertex of mine
        if (debug)
            fmt::print("In delete_low_edges in block with gid = {}, sender = {}, no edges from this block, exiting\n",
                       gid, sender_gid);
        return;  // we don't expect any edges from gid
    }

    // edges from neighbor come reversed, correct that
    std::transform(edges_from_sender.begin(), edges_from_sender.end(), edges_from_sender.begin(),
                   &reeber::reverse_amr_edge);

    if (debug) fmt::print("in FabTmtBlock::delete_low_edges, transform OK\n");

    // put edges into set to find the set difference
    std::set<AmrEdge> my_edges { iter->second.begin(), iter->second.end() };
    std::set<AmrEdge> neighbor_edges { edges_from_sender.begin(), edges_from_sender.end() };

#ifdef SEND_COMPONENTS
    for(Component& c : components_)
    {
        c.delete_low_edges(neighbor_edges);
    }
#endif

    int old_n_edges = my_edges.size();

    iter->second.clear();

    if (debug) fmt::print("In delete_low_edges in block with gid = {}, sender = {}, clear OK\n", gid, sender_gid);

    std::set_intersection(my_edges.begin(), my_edges.end(), neighbor_edges.begin(), neighbor_edges.end(),
                          std::back_inserter(iter->second));

    int new_n_edges = iter->second.size();

    if (iter->second.empty())
        gid_to_outgoing_edges_.erase(iter);
    if (debug)
        fmt::print("in Block::delete_low_edges, erase OK, old_n_edges = {}, new_n_edges = {}\n", old_n_edges,
                   new_n_edges);
}

template<class Real, unsigned D>
void FabTmtBlock<Real, D>::adjust_outgoing_edges()
{
    bool debug = false;

    size_t s = initial_edges_.size();
    initial_edges_.clear();
//    initial_edges_.reserve(s);

    for (const auto& gid_edge_vector_pair : gid_to_outgoing_edges_) {
        std::copy(gid_edge_vector_pair.second.begin(), gid_edge_vector_pair.second.end(),
                  std::back_inserter(initial_edges_));
    }

    std::set<int> neighbor_gids;
    for(const AmrEdge& e : initial_edges_)
    {
        neighbor_gids.insert(std::get<1>(e).gid);
    }

    int orig_link_old_size = original_link_gids_.size();

    original_link_gids_ = GidVector(neighbor_gids.begin(), neighbor_gids.end());
    new_receivers_ = neighbor_gids;

    if (debug) fmt::print("In adjust_outgoing_edges for gid = {}, old #edges = {}, new #edges = {}, old link size = {}, new link size = {}, new link = {}\n",
            gid, s, initial_edges_.size(), orig_link_old_size, original_link_gids_.size(), container_to_string(new_receivers_));

    gid_to_outgoing_edges_.clear();
}

template<class Real, unsigned D>
r::AmrVertexId FabTmtBlock<Real, D>::deepest(const AmrVertexId& v) const
{
    auto iter = vertex_to_deepest_.find(v);
    if (vertex_to_deepest_.cend() != iter)
        return iter->second;
    else
        throw std::runtime_error("Deepest not found for vertex");
}

// components_ are only computed once, return their roots
// TODO: cache this in ctor?
template<class Real, unsigned D>
std::vector<r::AmrVertexId> FabTmtBlock<Real, D>::get_original_deepest_vertices() const
{
    std::vector<AmrVertexId> result(components_.size());
    std::transform(components_.begin(), components_.end(), result.begin(), [](const Component& c) { return c.root_; });
    return result;
}


template<class Real, unsigned D>
std::vector<r::AmrVertexId> FabTmtBlock<Real, D>::get_current_deepest_vertices() const
{
    std::set<AmrVertexId> result_set;
    for(const auto& vertex_deepest_pair : vertex_to_deepest_)
    {
        result_set.insert(vertex_deepest_pair.second);
    }
    std::vector<AmrVertexId> result(result_set.begin(), result_set.end());
    return result;
}


template<class Real, unsigned D>
typename FabTmtBlock<Real, D>::Component& FabTmtBlock<Real, D>::find_component(const AmrVertexId& deepest_vertex)
{
    for (Component& cc : components_) {
        if (cc.root_ == deepest_vertex) {
            return cc;
        }
    }
    throw std::runtime_error("Connnected component not found");
}

template<class Real, unsigned D>
void FabTmtBlock<Real, D>::add_component_to_disjoint_sets(const AmrVertexId& deepest_vertex)
{
    //assert(components_disjoint_set_parent_.find(deepest_vertex) == components_disjoint_set_parent_.end());
    //assert(components_disjoint_set_size_.find(deepest_vertex) == components_disjoint_set_size_.end());

    components_disjoint_set_parent_[deepest_vertex] = deepest_vertex;
    components_disjoint_set_size_[deepest_vertex] = 1;
}


template<class Real, unsigned D>
void FabTmtBlock<Real, D>::create_component(const AmrVertexId& deepest_vertex, const AmrEdgeContainer& edges)
{
    bool debug = false;

    if (debug)
        fmt::print("Entered create_component, gid = {}, deepest_vertex = {}, #edges = {}\n", gid, deepest_vertex,
                   edges.size());

    set_deepest(deepest_vertex, deepest_vertex);

    if (debug)
        fmt::print("In create_component, gid = {}, deepest_vertex = {}, before emplace, size = {}\n", gid,
                   deepest_vertex, components_.size());

    components_.emplace_back(deepest_vertex, edges);

    if (debug)
        fmt::print("In create_component, gid = {}, deepest_vertex = {}, added to components, size = {}\n", gid,
                   deepest_vertex, components_.size());

    add_component_to_disjoint_sets(deepest_vertex);

}

template<class Real, unsigned D>
void FabTmtBlock<Real, D>::compute_connected_components(const VertexEdgesMap& vertex_to_outgoing_edges)
{
    bool debug = false;

    const auto& const_tree = mt_;

    if (debug) fmt::print("compute_connected_components called\n");

    VertexNeighborMap component_nodes;

    for (const auto& vert_neighb_pair : const_tree.nodes()) {

        component_nodes.clear();

        if (debug)
            fmt::print("in compute_connected_component, gid = {} processing vertex = {}\n", gid,
                       vert_neighb_pair.first);

        component_nodes[vert_neighb_pair.first] = vert_neighb_pair.second;
        Neighbor u = vert_neighb_pair.second;

        if (!deepest_computed(u->vertex)) {

            if (debug) fmt::print("in compute_connected_compponent, gid = {}, deepest not computed, traversing\n", gid);

            // nodes we traverse while going down,
            // they all should get deepest node
            std::vector<AmrVertexId> visited_neighbors;

            Neighbor u_ = u;
            Neighbor v = std::get<1>(u_->parent());
            Neighbor s = std::get<0>(u_->parent());
            visited_neighbors.push_back(u->vertex);
            visited_neighbors.push_back(v->vertex);
            visited_neighbors.push_back(s->vertex);

            while (u_ != v and not deepest_computed(v->vertex)) {

                visited_neighbors.push_back(u_->vertex);
                visited_neighbors.push_back(v->vertex);
                visited_neighbors.push_back(s->vertex);

                u_ = v;

                v = std::get<1>(u_->parent());
                s = std::get<0>(u_->parent());

                visited_neighbors.push_back(v->vertex);
                visited_neighbors.push_back(s->vertex);

            }

            AmrVertexId deepest_vertex = (deepest_computed(v)) ? deepest(v) : v->vertex;

            AmrEdgeContainer edges;
            for (const AmrVertexId& v : visited_neighbors) {
                auto find_iter = vertex_to_outgoing_edges.find(v);
                if (find_iter != vertex_to_outgoing_edges.end()) {
                    edges.insert(edges.end(), find_iter->second.begin(), find_iter->second.end());
                }
            }

            if (not deepest_computed(v)) {
                // create new component
                if (debug)
                    fmt::print("in compute_connected_compponent, gid = {}, creating new cc, deepest = {}\n", gid,
                               deepest_vertex);
                create_component(deepest_vertex, edges);

            } else {
                if (debug)
                    fmt::print("in compute_connected_compponent, gid = {}, adding edges to existing cc, deepest = {}\n",
                               gid, deepest_vertex);
                // add to existing components new block_ids
#ifdef SEND_COMPONENTS
                find_component(deepest_vertex).add_edges(edges);
#endif
            }

            for (const AmrVertexId& v : visited_neighbors) {
                if (debug)
                    fmt::print(
                            "in compute_connected_compponent, gid = {}, deepest = {}, setting to visited neighbor {}\n",
                            gid, deepest_vertex, v);
                set_deepest(v, deepest_vertex);
            }
        }
    }

    //assert(std::accumulate(const_tree.nodes().cbegin(), const_tree.nodes().cend(), true,
    //                       [this](const bool& prev, const typename VertexNeighborMap::value_type& vn) {
    //                           return prev and this->deepest_computed(vn.second);
    //                       }));
#ifdef SEND_COMPONENTS
    // copy nodes from local merge tree of block to merge trees of components
    for (const auto& vertex_deepest_pair : vertex_to_deepest_) {
        TripletMergeTree& mt = find_component(vertex_deepest_pair.second).merge_tree_;

        AmrVertexId u = vertex_deepest_pair.first;
        Neighbor mt_n_u = const_tree.nodes().at(u);
        Value val = mt_n_u->value;
        //assert(u == mt_n_u->vertex);

        AmrVertexId s = std::get<0>(mt_n_u->parent())->vertex;
        AmrVertexId v = std::get<1>(mt_n_u->parent())->vertex;

        Neighbor n_u, n_s, n_v;

        n_u = mt.add_or_update(u, val);
        n_s = mt.find_or_add(s, 0);
        n_v = mt.find_or_add(v, 0);
        mt.link(n_u, n_s, n_v);
    }
#endif
}


template<class Real, unsigned D>
bool FabTmtBlock<Real, D>::edge_exists(const AmrEdge& e) const
{
    return mt_.contains(std::get<0>(e)) and mt_.contains(std::get<1>(e));
}


template<class Real, unsigned D>
bool FabTmtBlock<Real, D>::edge_goes_out(const AmrEdge& e) const
{
    return mt_.contains(std::get<0>(e)) xor mt_.contains(std::get<1>(e));
}


template<class Real, unsigned D>
bool FabTmtBlock<Real, D>::are_components_connected(const AmrVertexId& deepest_a, const AmrVertexId& deepest_b)
{
    bool result = find_component_in_disjoint_sets(deepest_a) == find_component_in_disjoint_sets(deepest_b);
    // for debug - assume everything is connected
    return result;
}


template<class Real, unsigned D>
bool FabTmtBlock<Real, D>::is_component_connected_to_any_internal(const FabTmtBlock::AmrVertexId& deepest)
{
    for (const auto& cc : components_)
        if (are_components_connected(deepest, cc.root_))
            return true;
    // for debug - assume everything is connected
    //assert(false);
    return false;
}

template<class Real, unsigned D>
r::AmrVertexId FabTmtBlock<Real, D>::find_component_in_disjoint_sets(AmrVertexId x)
{
    bool debug = false;
    if (debug) fmt::print("Entered find_component_in_disjoint_sets, x = {}\n", x);
    while (components_disjoint_set_parent_.at(x) != x) {
        AmrVertexId next = components_disjoint_set_parent_[x];
        if (debug) fmt::print("in find_component_in_disjoint_sets, x = {}, next = {}\n", x, next);
        components_disjoint_set_parent_[x] = components_disjoint_set_parent_.at(next);
        x = next;
    }
    if (debug) fmt::print("Exiting find_component_in_disjoint_sets, x = {}\n", x);
    return x;
}

template<class Real, unsigned D>
void FabTmtBlock<Real, D>::connect_components(const FabTmtBlock::AmrVertexId& deepest_a,
                                              const FabTmtBlock::AmrVertexId& deepest_b)
{
    AmrVertexId a_root = find_component_in_disjoint_sets(deepest_a);
    AmrVertexId b_root = find_component_in_disjoint_sets(deepest_b);

    if (a_root == b_root)
        return;

    auto a_size = components_disjoint_set_size_.at(a_root);
    auto b_size = components_disjoint_set_size_.at(b_root);
    if (a_size < b_size) {
        std::swap(a_root, b_root);
        std::swap(a_size, b_size);
    }

    components_disjoint_set_parent_[b_root] = a_root;
    components_disjoint_set_size_[a_root] += b_size;
}


template<class Real, unsigned D>
int FabTmtBlock<Real, D>::is_done_simple(const std::vector<FabTmtBlock::AmrVertexId>& vertices_to_check)
{
    // check that all edge outgoing from the current region
    // do not start in any component of our original local tree

//        bool debug = (gid == 0);
    bool debug = false;

    if (debug) fmt::print("is_done_simple, gid = {}, #vertices = {}, round = {}\n", gid, vertices_to_check.size(), round_);
    for (const AmrVertexId& v : vertices_to_check) {
        AmrVertexId deepest_v = vertex_to_deepest_.at(v);
        if (is_component_connected_to_any_internal(deepest_v)) {
            if (debug) fmt::print("is_done_simple, gid = {}, v = {}, returning 0\n", gid, v);
            return 0;
        }
    }
    if (debug) fmt::print("is_done_simple, gid = {}, returning 1\n", gid);
    return 1;
}

#ifdef SEND_COMPONENTS
template<class Real, unsigned D>
int FabTmtBlock<Real, D>::are_all_components_done() const
{
    bool debug = false;
    if (debug) fmt::print("are_all_components_done, gid = {}\n", gid);

    for (const Component& c : components_)
        if (not c.is_done())
            return 0;

    if (debug) fmt::print("are_all_components_done, gid = {}, returning 1\n", gid);
    return 1;
}
#endif

template<class Real, unsigned D>
bool FabTmtBlock<Real, D>::gid_must_be_in_link(int gid) const
{
    for (const Component& c : components_)
        if (c.current_neighbors_.count(gid))
            return true;
    return false;
}


template<class Real, unsigned D>
void FabTmtBlock<Real, D>::save(const void* b, diy::BinaryBuffer& bb)
{
    const FabTmtBlock* block = static_cast<const FabTmtBlock*>(b);

    diy::save(bb, block->gid);
//        diy::save(bb, block->local_);
    diy::save(bb, block->mt_);
    diy::save(bb, block->original_tree_);
    //            diy::save(bb, block->components_);
    diy::save(bb, block->domain_);
    diy::save(bb, block->initial_edges_);
//    diy::save(bb, block->gid_to_outgoing_edges_);
    diy::save(bb, block->new_receivers_);
    diy::save(bb, block->processed_receivers_);
    diy::save(bb, block->original_link_gids_);
    diy::save(bb, block->negate_);
    diy::save(bb, block->vertex_to_deepest_);
    diy::save(bb, block->components_disjoint_set_parent_);
    diy::save(bb, block->components_disjoint_set_size_);
    diy::save(bb, block->round_);
}

template<class Real, unsigned D>
void FabTmtBlock<Real, D>::load(void* b, diy::BinaryBuffer& bb)
{
    FabTmtBlock* block = static_cast<FabTmtBlock*>(b);

    diy::load(bb, block->gid);
//        diy::load(bb, block->local_);
    diy::load(bb, block->mt_);
    diy::load(bb, block->original_tree_);
    //            diy::load(bb, block->components_);
    diy::load(bb, block->domain_);
    diy::load(bb, block->initial_edges_);
//    diy::load(bb, block->gid_to_outgoing_edges_);
    diy::load(bb, block->new_receivers_);
    diy::load(bb, block->processed_receivers_);
    diy::load(bb, block->original_link_gids_);
    diy::load(bb, block->negate_);
    diy::load(bb, block->vertex_to_deepest_);
    diy::load(bb, block->components_disjoint_set_parent_);
    diy::load(bb, block->components_disjoint_set_size_);
    diy::load(bb, block->round_);
}

