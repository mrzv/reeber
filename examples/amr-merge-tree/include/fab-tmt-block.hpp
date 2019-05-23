template<class Real, unsigned D>
bool FabTmtBlock<Real, D>::cmp(Real a, Real b) const
{
    if (negate_)
        return a > b;
    else
        return a < b;
}

template<class Real, unsigned D>
void FabTmtBlock<Real, D>::set_mask(const diy::Point<int, D>& v_mask,
        diy::AMRLink* l,
        const Real& rho,
        bool is_absolute_threshold)
{
    using Position = diy::Point<int, D>;

    int debug_gid = local_.gid();

    bool debug = false; //gid == 24316;
//    debug = (gid == 63) and (v_mask[0] == 2 and v_mask[1] == 2 and v_mask[2] == 65);

    auto v_local = local_.local_position_from_mask(v_mask);
    auto v_global = local_.global_position_from_local(v_local);

    bool is_in_core = local_.core_contains_global(v_global);

    bool is_on_boundary = local_.is_on_boundary(v_global);

    bool is_ghost = local_.is_outer(v_mask);


    // TODO: get rid of this!
//    if (is_in_core and is_on_boundary)
//    {
//        Real new_value = 0;
//        int n = 0;
//        for(auto internal_v_global : local_.inside_link(v_global))
//        {
//            auto internal_v = local_.local_position_from_global(internal_v_global);
////            if (debug) fmt::print("HERE: internal_v = {}, internal_v_global = {}, v_global = {}\n", internal_v, internal_v_global, v_global);
//            new_value += fab_(internal_v);
//            ++n;
//        }
//
//        if (n  == 0)
//            throw std::runtime_error("zero division in correction code");
//        new_value /= n;
//
//        if (debug)
//        {
//            Real old_value = fab_(v_local);
//            fab_(v_local) = new_value;
//            if (debug)
//            {
//                fmt::print(
//                        "gid = {}, in set_mask, v_mask = {}, old_value = {}, new_value = {}, is_ghost = {}, is_on_boundary = {}, is_in_core = {}\n",
//                        debug_gid, v_mask,
//                        old_value, new_value, is_ghost, is_on_boundary, is_in_core);
//            }
//        }
//    }

    Real value = std::numeric_limits<Real>::infinity();
    if (is_in_core)
    {
        value = fab_(v_local);
        //fmt::print("VALUE = {}, v_local = {}, gid = {}\n", value, v_local, gid);
    }

    bool is_low = false;
    if (is_absolute_threshold and is_in_core)
        is_low = not is_ghost and
                cmp(rho, value);   //(negate_ ? fab_(v_mask) < rho : fab_(v_mask) > rho);

    r::AmrVertexId v_idx;
    // does not matter here
    if (not is_ghost) v_idx = local_.get_vertex_from_global_position(local_.global_position_from_local(v_mask));

    if (debug)
    {
        if (n_debug_printed_core_ < 10 and not is_on_boundary and is_in_core)
        {
            n_debug_printed_core_++;
            fmt::print(
                    "PRINTING CORE gid = {}, in set_mask, v_mask = {}, v_idx = {}, value = {}, is_ghost = {}, is_on_boundary = {}, is_in_core = {}\n",
                    debug_gid, v_mask,
                    v_idx,
                    value, is_ghost, is_on_boundary, is_in_core);
        }
        if (n_debug_printed_bdry_ < 10 and is_on_boundary)
        {
            n_debug_printed_bdry_++;
            fmt::print(
                    "PRINTING BDRY gid = {}, in set_mask, v_mask = {}, v_idx = {}, value = {}, is_ghost = {}, is_on_boundary = {}, is_in_core = {}\n",
                    debug_gid, v_mask,
                    v_idx,
                    value, is_ghost, is_on_boundary, is_in_core);
        }
    }

    // initialization, actual mask to be set later
    if (is_ghost)
    {
        if (debug)
        { fmt::print("in set_mask, gid = {}, v_mask = {}, GHOST detected\n", debug_gid, v_mask); }
        local_.set_mask(v_mask, MaskedBox::GHOST);
    } else
    {
        local_.set_mask(v_mask, MaskedBox::ACTIVE);
    }

    const int v_ref = local_.refinement();
    const int v_level = local_.level();

    Position v_glob = wrap_point(v_mask + local_.mask_from(), domain_, v_ref, true);

//    if (debug)
//    {
//        fmt::print("in set_mask, gid = {}, unwrapped v_glob = {}, wrapped = {}, v_idx = {}, domain = [{} - {}], local = {}\n", local_.gid(),
//                   v_mask + local_.mask_from(), v_glob, v_idx, domain_.min, domain_.max, local_);
//    }


    bool mask_set{false};

    // loop over neighbouring blocks one level above
    for(int i = 0; i < l->size(); ++i)
        if (l->level(i) == v_level + 1 && neighbor_contains<D>(i, l, v_glob, v_ref))
        {
            if (debug)
            {
                fmt::print("Is masked ABOVE {} , is_ghost = {}, gid = {}, v_idx = {}\n", local_.gid(), is_ghost,
                        l->target(i).gid, v_idx);
            }
            mask_set = true;
            local_.set_mask(v_mask, l->target(i).gid);
            if (not is_ghost)
                n_masked_++;
            break;
        }

    // real cells can only be masked by blocks above, only ghost cells need to look at same level and below
    if (not mask_set and is_ghost)
    {
        // loop over neighbours on the same level
        for(int i = 0; i < l->size(); ++i)
        {
            if (l->level(i) == v_level && neighbor_contains<D>(i, l, v_glob, v_ref))
            {
                if (debug)
                {
                    fmt::print("Is masked SAME {} , is_ghost = {}, gid = {}, v_idx = {}\n", l->target(i).gid, is_ghost,
                            local_.gid(), v_idx);
                }
                mask_set = true;
                local_.set_mask(v_mask, l->target(i).gid);
                if (not is_ghost)
                    n_masked_++;
                break;
            }
        }
    }

    // loop over neighbours one level below. This is not actual masking,
    // but we save this info in ghost cells to get outgoing edges.
    if (not mask_set and is_ghost)
    {
        for(int i = 0; i < l->size(); ++i)
        {
            if (l->level(i) == v_level - 1 && neighbor_contains<D>(i, l, v_glob, v_ref))
            {
                if (debug)
                {
                    fmt::print("Is masked FAKE {} , is_ghost = {}, gid = {}, v_idx = {}\n", l->target(i).gid, is_ghost,
                            local_.gid(), v_idx);
                }
                mask_set = true;
                local_.set_mask(v_mask, l->target(i).gid);
                break;
            }
        }
    }

    if (not mask_set and is_low)
    {
        local_.set_mask(v_mask, MaskedBox::LOW);
    }

    if (is_ghost and is_in_core)
        throw std::runtime_error("Contradiction");

    if (not is_absolute_threshold and not mask_set and not is_ghost)
    {
        // we need to store local sum and local number of unmasked vertices in a block
        // and use this later to mark low vertices
        n_unmasked_++;
        auto v_local = local_.local_position_from_mask(v_mask);
        sum_ += value;
        if (debug)
            fmt::print(
                    "FAB INFO gid = {}, sum = {}, v_mask = {}, v_local = {}, f(v) = {}, index(v) = {} ,stride = {}, shape = {}, data = {}, size = {}\n",
                    gid, sum_, v_mask, v_local, fab_(v_local), fab_.index(v_local), fab_.stride_, fab_.shape(),
                    (void*) fab_.data(), sizeof(fab_.data()[fab_.index(v_local)]));
    }

    if (debug)
    {
        fmt::print("in set_mask, is_ghost = {}, final mask = {}, {}, gid = {}, v_mask = {},  v_idx = {}\n",
                is_ghost, local_.pretty_mask_value(v_mask), local_.mask(v_mask), local_.gid(), v_mask, v_idx);
    }

    if (is_ghost and local_.mask(v_mask) == gid)
    {
        fmt::print("in set_mask, is_ghost = {}, final mask = {}, {}, gid = {}, v_mask = {},  v_idx = {}\n",
                is_ghost, local_.pretty_mask_value(v_mask), local_.mask(v_mask), local_.gid(), v_mask, v_idx);
        fmt::print("Error, same gid in mask\n");
        throw std::runtime_error("Bad mask");
    }
}

template<class Real, unsigned D>
void FabTmtBlock<Real, D>::set_low(const diy::Point<int, D>& v_bounds,
        const Real& absolute_threshold)
{
    auto v_mask = local_.mask_position_from_local(v_bounds);
    if (local_.mask(v_mask) != MaskedBox::ACTIVE)
        return;
    bool is_low = cmp(absolute_threshold,
            fab_(v_bounds)); //   negate_ ? fab_(v_bounds) < absolute_threshold : fab_(v_bounds) > absolute_threshold;
    if (is_low)
    {
        n_low_++;
        local_.set_mask(v_mask, MaskedBox::LOW);
    } else
    {
        n_active_++;
//        if (gid == 0) fmt::print("HERE ACTIVE: {}\n", local_.global_position_from_local(v_bounds));
    }
}

template<class Real, unsigned D>
void FabTmtBlock<Real, D>::init(Real absolute_rho, diy::AMRLink* amr_link)
{
    bool debug = false; //gid == 1 or gid == 100;
    std::string debug_prefix = "In FabTmtBlock::init, gid = " + std::to_string(gid);

    diy::for_each(local_.bounds_shape(), [this, absolute_rho](const Vertex& v_bounds) {
        this->set_low(v_bounds, absolute_rho);
    });

    if (debug)
        fmt::print("{}, absolute_rho = {}, n_active = {}, n_masked = {}, total_size = {}, level = {}\n", debug_prefix,
                absolute_rho, n_active_,
                n_masked_, local_.core_shape()[0] * local_.core_shape()[1] * local_.core_shape()[2], level());

    reeber::compute_merge_tree2(current_merge_tree_, local_, fab_);
    current_merge_tree_.make_deep_copy(original_tree_);

    VertexEdgesMap vertex_to_outgoing_edges;
    compute_outgoing_edges(amr_link, vertex_to_outgoing_edges);

    compute_original_connected_components(vertex_to_outgoing_edges);
    sparsify_prune_original_tree();

    // TODO: delete this? we are going to overwrite this in adjust_outgoing_edges anyway
    for(int i = 0; i < amr_link->size(); ++i)
    {
        if (amr_link->target(i).gid != gid)
        {
            new_receivers_.insert(amr_link->target(i).gid);
            original_link_gids_.push_back(amr_link->target(i).gid);
        }
    }

    if (debug)
        fmt::print("{}, constructed, refinement = {}, level = {}, local = {},  #components = {}\n",
                debug_prefix, refinement(), level(), local_, components_.size());
    if (debug)
    {
        int n_edges = 0;
        for(auto& gid_edges : gid_to_outgoing_edges_)
        { n_edges += gid_edges.second.size(); }
        fmt::print("{},  constructed, tree.size = {}, new_receivers.size = {}, n_edges = {}\n",
                debug_prefix, current_merge_tree_.size(), new_receivers_.size(), n_edges);
    }

    assert(current_merge_tree_.size() >= original_tree_.size());
}

template<class Real, unsigned D>
void FabTmtBlock<Real, D>::sparsify_prune_original_tree()
{
    std::unordered_set<AmrVertexId> special;
    for(const AmrEdge& out_edge : get_all_outgoing_edges())
    {
        special.insert(std::get<0>(out_edge));
    }
//    r::remove_degree_two(original_tree_, [&special](AmrVertexId u) { return special.find(u) != special.end(); });
    r::sparsify(original_tree_, [&special](AmrVertexId u) { return special.find(u) != special.end(); });
}

template<unsigned D>
r::AmrEdgeContainer
get_vertex_edges(const diy::Point<int, D>& v_glob, const reeber::MaskedBox<D>& local, diy::AMRLink* l,
        const diy::DiscreteBounds& domain, bool debug = false)
{
    using Position = diy::Point<int, D>;

    r::AmrVertexId v_glob_idx = local.get_vertex_from_global_position(v_glob);

//    bool debug = v_glob_idx == r::AmrVertexId(0, 64) or v_glob_idx == r::AmrVertexId(1, 4486);

    if (debug)
    { fmt::print("get_vertex_edges called for {}, vertex {}, box = {}\n", v_glob, v_glob_idx, local); }

    r::AmrEdgeContainer result;

    for(const Position& neighb_v_glob : local.outer_edge_link(v_glob))
    {

        Position neighb_v_bounds = neighb_v_glob - local.bounds_from();

        Position wrapped_neighb_vert_glob = wrap_point(neighb_v_glob, domain, local.refinement());

        if (debug)
        {
            fmt::print(
                    "In get_vertex_edges, v_glob = {}, gid = {}, neighb_v_bounds = {}, neighb_v_glob = {}, wrapped_neighb_vert_glob = {}\n",
                    v_glob,
                    local.gid(), neighb_v_bounds, neighb_v_glob, wrapped_neighb_vert_glob);
        }

        Position neighb_v_mask = local.mask_position_from_local(neighb_v_bounds);
        int masking_gid = local.mask(neighb_v_mask);

        size_t link_idx = 0;
        bool link_idx_found = false;
        for(; link_idx < (size_t) l->size(); ++link_idx)
        {
            if (l->target(link_idx).gid == masking_gid)
            {
                link_idx_found = true;
                break;
            }
        }

        //assert(link_idx_found);
        if (not link_idx_found)
        {
            fmt::print("Error here: masking_gid = {}\n", local.pretty_mask_value(masking_gid));
            throw std::runtime_error("masking_gid not found in link");
        }

        auto nb_level = l->level(link_idx);
        Position nb_from = point_from_dynamic_point<D>(l->bounds(link_idx).min);
        Position nb_to = point_from_dynamic_point<D>(l->bounds(link_idx).max);
        // TODO: vector refinement
        int nb_refinement = l->refinement(link_idx)[0];

        if (debug)
        {
            fmt::print(
                    "In get_vertex_edges, v_glob = {}, gid = {}, masked by masking_gid= {}, glob_coord = {}, nb_from = {}, nb_to = {}\n",
                    v_glob,
                    local.gid(), masking_gid, wrapped_neighb_vert_glob, nb_from, nb_to);
        }

        //assert(abs(nb_level - local.level()) <= 1);

        if (nb_level <= local.level())
        {
            // neighbour is on the same level or below me, neighbouring vertex corresponds
            // to the unique vertex in a neighbouring block
            size_t neighb_vertex_idx = get_vertex_id(wrapped_neighb_vert_glob, local.refinement(), link_idx, l,
                    local.mask_grid().c_order());

            result.emplace_back(v_glob_idx, reeber::AmrVertexId{masking_gid, neighb_vertex_idx});

            if (debug)
            {
                fmt::print("In get_vertex_edges, v_glob = {}, masking_gid = {}, Added edge to idx = {}\n", v_glob,
                        masking_gid,
                        neighb_vertex_idx);
            }

        } else if (nb_level > local.level())
        {
            Position masking_box_from, masking_box_to;
            std::tie(masking_box_from, masking_box_to) = refine_vertex(neighb_v_glob, local.refinement(),
                    nb_refinement);
            if (debug)
                fmt::print(
                        "In get_vertex_edges, v_glob = {}, masking_gid = {}, masking_box_from = {}, masking_box_to = {}\n",
                        v_glob, local.gid(), masking_box_from, masking_box_to);

            reeber::Box<D> masking_box(masking_box_from, masking_box_to);
            reeber::Box<D> covering_box(masking_box_from - Position::one(), masking_box_to + Position::one());

            for(const Position& masking_position : masking_box.positions())
            {
                for(Position covering_position_glob : covering_box.position_link(masking_position))
                {
                    Position covering_position_coarsened = wrap_point(
                            coarsen_point(covering_position_glob, nb_refinement,
                                    local.refinement()),
                            domain, local.refinement());
                    if (debug)
                    {
                        fmt::print(
                                "In get_vertex_edges, v_glob = {}, masking_gid = {}, cov_vert = {}, cov_vert_loc = {}\n",
                                v_glob, local.gid(), covering_position_glob,
                                covering_position_coarsened);
                    }

                    if (covering_position_coarsened == v_glob)
                    {

                        Position masking_position_global = wrap_point(masking_position, domain, nb_refinement);
                        r::AmrVertexId masking_vertex_idx{masking_gid,
                                                          get_vertex_id(masking_position_global, nb_refinement,
                                                                  link_idx, l,
                                                                  local.mask_grid().c_order())};

                        if (debug)
                        {
                            fmt::print(
                                    "IN get_vertex_edges, v_glob = {}, masking_gid = {}, masking_position = {}, nb_from = {}, nb_to = {}, adding edge {}\n",
                                    v_glob, local.gid(), masking_position, nb_from, nb_to,
                                    r::AmrEdge{v_glob_idx, masking_vertex_idx});
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
    bool debug = false; //gid == 0 or gid == 1;
    if (debug) fmt::print("Enter compute_outgoing_edges, gid = {}\n", gid);

    auto receivers = link_unique(l, gid);
    if (debug)
        fmt::print("In compute_outgoing_edges for block = {}, link size = {}, unique = {}\n", gid, l->size(),
                receivers.size());

    for(const Vertex& v_glob : local_.active_global_positions())
    {
        auto v_idx = local_.get_vertex_from_global_position(v_glob);
//        debug = (gid == 0 or gid == 1) and (v_idx == 64 or v_idx == 4486);
//        debug = (v_idx.gid == 0 and v_idx.vertex == 64) or (v_idx.gid == 1 and v_idx.vertex == 4486);
//        debug = (v_glob[0] == 0 and v_glob[1] == 0 and v_glob[2] == 6);

        AmrEdgeContainer out_edges = get_vertex_edges(v_glob, local_, l, domain(), debug);

        if (debug)
        {
            for(auto&& e : out_edges)
            {
                fmt::print("outogoing edge e = {} {}\n", std::get<0>(e), std::get<1>(e));
            }
        }

        if (not out_edges.empty())
        {
            vertex_to_outgoing_edges[local_.get_vertex_from_global_position(v_glob)] = out_edges;
            std::copy(out_edges.begin(), out_edges.end(), std::back_inserter(initial_edges_));
            for(const AmrEdge& e : out_edges)
            {
                assert(std::get<0>(e).gid == gid);
                if (debug) fmt::print("in compute_outgoing_edges, adding {} to gid_to_outoging_edges\n", e);
                gid_to_outgoing_edges_[std::get<1>(e).gid].push_back(e);
            }
        }
    }
}

//template<class Real, unsigned D>
//void FabTmtBlock<Real, D>::adjust_original_gids(int sender_gid, FabTmtBlock::GidVector& edges_from_sender)
//{
//    auto iter = std::find(original_link_gids_.begin(), original_link_gids_.end(), sender_gid);
//    bool i_talk_to_sender = iter != original_link_gids_.end();
//    bool sender_talks_to_me =
//            std::find(edges_from_sender.begin(), edges_from_sender.end(), gid) != edges_from_sender.end();
//    if (i_talk_to_sender and not sender_talks_to_me)
//        original_link_gids_.erase(iter);
//}

// delete edges from this block that end in a low vertex of a neighbor block.
// edges_from_gid: outogoing eddes we receive from a neigbor block
// subtract them from the edges we keep in gid_to_outgoing_edges_
template<class Real, unsigned D>
void FabTmtBlock<Real, D>::delete_low_edges(int sender_gid, FabTmtBlock::AmrEdgeContainer& edges_from_sender)
{
    bool debug = false; //gid == 1 or gid == 100;

    auto iter = gid_to_outgoing_edges_.find(sender_gid);
    if (iter == gid_to_outgoing_edges_.end())
    {
        // all edges from neighbor must end in low vertex of mine
        if (debug)
            fmt::print("In delete_low_edges in block with gid = {}, sender = {}, no edges from this block, exiting\n",
                    gid, sender_gid);
        return;  // we don't expect any edges from gid
    }

    // edges from neighbor come reversed, correct that
    std::transform(edges_from_sender.begin(), edges_from_sender.end(), edges_from_sender.begin(),
            &reeber::reverse_amr_edge);

//    if (debug) fmt::print("in FabTmtBlock::delete_low_edges, transform OK, n_edges_from_sender = {}\n", edges_from_sender.size());

    // put edges into set to find the set difference
    std::set<AmrEdge> my_edges{iter->second.begin(), iter->second.end()};
    std::set<AmrEdge> neighbor_edges{edges_from_sender.begin(), edges_from_sender.end()};

    int old_n_edges = my_edges.size();

    iter->second.clear();

//    if (debug) fmt::print("In delete_low_edges in block with gid = {}, sender = {}, clear OK, my_edges.size = {}, neighbor_edges.size = {}\n", gid, sender_gid, my_edges.size(), neighbor_edges.size());

    //if (debug)
    //{
    //    for(auto e : my_edges)
    //    {
    //        fmt::print("my gid = {}, sender_gid = {}, my_edge = {}\n", gid, sender_gid, e);
    //    }
    //    for(auto e : neighbor_edges)
    //    {
    //        fmt::print("my gid = {}, sender_gid = {}, edge_from_sender = {}\n", gid, sender_gid, e);
    //    }
    // }

    std::set_intersection(my_edges.begin(), my_edges.end(), neighbor_edges.begin(), neighbor_edges.end(),
            std::back_inserter(iter->second));

    int new_n_edges = iter->second.size();

    if (iter->second.empty())
        gid_to_outgoing_edges_.erase(iter);
    if (debug)
        fmt::print(
                "in Block::delete_low_edges, gid = {}, sender_gid = {}, erase OK, old_n_edges = {}, new_n_edges = {}\n",
                gid, sender_gid, old_n_edges, new_n_edges);
}

template<class Real, unsigned D>
void FabTmtBlock<Real, D>::adjust_outgoing_edges()
{
    bool debug = false;

    size_t s = initial_edges_.size();
    initial_edges_.clear();

    for(const auto& gid_edge_vector_pair : gid_to_outgoing_edges_)
    {
        std::copy(gid_edge_vector_pair.second.begin(), gid_edge_vector_pair.second.end(),
                std::back_inserter(initial_edges_));
    }
    std::sort(initial_edges_.begin(), initial_edges_.end());

    std::set<int> neighbor_gids;
    for(const AmrEdge& e : initial_edges_)
    {
        neighbor_gids.insert(std::get<1>(e).gid);
    }

    int orig_link_old_size = original_link_gids_.size();

    if (debug)
        fmt::print(
                "In adjust_outgoing_edges for gid = {}, old #edges = {}, new #edges = {}, old link size = {}, new link size = {}, new link = {}\n",
                gid, s, initial_edges_.size(), orig_link_old_size, original_link_gids_.size(),
                container_to_string(new_receivers_));

#ifdef AMR_MT_SEND_COMPONENTS
    for(Component& c : components_)
        c.set_edges(initial_edges_, original_vertex_to_deepest_);
#endif

    gid_to_outgoing_edges_.clear();

    // remove from original_vertex_to_deepest_ map
    // inner vertices, we don't need this information any more
    std::unordered_set<AmrVertexId> edge_vertices;
    for(const auto& e : initial_edges_)
    {
        assert(std::get<0>(e).gid == gid);
        edge_vertices.insert(std::get<0>(e));
    }

    for(auto iter = original_vertex_to_deepest_.cbegin(); iter != original_vertex_to_deepest_.cend();)
    {
        iter = (edge_vertices.count(iter->first) == 1) ? std::next(iter) : original_vertex_to_deepest_.erase(iter);
    }
}

template<class Real, unsigned D>
r::AmrVertexId FabTmtBlock<Real, D>::original_deepest(const AmrVertexId& v) const
{
    auto iter = original_vertex_to_deepest_.find(v);
    if (original_vertex_to_deepest_.cend() != iter)
        return iter->second;
    else
        throw std::runtime_error("Deepest not found for vertex");
}

template<class Real, unsigned D>
r::AmrVertexId FabTmtBlock<Real, D>::final_deepest(const AmrVertexId& v) const
{
    auto iter = final_vertex_to_deepest_.find(v);
    if (final_vertex_to_deepest_.cend() != iter)
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


//template<class Real, unsigned D>
//void FabTmtBlock<Real, D>::add_component_to_disjoint_sets(const AmrVertexId& deepest_vertex)
//{
//    //assert(components_disjoint_set_parent_.find(deepest_vertex) == components_disjoint_set_parent_.end());
//    //assert(components_disjoint_set_size_.find(deepest_vertex) == components_disjoint_set_size_.end());
//
//    components_disjoint_set_parent_[deepest_vertex] = deepest_vertex;
//    components_disjoint_set_size_[deepest_vertex] = 1;
//}


template<class Real, unsigned D>
void FabTmtBlock<Real, D>::create_component(const AmrVertexId& deepest_vertex)
{
    bool debug = false;
    if (debug) fmt::print("Entered create_component, gid = {}, deepest_vertex = {}\n", gid, deepest_vertex);

    set_original_deepest(deepest_vertex, deepest_vertex);

    if (debug)
        fmt::print("In create_component, gid = {}, deepest_vertex = {}, before emplace, size = {}\n", gid,
                deepest_vertex, components_.size());

    components_.emplace_back(deepest_vertex);

    if (debug)
        fmt::print("In create_component, gid = {}, deepest_vertex = {}, added to components, size = {}\n", gid,
                deepest_vertex, components_.size());

//    add_component_to_disjoint_sets(deepest_vertex);
}

template<class Real, unsigned D>
void FabTmtBlock<Real, D>::compute_original_connected_components(const VertexEdgesMap& vertex_to_outgoing_edges)
{
    bool debug = false;

    const auto& const_tree = current_merge_tree_;

    if (debug) fmt::print("compute_original_connected_components called\n");

    VertexNeighborMap component_nodes;

    for(const auto& vert_neighb_pair : const_tree.nodes())
    {

        component_nodes.clear();

        if (debug)
            fmt::print("in compute_connected_component, gid = {} processing vertex = {}\n", gid,
                    vert_neighb_pair.first);

        component_nodes[vert_neighb_pair.first] = vert_neighb_pair.second;
        Neighbor u = vert_neighb_pair.second;

        if (!original_deepest_computed(u->vertex))
        {

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

            while(u_ != v and not original_deepest_computed(v->vertex))
            {

                visited_neighbors.push_back(u_->vertex);
                visited_neighbors.push_back(v->vertex);
                visited_neighbors.push_back(s->vertex);

                u_ = v;

                v = std::get<1>(u_->parent());
                s = std::get<0>(u_->parent());

                visited_neighbors.push_back(v->vertex);
                visited_neighbors.push_back(s->vertex);

            }

            AmrVertexId deepest_vertex = (original_deepest_computed(v)) ? original_deepest(v) : v->vertex;

            AmrEdgeContainer edges;
            for(const AmrVertexId& v : visited_neighbors)
            {
                auto find_iter = vertex_to_outgoing_edges.find(v);
                if (find_iter != vertex_to_outgoing_edges.end())
                {
                    edges.insert(edges.end(), find_iter->second.begin(), find_iter->second.end());
                }
            }

//            if (not original_deepest_computed(v))
//            {
//                if (debug) fmt::print("in compute_connected_compponent, gid = {}, creating new cc, deepest = {}\n", gid, deepest_vertex);
//                create_component(deepest_vertex);
//            }

            for(const AmrVertexId& v : visited_neighbors)
            {
                if (debug)
                    fmt::print(
                            "in compute_connected_compponent, gid = {}, deepest = {}, setting to visited neighbor {}\n",
                            gid, deepest_vertex, v);
                set_original_deepest(v, deepest_vertex);
            }
        }
    }

    //assert(std::accumulate(const_tree.nodes().cbegin(), const_tree.nodes().cend(), true,
    //                       [this](const bool& prev, const typename VertexNeighborMap::value_type& vn) {
    //                           return prev and this->original_deepest_computed(vn.second);
    //                       }));
#ifdef AMR_MT_SEND_COMPONENTS
    current_vertex_to_deepest_ = original_vertex_to_deepest_;

    // copy nodes from local merge tree of block to merge trees of components
    for (const auto& vertex_deepest_pair : original_vertex_to_deepest_)
    {
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

    for(const auto& vertex_deepest_pair : original_vertex_to_deepest_)
    {
        original_deepest_.insert(vertex_deepest_pair.second);
    }
    current_deepest_ = original_deepest_;
}

template<class Real, unsigned D>
void FabTmtBlock<Real, D>::compute_final_connected_components()
{
    bool debug = true;

    const auto& const_tree = current_merge_tree_;

    if (debug) fmt::print("compute_final_connected_components called\n");

    for(const auto& vert_neighb_pair : const_tree.nodes())
    {
//        debug = (vert_neighb_pair.second->vertex == reeber::AmrVertexId { 1, 32005});

        if (debug)
            fmt::print("in compute_final_connected_component, gid = {} processing vertex = {}\n", gid,
                    vert_neighb_pair.first);

        Neighbor u = vert_neighb_pair.second;

        if (!final_deepest_computed(u->vertex))
        {
            if (debug)
                fmt::print("in compute_final_connected_compponent, gid = {}, deepest not computed, traversing\n", gid);
            // nodes we traverse while going down,
            // they all should get deepest node
            std::vector<AmrVertexId> visited_neighbors;

            Neighbor u_ = u;
            Neighbor v = std::get<1>(u_->parent());
            Neighbor s = std::get<0>(u_->parent());
            visited_neighbors.push_back(u->vertex);
            visited_neighbors.push_back(v->vertex);
            visited_neighbors.push_back(s->vertex);

            while(u_ != v and not final_deepest_computed(v->vertex))
            {

                // nodes we traverse while going down,
                visited_neighbors.push_back(u_->vertex);
                visited_neighbors.push_back(v->vertex);
                visited_neighbors.push_back(s->vertex);

                u_ = v;

                v = std::get<1>(u_->parent());
                s = std::get<0>(u_->parent());

                visited_neighbors.push_back(v->vertex);
                visited_neighbors.push_back(s->vertex);

            }

            AmrVertexId deepest_vertex = (final_deepest_computed(v)) ? final_deepest(v) : v->vertex;

            if (debug)
                fmt::print("in compute_final_connected_compponent, v = {}, deepest vertex = {}\n", v->vertex,
                        deepest_vertex);

            if (not final_deepest_computed(v))
            {
                final_vertex_to_deepest_[deepest_vertex] = deepest_vertex;
            }

            for(const AmrVertexId& v : visited_neighbors)
            {
                if (debug)
                    fmt::print(
                            "in compute_final_connected_compponent, gid = {}, deepest = {}, setting to visited neighbor {}\n",
                            gid, deepest_vertex, v);
                final_vertex_to_deepest_[v] = deepest_vertex;
            }
        }
    }
}

template<class Real, unsigned D>
bool FabTmtBlock<Real, D>::edge_exists(const AmrEdge& e) const
{
    return current_merge_tree_.contains(std::get<0>(e)) and current_merge_tree_.contains(std::get<1>(e));
}

template<class Real, unsigned D>
bool FabTmtBlock<Real, D>::edge_goes_out(const AmrEdge& e) const
{
    return current_merge_tree_.contains(std::get<0>(e)) xor current_merge_tree_.contains(std::get<1>(e));
}


//template<class Real, unsigned D>
//bool FabTmtBlock<Real, D>::are_components_connected(const AmrVertexId& deepest_a, const AmrVertexId& deepest_b)
//{
//    bool result = find_component_in_disjoint_sets(deepest_a) == find_component_in_disjoint_sets(deepest_b);
//    // for debug - assume everything is connected
//    return result;
//}


//template<class Real, unsigned D>
//bool FabTmtBlock<Real, D>::is_component_connected_to_any_internal(const FabTmtBlock::AmrVertexId& deepest)
//{
//    for (const auto& cc : components_)
//        if (are_components_connected(deepest, cc.root_))
//            return true;
//    // for debug - assume everything is connected
//    //assert(false);
//    return false;
//}

//template<class Real, unsigned D>
//r::AmrVertexId FabTmtBlock<Real, D>::find_component_in_disjoint_sets(AmrVertexId x)
//{
//    bool debug = false;
//    if (debug) fmt::print("Entered find_component_in_disjoint_sets, x = {}\n", x);
//    while (components_disjoint_set_parent_.at(x) != x)
//    {
//        AmrVertexId next = components_disjoint_set_parent_[x];
//        if (debug) fmt::print("in find_component_in_disjoint_sets, x = {}, next = {}\n", x, next);
//        components_disjoint_set_parent_[x] = components_disjoint_set_parent_.at(next);
//        x = next;
//    }
//    if (debug) fmt::print("Exiting find_component_in_disjoint_sets, x = {}\n", x);
//    return x;
//}

//template<class Real, unsigned D>
//void FabTmtBlock<Real, D>::connect_components(const FabTmtBlock::AmrVertexId& deepest_a,
//                                              const FabTmtBlock::AmrVertexId& deepest_b)
//{
//    AmrVertexId a_root = find_component_in_disjoint_sets(deepest_a);
//    AmrVertexId b_root = find_component_in_disjoint_sets(deepest_b);
//
//    if (a_root == b_root)
//        return;
//
//    auto a_size = components_disjoint_set_size_.at(a_root);
//    auto b_size = components_disjoint_set_size_.at(b_root);
//    if (a_size < b_size)
//    {
//        std::swap(a_root, b_root);
//        std::swap(a_size, b_size);
//    }
//
//    components_disjoint_set_parent_[b_root] = a_root;
//    components_disjoint_set_size_[a_root] += b_size;
//}


template<class Real, unsigned D>
Real FabTmtBlock<Real, D>::scaling_factor() const
{
    Real result = 1;
    for(unsigned i = 0; i < D; ++i)
        result /= refinement();
    return result;
}

template<class Real, unsigned D>
int FabTmtBlock<Real, D>::is_done_simple(const std::vector<FabTmtBlock::AmrVertexId>& vertices_to_check)
{
    // check that all edge outgoing from the current region
    // do not start in any component of our original local tree

//        bool debug = (gid == 0);
    bool debug = true;

    if (n_active_ == 0)
        return 1;

    if (debug)
        fmt::print("Enter is_done_simple, gid = {}, #vertices = {}, round = {}\n",
                gid, vertices_to_check.size(), round_);
    std::set<AmrVertexId> cur_deepest_of_original_deepest;
    for(const AmrVertexId v : original_deepest_)
    {
        if (current_merge_tree_.contains(v))
        {
            Neighbor nv = current_merge_tree_[v];
            AmrVertexId cur_deepest = current_merge_tree_.find_deepest(nv)->vertex;
            cur_deepest_of_original_deepest.insert(cur_deepest);
        } else
        {
            fmt::print("ALARM is_done_simple, gid = {}, orginal deepst =  {} not found in tree\n", gid, v);
            throw std::runtime_error("Here");
        }
    }

    if (debug)
        fmt::print("is_done_simple, gid = {}, #vertices = {}, round = {}, cur_deepest_of_original_deepest = {}\n", gid,
                vertices_to_check.size(), round_, cur_deepest_of_original_deepest.size());
    for(const AmrVertexId& v : vertices_to_check)
    {
        Neighbor nv = current_merge_tree_[v];
        AmrVertexId dnv = current_merge_tree_.find_deepest(nv)->vertex;

        if (cur_deepest_of_original_deepest.count(dnv))
            return 0;
    }
    if (debug) fmt::print("is_done_simple, gid = {}, returning 1\n", gid);
    return 1;
}

template<class Real, unsigned D>
std::pair<Real, size_t> FabTmtBlock<Real, D>::get_local_stats() const
{
    std::pair<Real, size_t> result{0, 0};
    for(const auto& vertex_node_pair : get_merge_tree().nodes())
    {
        AmrVertexId vertex = vertex_node_pair.first;
        Neighbor u = vertex_node_pair.second;

        // ingore remains of degree 2 vertices in map
        if (vertex != u->vertex)
            continue;
        // ignore non-local vertices
        if (vertex.gid != gid)
            continue;

        result.second++;
        result.first += u->value;

        for(const auto& val_vertex_pair : u->vertices)
        {
            if (val_vertex_pair.second.gid != gid)
                continue;
            result.second++;
            result.first += val_vertex_pair.first;
        }
    }
    return result;
}

template<class Real, unsigned D>
void FabTmtBlock<Real, D>::compute_local_integral()
{
    local_integral_.clear();

    Real sf = scaling_factor();

    const bool negate = negate_;
    const auto& nodes = get_merge_tree().nodes();
    for(const auto& vertex_node_pair : nodes)
    {
        AmrVertexId current_vertex = vertex_node_pair.first;

        assert(current_vertex == vertex_node_pair.second->vertex);

        // save only information about local vertices
        if (current_vertex.gid != gid)
            continue;

        Node* current_node = vertex_node_pair.second;
        AmrVertexId root = final_vertex_to_deepest_[current_vertex];

        Real root_value = nodes.at(root)->value;
        local_integral_[root] += sf * current_node->value;
        for(const auto& value_vertex_pair : current_node->vertices)
        {
            assert(value_vertex_pair.second.gid == gid);
            local_integral_[root] += sf * value_vertex_pair.first;
        }
    }
}

#ifdef AMR_MT_SEND_COMPONENTS

template<class Real, unsigned D>
std::vector<r::AmrVertexId> FabTmtBlock<Real, D>::get_current_deepest_vertices() const
{
    std::set<AmrVertexId> result_set;
    for (const auto& vertex_deepest_pair : current_vertex_to_deepest_)
    {
        result_set.insert(vertex_deepest_pair.second);
    }
    std::vector<AmrVertexId> result(result_set.begin(), result_set.end());
    return result;
}

template<class Real, unsigned D>
typename FabTmtBlock<Real, D>::Component& FabTmtBlock<Real, D>::find_component(const AmrVertexId& deepest_vertex)
{
    for (Component& cc : components_)
    {
        if (cc.root_ == deepest_vertex)
        {
            return cc;
        }
    }
    throw std::runtime_error("Connnected component not found");
}

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
void FabTmtBlock<Real, D>::save(const void* b, diy::BinaryBuffer& bb)
{
    const FabTmtBlock* block = static_cast<const FabTmtBlock*>(b);

    diy::save(bb, block->gid);
//        diy::save(bb, block->local_);
    diy::save(bb, block->current_merge_tree_);
    diy::save(bb, block->original_tree_);
    //            diy::save(bb, block->components_);
    diy::save(bb, block->domain_);
    diy::save(bb, block->initial_edges_);
//    diy::save(bb, block->gid_to_outgoing_edges_);
    diy::save(bb, block->new_receivers_);
    diy::save(bb, block->processed_receivers_);
    diy::save(bb, block->original_link_gids_);
    diy::save(bb, block->negate_);
    diy::save(bb, block->original_vertex_to_deepest_);
//    diy::save(bb, block->components_disjoint_set_parent_);
//    diy::save(bb, block->components_disjoint_set_size_);
    diy::save(bb, block->round_);
}

template<class Real, unsigned D>
void FabTmtBlock<Real, D>::load(void* b, diy::BinaryBuffer& bb)
{
    FabTmtBlock* block = static_cast<FabTmtBlock*>(b);

    diy::load(bb, block->gid);
//        diy::load(bb, block->local_);
    diy::load(bb, block->current_merge_tree_);
    diy::load(bb, block->original_tree_);
    //            diy::load(bb, block->components_);
    diy::load(bb, block->domain_);
    diy::load(bb, block->initial_edges_);
//    diy::load(bb, block->gid_to_outgoing_edges_);
    diy::load(bb, block->new_receivers_);
    diy::load(bb, block->processed_receivers_);
    diy::load(bb, block->original_link_gids_);
    diy::load(bb, block->negate_);
    diy::load(bb, block->original_vertex_to_deepest_);
//    diy::load(bb, block->components_disjoint_set_parent_);
//    diy::load(bb, block->components_disjoint_set_size_);
    diy::load(bb, block->round_);
}

