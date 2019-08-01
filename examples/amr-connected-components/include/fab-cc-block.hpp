template<class Real, unsigned D>
FabComponentBlock<Real, D>::FabComponentBlock(diy::GridRef<Real, D>& fab_grid,
        std::vector<std::string>& extra_names,
        std::vector<diy::GridRef<Real, D>>& extra_grids,
        int _ref,
        int _level,
        const diy::DiscreteBounds& _domain,
        const diy::DiscreteBounds& bounds,
        const diy::DiscreteBounds& core,
        int _gid,
        diy::AMRLink* amr_link,
        Real rho,                                           // threshold for LOW value
        bool _negate,
        bool is_absolute_threshold)
        :
        gid(_gid),
        local_(project_point<D>(core.min), project_point<D>(core.max), project_point<D>(bounds.min),
                project_point<D>(bounds.max), _ref, _level, gid, fab_grid.c_order()),
        fab_(fab_grid.data(), fab_grid.shape(), fab_grid.c_order()),
        domain_(_domain),
        negate_(_negate),
        merge_tree_(negate_)
#ifdef REEBER_EXTRA_INTEGRAL
, extra_names_(extra_names),
extra_grids_(extra_grids)
#endif
{
    assert(extra_names.size() == extra_grids.size());
    bool debug = false;

    std::string debug_prefix = "FabComponentBlock ctor, gid = " + std::to_string(gid);

    if (debug) fmt::print("{} setting mask\n", debug_prefix);

    diy::for_each(local_.mask_shape(), [this, amr_link, rho, is_absolute_threshold](const Vertex& v) {
        this->set_mask(v, amr_link, rho, is_absolute_threshold);
    });

//    for(int i = 0; i < amr_link->size(); ++i)
//    {
//        max_gid_ = std::max(max_gid_, amr_link->target(i).gid);
//    }
//    local_.check_mask_validity(max_gid_);

    if (debug) fmt::print("{}, max_gid_ = {}, is_absolute_threshold = {}, valid mask\n", debug_prefix, max_gid_, is_absolute_threshold);

    if (is_absolute_threshold)
    {
        if (debug) fmt::print("{}, calling init(rho = {}, link)\n", debug_prefix, rho);
        init(rho, amr_link, false); // false means: don't call set_low, already did in set_mask
    }

}

template<class Real, unsigned D>
bool FabComponentBlock<Real, D>::cmp(Real a, Real b) const
{
    if (negate_)
    {
        return a > b;
    } else
    {
        return a < b;
    }
}

template<class Real, unsigned D>
void FabComponentBlock<Real, D>::set_mask(const diy::Point<int, D>& v_mask,
        diy::AMRLink* l,
        const Real& rho,
        bool is_absolute_threshold)
{
    using Position = diy::Point<int, D>;

    int debug_gid = local_.gid();

    bool debug = false;

    bool is_ghost = local_.is_outer(v_mask);

    size_t fab_index = local_.local_position_to_vertex(local_.local_position_from_mask(v_mask)).vertex;

    bool is_low = false;
    if (is_absolute_threshold)
    {
        is_low = not is_ghost and cmp(rho, fab_(fab_index));   //(negate_ ? fab_(v_mask) < rho : fab_(v_mask) > rho);
    }

    if (debug)
    {
        fmt::print("gid = {}, in set_mask, v_mask = {},  value = {}, is_ghost = {}\n", debug_gid, v_mask,
                fab_(fab_index), is_ghost);
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

    if (debug)
    {
        fmt::print(
                "in set_mask, gid = {}, unwrapped v_glob = {}, wrapped = {}, domain = [{} - {}], local = {}\n",
                local_.gid(),
                v_mask + local_.mask_from(), v_glob, domain_.min, domain_.max, local_);
    }

    bool mask_set{false};

    // loop over neighbouring blocks one level above
    for(int i = 0; i < l->size(); ++i)
        if (l->level(i) == v_level + 1 && neighbor_contains<D>(i, l, v_glob, v_ref))
        {
            if (debug)
            {
                fmt::print("Is masked ABOVE {} , is_ghost = {}, gid = {}\n", local_.gid(), is_ghost,
                        l->target(i).gid);
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
                    fmt::print("Is masked SAME {} , is_ghost = {}, gid = {}\n", l->target(i).gid, is_ghost,
                            local_.gid());
                }
                mask_set = true;
                local_.set_mask(v_mask, l->target(i).gid);
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
                    fmt::print("Is masked SAME {} , is_ghost = {}, gid = {}\n", l->target(i).gid, is_ghost,
                            local_.gid());
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
        n_low_++;
        sum_low_ += fab_(fab_index);
    } else if (not mask_set and is_absolute_threshold)
    {
        n_active_++;
        sum_active_ += fab_(fab_index);
    }

    if (not mask_set and not is_ghost)
    {
        // Important: we need to store local sum and local number of unmasked vertices in a block
        // and use this later to mark low vertices
        n_unmasked_++;
        sum_ += fab_(fab_index);
    }

    if (debug)
    {
        fmt::print("in set_mask, is_ghost = {}, final mask = {}, {}, gid = {}, v_mask = {}\n",
                is_ghost, local_.pretty_mask_value(v_mask), local_.mask(v_mask), local_.gid(), v_mask);
    }

    if (is_ghost and local_.mask(v_mask) == gid)
    {
        fmt::print("in set_mask, is_ghost = {}, final mask = {}, {}, gid = {}, v_mask = {}\n",
                is_ghost, local_.pretty_mask_value(v_mask), local_.mask(v_mask), local_.gid(), v_mask);
        fmt::print("Error, same gid in mask\n");
        throw std::runtime_error("Bad mask");
    }
}

template<class Real, unsigned D>
void FabComponentBlock<Real, D>::set_low(const diy::Point<int, D>& p_core,
        const Real& absolute_threshold)
{
    bool is_debug = false;
    std::string debug_prefix = fmt::format("in FabComponentBlock::set_low, p_core = {}, absolute_threshold = {} ",
            p_core, absolute_threshold);

    if (not local_.is_in_core(p_core))
    {
        if (is_debug) fmt::print("{}, not in core, exiting\n", debug_prefix);
        return;
    }

    Vertex p_mask = p_core + Vertex::one();
    if (local_.mask(p_mask) != MaskedBox::ACTIVE)
    {
        if (is_debug)
        {
            fmt::print("{}, p_mask = {}, mask_value = {} not ACTIVE, exiting\n", debug_prefix, p_mask,
                    local_.pretty_mask_value(p_mask));
        }
        return;
    }

    Vertex p_bounds = p_core + local_.ghost_adjustment();
    bool is_low = cmp(absolute_threshold,
            fab_(p_bounds)); //   negate_ ? fab_(v_bounds) < absolute_threshold : fab_(v_bounds) > absolute_threshold;

    if (is_low)
    {
        local_.set_mask(p_mask, MaskedBox::LOW);
        n_low_++;
        sum_low_ += fab_(p_bounds);
    } else
    {
        n_active_++;
        sum_active_ += fab_(p_bounds);
    }


    if (is_debug)
    {
        fmt::print(
                "{}, fab_(p_bounds = {}) = {}, is_low = {}, n_low = {}, n_active = {}, sum_low_ = {}, sum_active_ = {}\n",
                debug_prefix, p_bounds, fab_(p_bounds), is_low, sum_low_, sum_active_);
    }
}

template<class Real, unsigned D>
void FabComponentBlock<Real, D>::init(Real absolute_rho, diy::AMRLink* amr_link, bool must_set_low)
{
    bool debug = false;
    std::string debug_prefix = "In FabComponentBlock::init, gid = " + std::to_string(gid);

    if (debug) fmt::print("{}, absolute_rho = {}\n", debug_prefix, absolute_rho);

    if (must_set_low)
    {
        diy::for_each(local_.core_shape(), [this, absolute_rho](const Vertex& v_core) {
            this->set_low(v_core, absolute_rho);
        });
    }

    if (debug) fmt::print("{}, low set\n", debug_prefix, absolute_rho);

    reeber::compute_merge_tree2(merge_tree_, local_, fab_);

    if (debug) fmt::print("{}, rho = {}, merge tree computed, size = {}\n", debug_prefix, absolute_rho, merge_tree_.size());

    VertexEdgesMap vertex_to_outgoing_edges;
    compute_outgoing_edges(amr_link, vertex_to_outgoing_edges);

    compute_original_connected_components(vertex_to_outgoing_edges);

    if (debug)
    {
        fmt::print("{}, constructed, refinement = {}, level = {}, local = {}, domain.max = {}, #components = {}\n",
                debug_prefix, refinement(), level(), local_, domain().max, components_.size());
    }

    sparsify_prune_original_tree(vertex_to_outgoing_edges);

    if (debug)
    {
        int n_edges = 0;
        for(auto& gid_edges : vertex_to_outgoing_edges)
        { n_edges += gid_edges.second.size(); }
        fmt::print(
                "{},  constructed, tree.size = {}, n_edges = {},\n"
                " n_active = {}, n_low = {}, n_masked = {}, n_unmasked = {}, total_core_size = {}, diff_check = {}, diff_check_1 = {}\n"
                " sum_ = {}, sum_active_ = {}, sum_low_ = {}\n",
                debug_prefix, merge_tree_.size(), n_edges, n_active_, n_low_, n_masked_, n_unmasked_,
                local_.core_shape()[0] * local_.core_shape()[1] * local_.core_shape()[2],
                n_low_ + n_active_ - n_unmasked_, n_low_ + n_active_ + n_masked_
                        - local_.core_shape()[0] * local_.core_shape()[1] * local_.core_shape()[2],
                sum_, sum_active_, sum_low_);
    }

}

template<class Real, unsigned D>
void FabComponentBlock<Real, D>::sparsify_prune_original_tree(const VertexEdgesMap& vertex_to_outgoing_edges)
{
#ifndef REEBER_NO_SPARSIFICATION
    bool debug = false;
    size_t old_size = 0, new_size = 0;

    for(Component& c : components_)
    {
        old_size += c.tree_.size();
        c.sparsify(vertex_to_outgoing_edges);
        new_size += c.tree_.size();
    }

    if (debug)
    {
        fmt::print("gid = {}, MT size before sparsification = {}, after = {}\n", gid, old_size, new_size);
    }
#endif
}

template<unsigned D>
r::AmrEdgeContainer
get_vertex_edges(const diy::Point<int, D>& v_glob, const reeber::MaskedBox<D>& local, diy::AMRLink* l,
        const diy::DiscreteBounds& domain)
{
    using Position = diy::Point<int, D>;

    r::AmrVertexId v_glob_idx = local.get_vertex_from_global_position(v_glob);

    bool debug = false;

    if (debug)
    { fmt::print("get_vertex_edges called for {}, vertex {}, box = {}\n", v_glob, v_glob_idx, local); }

    r::AmrEdgeContainer result;

    for(const Position& neighb_v_glob : local.outer_edge_link(v_glob))
    {

        Position neighb_v_bounds = neighb_v_glob - local.bounds_from();

        if (debug)
        {
            fmt::print("In get_vertex_edges, v_glob = {}, gid = {}, neighb_v_bounds = {}, neighb_v_glob = {}\n", v_glob,
                    local.gid(), neighb_v_bounds, neighb_v_glob);
        }

        Position wrapped_neighb_vert_glob = wrap_point(neighb_v_glob, domain, local.refinement());

        Position neighb_v_mask = local.mask_position_from_local(neighb_v_bounds);
        int masking_gid = local.mask(neighb_v_mask);

        if (masking_gid == r::MaskedBox<D>::UNINIT)
        {
            fmt::print("v_glob = {}, v_glob_idx = {}, neighb_v_mask = {}\n", v_glob, v_glob_idx, masking_gid);
            throw std::runtime_error("bad masking gid");
        }

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
                    "In get_vertex_edges, v_glob = {}, gid = {}, masked by masking_gid= {}, nb_level = {}, my_level = {}, glob_coord = {}, nb_from = {}, nb_to = {}\n",
                    v_glob,
                    local.gid(), masking_gid, nb_level, local.level(), wrapped_neighb_vert_glob, nb_from, nb_to);
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
void FabComponentBlock<Real, D>::compute_outgoing_edges(diy::AMRLink* l, VertexEdgesMap& vertex_to_outgoing_edges)
{
    bool debug = false;
    for(const Vertex& v_glob : local_.active_global_positions())
    {
        AmrEdgeContainer out_edges = get_vertex_edges(v_glob, local_, l, domain());
        if (debug) for(auto&& e : out_edges) fmt::print("outogoing edge e = {} {}\n", std::get<0>(e), std::get<1>(e));

        if (not out_edges.empty())
        {
            vertex_to_outgoing_edges[local_.get_vertex_from_global_position(v_glob)] = out_edges;
            for(const AmrEdge& e : out_edges)
            {
                assert(std::get<0>(e).gid == gid);
                if (debug) fmt::print("in compute_outgoing_edges, adding {} to gid_to_outoging_edges\n", e);
                gid_to_outgoing_edges_[std::get<1>(e).gid].push_back(e);
            }
        }
    }
}

// delete edges from this block that end in a low vertex of a neighbor block.
// edges_from_gid: outogoing eddes we receive from a neigbor block
// subtract them from the edges we keep in gid_to_outgoing_edges_
template<class Real, unsigned D>
void FabComponentBlock<Real, D>::delete_low_edges(int sender_gid, AmrEdgeContainer& edges_from_sender,
        const VertexVertexMap& received_vertex_to_deepest)
{
    bool debug = false;

    auto iter = gid_to_outgoing_edges_.find(sender_gid);
    if (iter == gid_to_outgoing_edges_.end())
    {
        // all edges from neighbor must end in low vertex of mine
        if (debug)
        {
            fmt::print("In delete_low_edges in block with gid = {}, sender = {}, no edges from this block, exiting\n",
                    gid, sender_gid);
        }
        return;  // we don't expect any edges from gid
    }

    // edges from neighbor come reversed, correct that
    std::transform(edges_from_sender.begin(), edges_from_sender.end(), edges_from_sender.begin(),
            &reeber::reverse_amr_edge);

    if (debug) fmt::print("in FabComponentBlock::delete_low_edges, transform OK\n");

    // put edges into set to find the set difference
    std::set<AmrEdge> my_edges{iter->second.begin(), iter->second.end()};
    std::set<AmrEdge> neighbor_edges{edges_from_sender.begin(), edges_from_sender.end()};

    int old_n_edges = my_edges.size();

    iter->second.clear();

    if (debug) fmt::print("In delete_low_edges in block with gid = {}, sender = {}, clear OK\n", gid, sender_gid);

    std::set_intersection(my_edges.begin(), my_edges.end(), neighbor_edges.begin(), neighbor_edges.end(),
            std::back_inserter(iter->second));

    int new_n_edges = iter->second.size();
    assert(old_n_edges >= new_n_edges);
    if (old_n_edges < new_n_edges)
    {
        throw std::runtime_error("Bug in delete_low_edges, cannot have more edges after deletion");
    }
    if (debug)
    {
        fmt::print("in delete_low_edges, erase OK, old_n_edges = {}, new_n_edges = {}\n", old_n_edges, new_n_edges);
    }

    if (iter->second.empty())
    {
        gid_to_outgoing_edges_.erase(iter);
        if (debug) fmt::print("remove container from gid_to_outgoing_edges and exit\n");
        return;
    }

    for(auto&& e : iter->second)
    {
        if (debug) fmt::print("gid = {}, loop over common edges, consider e = {}\n", gid, e);

        AmrVertexId my_deepest = vertex_to_deepest_.at(std::get<0>(e));
        AmrVertexId other_deepest = received_vertex_to_deepest.at(std::get<1>(e));

        Component& my_component = get_component_by_deepest(my_deepest);
        my_component.add_current_neighbor(other_deepest);

        if (debug)
        {
            fmt::print("gid = {}, added to current_neighbors of  my_deepest = {} and other_deeepest= {}\n", gid,
                    my_deepest, other_deepest);
        }
    }
}

template<class Real, unsigned D>
void FabComponentBlock<Real, D>::adjust_outgoing_edges()
{
//    int total_n_edges = std::accumulate(gid_to_outgoing_edges_.begin(), gid_to_outgoing_edges_.end(), 0, [](int a, const auto& x) { return a + x.second.size(); });
//    std::vector<AmrEdge> all_edges;
//    all_edges.reserve(total_n_edges);
    for(const auto& int_set_pair : gid_to_outgoing_edges_)
    {
        for(const AmrEdge& e : int_set_pair.second)
        {
            AmrVertexId our_vertex = std::get<0>(e);
            AmrVertexId deepest = vertex_to_deepest_.at(our_vertex);
            Component& c = get_component_by_deepest(deepest);
            c.add_edge(e);
//            all_edges.push_back(e);
        }
    }

//    std::unordered_set<AmrVertexId> processed_components;
//    for(Component& c : components_)
//    {
//        if (processed_components.count(c.original_deepest()))
//        {
//            continue;
//        }
//        processed_components.insert(c.original_deepest());
//
//        std::unordered_set<AmrVertexId> all_neighbors;
//        auto all_comps = component_of(c.original_deepest());
//
//        for(const auto& d : all_comps)
//        {
//            if (d.gid != gid)
//            {
//                continue;
//            }
//            Component& x = get_component_by_deepest(d);
//            all_neighbors.insert(x.current_neighbors().begin(), x.current_neighbors().end());
//        }
//
//        for(const auto& d : all_comps)
//        {
//            if (d.gid != gid)
//            {
//                continue;
//            }
//            Component& x = get_component_by_deepest(d);
//            x.set_current_neighbors(all_neighbors);
//            processed_components.insert(d);
//        }
//    }

}

template<class Real, unsigned D>
bool FabComponentBlock<Real, D>::is_deepest_computed(const AmrVertexId& v) const
{
    return vertex_to_deepest_.count(v);
}
//
//template<class Real, unsigned D>
//void FabComponentBlock<Real, D>::set_deepest(const AmrVertexId& v, const AmrVertexId& deepest)
//{
//    vertex_to_deepest_[v] = deepest;
//}

template<class Real, unsigned D>
void FabComponentBlock<Real, D>::compute_original_connected_components(
        const FabComponentBlock::VertexEdgesMap& vertex_to_outgoing_edges)
{
#ifdef REEBER_EXTRA_INTEGRAL
    local_integral_.clear();
#endif

#ifdef DO_DETAILED_TIMING
    dlog::Timer timer;
    dlog::Timer copy_nodes_timer;
#endif

    Real sf = scaling_factor();

    [[maybe_unused]] bool debug = false;

    const TripletMergeTree& const_tree = merge_tree_;
    std::unordered_set<AmrVertexId> processed_deepest;

    for(const auto& vertex_neighbor_pair : const_tree.nodes())
    {

        AmrVertexId u = vertex_neighbor_pair.first;
        Neighbor n = vertex_neighbor_pair.second;

        Neighbor deepest_neighbor = merge_tree_.find_deepest(n);
        AmrVertexId deepest_vertex = deepest_neighbor->vertex;

        vertex_to_deepest_[u] = deepest_vertex;

        // do we need this?
        for(auto vv : n->vertices)
        {
            vertex_to_deepest_[vv.second] = deepest_vertex;
        }

#ifdef REEBER_EXTRA_INTEGRAL
        local_integral_[deepest_vertex]["n_vertices"] += 1 + n->vertices.size();
        local_integral_[deepest_vertex]["n_vertices_sf"] += sf * (1 + n->vertices.size());
        local_integral_[deepest_vertex]["function_value"] += sf * fab_(u);
        for(size_t i = 0; i < extra_names_.size(); ++i)
        {
            auto density_iter = std::find(extra_names_.begin(), extra_names_.end(), "density");

            size_t density_idx = density_iter - extra_names_.begin();

            if (extra_names_[i] == "xmom" or extra_names_[i] == "ymom" or extra_names_[i] == "zmom")
            {
                // we need velocities - divide momentum by density; assumes density is the first field
                if (density_iter == extra_names_.end())
                {
                    throw std::runtime_error("density not found");
                }

                local_integral_[deepest_vertex][extra_names_.at(i)] += sf * (extra_grids_.at(i)(u)
                        / extra_grids_[density_idx](u));
            } else
            {
                local_integral_[deepest_vertex][extra_names_.at(i)] += sf * extra_grids_.at(i)(u);
            }

            for(auto vvv : n->vertices)
            {
                AmrVertexId vv = vvv.second;
                if (extra_names_[i] == "xmom" or extra_names_[i] == "ymom" or extra_names_[i] == "zmom")
                {
                    local_integral_[deepest_vertex][extra_names_.at(i)] += sf * (extra_grids_.at(i)(vv)
                            / extra_grids_[density_idx](vv));
                } else
                {
                    local_integral_[deepest_vertex][extra_names_.at(i)] += sf * extra_grids_.at(i)(vv);
                }
            }
        }
#endif

        if (processed_deepest.count(deepest_vertex) == 0)
        {
            // we encounter this deepest vertex for the first time
            AmrVertexId deepest_value = deepest_neighbor->vertex;
            // local integral is still incomplete and will be set below,
            // after all active vertices were processed
            components_.emplace_back(negate_, deepest_vertex, deepest_value);
            processed_deepest.insert(deepest_vertex);
        }

#ifdef DO_DETAILED_TIMING
        copy_nodes_timer.restart();
#endif

        TripletMergeTree& mt = get_component_by_deepest(deepest_vertex).tree_;
        Value val = n->value;

        AmrVertexId s = std::get<0>(n->parent())->vertex;
        AmrVertexId v = std::get<1>(n->parent())->vertex;

        Neighbor n_u, n_s, n_v;

        n_u = mt.add_or_update(u, val);
        n_s = mt.find_or_add(s, 0);
        n_v = mt.find_or_add(v, 0);
        mt.link(n_u, n_s, n_v);

#ifdef DO_DETAILED_TIMING
        copy_nodes_time += copy_nodes_timer.elapsed();
#endif
    }

#ifdef DO_DETAILED_TIMING
    compute_components_time += timer.elapsed();
#endif

#ifdef REEBER_EXTRA_INTEGRAL
    for(auto& deepest_extra_values_pair : local_integral_)
    {
        AmrVertexId deepest = deepest_extra_values_pair.first;
        get_component_by_deepest(deepest).set_extra_values(deepest_extra_values_pair.second);

        if (debug)
        {
            for(auto ev : deepest_extra_values_pair.second)
            {
                LOG_SEV(info) << "In compute_original_components, gid = " << gid << ", field " << ev.first << " = "
                              << ev.second;
            }
        }

    }
#endif
}

// fill in vertex_to_deepest_ map with correct values
template<class Real, unsigned D>
void FabComponentBlock<Real, D>::compute_final_connected_components()
{
    bool debug = false;

    vertex_to_deepest_.clear();

    const auto& const_tree = merge_tree_;

    if (debug) fmt::print("compute_final_connected_components called, gid = {}\n", gid);

    for(const auto& vert_neighb_pair : const_tree.nodes())
    {
//        debug = (vert_neighb_pair.second->vertex == reeber::AmrVertexId { 1, 32005});

        if (debug)
            fmt::print("in compute_final_connected_component, gid = {} processing vertex = {}\n", gid,
                    vert_neighb_pair.first);

        Neighbor u = vert_neighb_pair.second;

        if (!is_deepest_computed(u->vertex))
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

            while(u_ != v and not is_deepest_computed(v->vertex))
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

            AmrVertexId deepest_vertex = (is_deepest_computed(v->vertex)) ? vertex_to_deepest_.at(v->vertex)
                                                                          : v->vertex;

            if (debug)
                fmt::print("in compute_final_connected_compponent, v = {}, deepest vertex = {}\n", v->vertex,
                        deepest_vertex);

            if (not is_deepest_computed(v->vertex))
            {
                vertex_to_deepest_[deepest_vertex] = deepest_vertex;
            }

            for(const AmrVertexId& v : visited_neighbors)
            {
                if (debug)
                    fmt::print(
                            "in compute_final_connected_compponent, gid = {}, deepest = {}, setting to visited neighbor {}\n",
                            gid, deepest_vertex, v);
                vertex_to_deepest_[v] = deepest_vertex;
            }
        }
    }
    if (debug) fmt::print("Exit compute_final_connected_components called, gid = {}\n", gid);
}

template<class Real, unsigned D>
void FabComponentBlock<Real, D>::update_connectivity(const AmrVertexContainer& deepest)
{
    bool debug = false;
    for(const Component& c : components_)
    {
        auto original_deepest = c.original_deepest();
        if (merge_tree_.contains(original_deepest))
        {
            auto current_deepest = merge_tree_.find_deepest(merge_tree_[original_deepest])->vertex;
            vertex_to_deepest_[original_deepest] = current_deepest;
        }
    }

    for(AmrVertexId v : deepest)
    {
        if (merge_tree_.contains(v))
        {
            auto current_deepest = merge_tree_.find_deepest(merge_tree_[v])->vertex;
            if (debug)
                fmt::print("in update_connectivity, gid = {}, v = {}, current_deepest = {}\n", gid, v, current_deepest);
            vertex_to_deepest_[v] = current_deepest;
        }
    }
}

template<class Real, unsigned D>
Real FabComponentBlock<Real, D>::scaling_factor() const
{
    Real result = 1;
    for(unsigned i = 0; i < D; ++i)
    {
        result /= refinement();
    }
    return result;
}

template<class Real, unsigned D>
int FabComponentBlock<Real, D>::are_all_components_done() const
{
    bool debug = false;
    if (debug) fmt::print("are_all_components_done, gid = {}\n", gid);

    for(const Component& c : components_)
    {
        if (not c.is_done_sending())
        {
            return 0;
        }
    }

    if (debug) fmt::print("are_all_components_done, gid = {}, returning 1\n", gid);
    return 1;
}

//TODO: fix
template<class Real, unsigned D>
bool FabComponentBlock<Real, D>::check_symmetry(int other_gid,
        const std::vector<FabComponentBlock::Component>& received_components)
{
//    bool debug = false;
    // check that all edges (outer_vertex -> my_vertex) have corresponding (my_vertex->outer_vertex) edge
//    for(const Component& rc : received_components)
//    {
//        for(int received_deepest_gid : rc.current_neighbors())
//        {
//            if (received_deepest_gid == gid)
//            {
//                Component& my_component = get_component_by_deepest(received_deepest);
//                if (my_component.current_neighbors().count(received_deepest) == 0)
//                {
//                    throw std::runtime_error("Asymmetry-1");
////                    return false;
//                }
//            }
//        }
//    }

//    for(const Component& my_component : components_)
//    {
//        AmrVertexId my_deepest = my_component.original_deepest();
//        for(AmrVertexId other_deepest : my_component.current_neighbors())
//        {
//            if (other_deepest.gid == other_gid)
//            {
//                auto iter = find_if(received_components.begin(), received_components.end(),
//                                    [other_deepest](const Component& other_component) {
//                                        return other_component.original_deepest() == other_deepest;
//                                    });
//                if (iter == received_components.end())
//                {
//                    fmt::print("Cannot find component {}, this = {}", other_deepest, container_to_string(components_));
//                    throw std::runtime_error("Asymmetry-2");
////                    return false;
//                }
//                if (iter->current_neighbors().count(my_deepest) == 0)
//                {
//                    fmt::print("Error, my_deepest = {}, other_deepest = {}, my_component = {}", my_deepest, other_deepest, my_component);
//                    throw std::runtime_error("Asymmetry-3");
////                    return false;
//                }
//            }
//        }
//    }

//    if (debug) fmt::print("gid = {}, symmetry OK\n", gid);
    return true;
}

template<class Real, unsigned D>
void FabComponentBlock<Real, D>::sanity_check_fin() const
{
//    std::unordered_set<AmrVertexId> components;
//    // check that all received components have integral value
//    for(const auto& x : disjoint_sets_.size_)
//    {
//        auto v = x.first;
//        components.insert(v);
//        assert(original_integral_values_.count(v));
//    }
//
//    for(const auto v : components)
//    {
//        assert(std::any_of(components_.begin(), components_.end(), [v, this](const Component& c) {
//            return this->disjoint_sets_.are_connected(v, c.original_deepest());
//        }));
//    }
}

template<class Real, unsigned D>
void FabComponentBlock<Real, D>::compute_local_integral()
{
#ifdef REEBER_EXTRA_INTEGRAL
    bool debug = false;

    if (debug) fmt::print("Enter compute_local_integral, gid = {}\n", gid);

    for(auto li_iter = local_integral_.begin(); li_iter != local_integral_.end();)
    {

        AmrVertexId v = li_iter->first;

        if (vertex_to_deepest_.count(v) == 0)
        {
            fmt::print("ERROR HERE, v not found, v= {}, gid = {}\n", v, gid);
            LOG_SEV_IF(true, info) << "ERROR HERE, v not found, v= " << v << ", gid = " << gid;
            dlog::flush();
        }

        AmrVertexId root = vertex_to_deepest_.at(v);

        LOG_SEV_IF(debug, info) << "In loop over local_integral, gid = " << gid << ", v = " << v << ", root = " << root;

        // if deepest vertex belongs to another block, skip it
        if (root != v)
        {
            if (local_integral_.count(root) == 0)
            {
                fmt::print("ERROR HERE, root not found, v= {}, root= {}, gid = {}\n", v, root, gid);
                LOG_SEV_IF(true, info) << "ERROR HERE, root not found, v= " << v << ", root = " << root << ", gid = "
                                       << gid;
                dlog::flush();
            }

            assert(local_integral_.count(root));
            if (root.gid == gid)
            {
                for(const auto& field_sum : li_iter->second)
                {
                    if (local_integral_.at(root).count(field_sum.first) == 0)
                    {
                        fmt::print("ERROR HERE, field not found, v= {}, root= {}, filed = {}, gid = {}\n", v, root,
                                field_sum.first, gid);
                        LOG_SEV(error) << "ERROR HERE, field not found, v= " << v << ", root = " << root
                                       << ", field = " << field_sum.first
                                       << ", gid = " << gid;
                        dlog::flush();
                    }

                    local_integral_.at(root).at(field_sum.first) += field_sum.second;
                    LOG_SEV_IF(debug, info) << "In loop over local_integral, gid = " << gid << ", v = " << v
                                            << ", root = " << root << ", adding to " << field_sum.first << ", value "
                                            << field_sum.second << ", new value = "
                                            << local_integral_.at(root).at(field_sum.first);
                }
            }
            li_iter = local_integral_.erase(li_iter);
        } else
        {
            ++li_iter;
        }
    }

    for(auto root_field_values : local_integral_)
    {
        for(auto field_val_pair : root_field_values.second)
        {
            LOG_SEV_IF(debug, info) << "Computed local_integral, gid = " << gid << ", root = "
                                    << root_field_values.first
                                    << ", field = " << field_val_pair.first
                                    << ", value " << field_val_pair.second;
        }
    }

#endif
}

template<class Real, unsigned D>
void FabComponentBlock<Real, D>::sparsify_local_tree(const AmrVertexSet& keep)
{
#ifndef REEBER_NO_SPARSIFICATION
    bool debug = false; //(gid % 100 == 1);
    [[maybe_unused]] auto old_size = merge_tree_.size();
    r::sparsify(merge_tree_,
            [this, &keep](AmrVertexId u) {
#ifdef REEBER_EXTRA_INTEGRAL
                return (u.gid == this->gid) or (this->local_integral_.find(u) != this->local_integral_.end()) or keep.find(u) != keep.end();
#else
                return u.gid == this->gid or keep.find(u) != keep.end();
#endif
            });
    [[maybe_unused]] auto new_size = merge_tree_.size();

    if (debug)
    { fmt::print("gid = {}, before sparsification {}, after {}\n", gid, old_size, new_size); }
#endif
}

template<class Real, unsigned D>
void FabComponentBlock<Real, D>::save(const void* b, diy::BinaryBuffer& bb)
{
    const FabComponentBlock* block = static_cast<const FabComponentBlock*>(b);

    diy::save(bb, block->gid);
    diy::save(bb, block->domain_);
    diy::save(bb, block->negate_);
    diy::save(bb, block->round_);
    diy::save(bb, block->vertex_to_deepest_);
    diy::save(bb, block->merge_tree_);
    diy::save(bb, block->local_diagrams_);
    diy::save(bb, block->components_);
}

template<class Real, unsigned D>
void FabComponentBlock<Real, D>::load(void* b, diy::BinaryBuffer& bb)
{
    FabComponentBlock* block = static_cast<FabComponentBlock*>(b);

    diy::load(bb, block->gid);
    diy::load(bb, block->domain_);
    diy::load(bb, block->negate_);
    diy::load(bb, block->round_);
    diy::load(bb, block->vertex_to_deepest_);
    diy::load(bb, block->merge_tree_);
    diy::load(bb, block->local_diagrams_);
    diy::load(bb, block->components_);
}

