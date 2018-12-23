
template<class Real, unsigned D>
bool FabComponentBlock<Real, D>::cmp(Real a, Real b) const
{
    if (negate_)
        return a > b;
    else
        return a < b;
}


template<class Real, unsigned D>
void FabComponentBlock<Real, D>::set_mask(const diy::Point<int, D>& v_bounds,
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
        is_low = not is_ghost and
                 not cmp(fab_(v_bounds), rho);   //(negate_ ? fab_(v_bounds) < rho : fab_(v_bounds) > rho);

    r::AmrVertexId v_idx;
    if (not is_ghost) v_idx = local_.get_vertex_from_global_position(local_.global_position_from_local(v_bounds));

    if (debug)
    {
        fmt::print("gid = {}, in set_mask, v_bounds = {}, v_idx = {}, value = {}, is_ghost = {}\n", debug_gid, v_bounds,
                   v_idx,
                   fab_(v_bounds), is_ghost);
    }

    // initialization, actual mask to be set later
    if (is_ghost)
    {
        if (debug)
        { fmt::print("in set_mask, gid = {}, v_bounds = {}, GHOST detected\n", debug_gid, v_bounds); }
        local_.set_mask(v_bounds, MaskedBox::GHOST);
    } else
    {
        local_.set_mask(v_bounds, MaskedBox::ACTIVE);
    }

    const int v_ref = local_.refinement();
    const int v_level = local_.level();

    Position v_glob = wrap_point(v_bounds + local_.bounds_from(), domain_, v_ref);

    if (debug)
    {
        fmt::print("in set_mask, gid = {}, unwrapped v_glob = {}, wrapped = {}, v_idx = {}\n", local_.gid(),
                   v_bounds + local_.bounds_from(), v_glob, v_idx);
    }


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
            local_.set_mask(v_bounds, l->target(i).gid);
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
                local_.set_mask(v_bounds, l->target(i).gid);
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
                    fmt::print("Is masked SAME {} , is_ghost = {}, gid = {}, v_idx = {}\n", l->target(i).gid, is_ghost,
                               local_.gid(), v_idx);
                }
                mask_set = true;
                local_.set_mask(v_bounds, l->target(i).gid);
                break;
            }
        }
    }

    if (not mask_set and is_low)
    {
        local_.set_mask(v_bounds, MaskedBox::LOW);
    }

    if (not is_absolute_threshold and not mask_set and not is_ghost)
    {
        // we need to store local sum and local number of unmasked vertices in a block
        // and use this later to mark low vertices
        n_unmasked_++;
        sum_ += fab_(v_bounds);
    }

    if (debug)
    {
        fmt::print("in set_mask, final mask = {}, gid = {}, v_bounds = {},  v_idx = {}\n",
                   local_.pretty_mask_value(v_bounds), local_.gid(), v_bounds, v_idx);
    }
}

template<class Real, unsigned D>
void FabComponentBlock<Real, D>::set_low(const diy::Point<int, D>& v_bounds,
                                         const Real& absolute_threshold)
{
    if (local_.mask(v_bounds) != MaskedBox::ACTIVE)
        return;
    bool is_low = !cmp(fab_(v_bounds),
                       absolute_threshold); //   negate_ ? fab_(v_bounds) < absolute_threshold : fab_(v_bounds) > absolute_threshold;
    if (is_low)
        local_.set_mask(v_bounds, MaskedBox::LOW);
}

template<class Real, unsigned D>
void FabComponentBlock<Real, D>::init(Real absolute_rho, diy::AMRLink* amr_link)
{
    bool debug = false;
    std::string debug_prefix = "In FabComponentBlock::init, gid = " + std::to_string(gid);

    diy::for_each(local_.mask_shape(), [this, absolute_rho](const Vertex& v) {
        this->set_low(v, absolute_rho);
    });

    VertexEdgesMap vertex_to_outgoing_edges;
    compute_outgoing_edges(amr_link, vertex_to_outgoing_edges);
    compute_original_connected_components(vertex_to_outgoing_edges);

    // TODO: delete this? we are going to overwrite this in adjust_outgoing_edges anyway
//    for(int i = 0; i < amr_link->size(); ++i)
//    {
//        if (amr_link->target(i).gid != gid)
//        {
//            new_receivers_.insert(amr_link->target(i).gid);
//            original_link_gids_.push_back(amr_link->target(i).gid);
//        }
//    }

    if (debug)
        fmt::print("{}, constructed, refinement = {}, level = {}, local = {}, domain.max = {}, #components = {}\n",
                   debug_prefix, refinement(), level(), local_, domain().max, components_.size());
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

        int gid = local.mask(neighb_v_bounds);

        size_t link_idx = 0;
        bool link_idx_found = false;
        for(; link_idx < (size_t) l->size(); ++link_idx)
        {
            if (l->target(link_idx).gid == gid)
            {
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

        if (debug)
        {
            fmt::print(
                    "In get_vertex_edges, v_glob = {}, gid = {}, masked by gid= {}, glob_coord = {}, nb_from = {}, nb_to = {}\n",
                    v_glob,
                    local.gid(), gid, wrapped_neighb_vert_glob, nb_from, nb_to);
        }

        //assert(abs(nb_level - local.level()) <= 1);

        if (nb_level <= local.level())
        {
            // neighbour is on the same level or below me, neighbouring vertex corresponds
            // to the unique vertex in a neighbouring block
            size_t neighb_vertex_idx = get_vertex_id(wrapped_neighb_vert_glob, local.refinement(), link_idx, l,
                                                     local.mask_grid().c_order());

            result.emplace_back(v_glob_idx, reeber::AmrVertexId{gid, neighb_vertex_idx});

            if (debug)
            {
                fmt::print("In get_vertex_edges, v_glob = {}, gid = {}, Added edge to idx = {}\n", v_glob, local.gid(),
                           neighb_vertex_idx);
            }

        } else if (nb_level > local.level())
        {
            Position masking_box_from, masking_box_to;
            std::tie(masking_box_from, masking_box_to) = refine_vertex(neighb_v_glob, local.refinement(),
                                                                       nb_refinement);
            if (debug)
                fmt::print("In get_vertex_edges, v_glob = {}, gid = {}, masking_box_from = {}, masking_box_to = {}\n",
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
                        fmt::print("In get_vertex_edges, v_glob = {}, gid = {}, cov_vert = {}, cov_vert_loc = {}\n",
                                   v_glob, local.gid(), covering_position_glob,
                                   covering_position_coarsened);
                    }

                    if (covering_position_coarsened == v_glob)
                    {

                        Position masking_position_global = wrap_point(masking_position, domain, nb_refinement);
                        r::AmrVertexId masking_vertex_idx{gid, get_vertex_id(masking_position_global, nb_refinement,
                                                                             link_idx, l,
                                                                             local.mask_grid().c_order())};

                        if (debug)
                        {
                            fmt::print(
                                    "IN get_vertex_edges, v_glob = {}, gid = {}, masking_position = {}, nb_from = {}, nb_to = {}, adding edge {}\n",
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
//        debug = (v_glob[0] == 0 and v_glob[1] == 0 and v_glob[2] == 6);

        AmrEdgeContainer out_edges = get_vertex_edges(v_glob, local_, l, domain());

        if (debug)
            for(auto&& e : out_edges)
                fmt::print("outogoing edge e = {} {}\n", std::get<0>(e), std::get<1>(e));

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
        if (debug) fmt::print("In delete_low_edges in block with gid = {}, sender = {}, no edges from this block, exiting\n", gid, sender_gid);
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
        throw std::runtime_error("Bug in delete_low_edges, cannot have more edges after deletion");
    if (debug) fmt::print("in delete_low_edges, erase OK, old_n_edges = {}, new_n_edges = {}\n", old_n_edges, new_n_edges);

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

        disjoint_sets_.make_component_if_not_exists(other_deepest);
        disjoint_sets_.unite_components(my_deepest, other_deepest);

        Component& my_component = get_component_by_deepest(my_deepest);
        my_component.current_neighbors_.insert(other_deepest);

        if (debug) fmt::print("gid = {}, added to current_neighbors of  my_deepest = {} and other_deeepest= {}\n", gid, my_deepest, other_deepest);
    }
}

template<class Real, unsigned D>
void FabComponentBlock<Real, D>::adjust_outgoing_edges()
{
    // TODO: rename procedure
    for(Component& c : components_)
    {
        for(AmrVertexId v : c.current_neighbors_)
        {
            c.current_gids_.insert(v.gid);
        }
    }
}


template<class Real, unsigned D>
void FabComponentBlock<Real, D>::compute_original_connected_components(
        const FabComponentBlock::VertexEdgesMap& vertex_to_outgoing_edges)
{
    bool debug = false;

    auto vertices_ = local_.vertices();

    std::vector<AmrVertexId> vertices(std::begin(vertices_), std::end(vertices_));


    std::unordered_set<AmrVertexId> global_processed_vertices;
    for(const auto& v : vertices)
    {
        if (global_processed_vertices.count(v))
            continue;

        std::unordered_set<AmrVertexId> component_vertices;
        AmrEdgeContainer component_edges;
        std::stack<AmrVertexId> vertices_to_process;
        vertices_to_process.push(v);

        Real deepest_value = fab_(v);
        Real integral_val = 0.0;
        AmrVertexId deepest = v;

        while(!vertices_to_process.empty())
        {
            AmrVertexId u = vertices_to_process.top();
            vertices_to_process.pop();


            if (global_processed_vertices.count(u))
                continue;

            component_vertices.insert(u);
            global_processed_vertices.insert(u);
            if (vertex_to_outgoing_edges.count(u))
                std::copy(vertex_to_outgoing_edges.at(u).begin(), vertex_to_outgoing_edges.at(u).end(),
                          std::back_inserter(component_edges));

            Real u_val = fab_(u);

            integral_val += u_val;


            if (cmp(u_val, deepest_value))
            {
                deepest_value = u_val;
                deepest = u;
            }

            for(auto&& w : local_.link(u))
            {
                vertices_to_process.push(w);
            }
        }

        for(auto&& component_vertex : component_vertices)
        {
            vertex_to_deepest_[component_vertex] = deepest;
        }

        integral_val *= scaling_factor();

        components_.emplace_back(negate_, deepest, deepest_value);

        disjoint_sets_.make_component(deepest);
        original_integral_values_[deepest] = integral_val;
    }


    std::vector<AmrVertexId> vv_check { {7, 34947}, {1, 32005}, {2, 12113}, {7, 34947}, {1, 32005},
                                    {2, 12113}, {3, 4773}, {1, 16980}, {1, 32005}, {3, 4773},
                                    {1, 16980}, {0, 33913}, {1, 34350} };

    for(AmrVertexId v : vv_check)
    {
        if (v.gid != gid)
            continue;

//        fmt::print("DEBUG v = {}, deepest = {}\n", v, vertex_to_deepest_.at(v));

        auto deepest = vertex_to_deepest_.at(v);
        Component& c = get_component_by_deepest(deepest);
//        fmt::print("DEBUG c = {}\n", c);
    }

}

template<class Real, unsigned D>
bool FabComponentBlock<Real, D>::edge_exists(const AmrEdge& e) const
{
    return disjoint_sets_.has_component(std::get<0>(e)) and disjoint_sets_.has_component(std::get<1>(e));
}


template<class Real, unsigned D>
bool FabComponentBlock<Real, D>::edge_goes_out(const AmrEdge& e) const
{
    return disjoint_sets_.has_component(std::get<0>(e)) xor disjoint_sets_.has_component(std::get<1>(e));
}

//template<class Real, unsigned D>
//bool FabComponentBlock<Real, D>::is_component_connected_to_any_internal(const AmrVertexId& deepest)
//{
//    for(const auto& cc : components_)
//        if (disjoint_sets_.are_connected(deepest, cc.root_))
//            return true;
//    // for debug - assume everything is connected
//    //assert(false);
//    return false;
//}


template<class Real, unsigned D>
Real FabComponentBlock<Real, D>::scaling_factor() const
{
    Real result = 1;
    for(unsigned i = 0; i < D; ++i)
        result /= refinement();
    return result;
}


//template<class Real, unsigned D>
//int FabComponentBlock<Real, D>::is_done_simple(const std::vector<FabComponentBlock::AmrVertexId>& vertices_to_check)
//{
//    // check that all edge outgoing from the current region
//    // do not start in any component of our original local tree
//
////        bool debug = (gid == 0);
//    bool debug = false;
//
//    if (debug)
//        fmt::print("is_done_simple, gid = {}, #vertices = {}, round = {}\n", gid, vertices_to_check.size(), round_);
//    for(const AmrVertexId& v : vertices_to_check)
//    {
//        if (is_component_connected_to_any_internal(v))
//        {
//            if (debug) fmt::print("is_done_simple, gid = {}, v = {}, returning 0\n", gid, v);
//            return 0;
//        }
//    }
//    if (debug) fmt::print("is_done_simple, gid = {}, returning 1\n", gid);
//    return 1;
//}



template<class Real, unsigned D>
int FabComponentBlock<Real, D>::get_n_components_for_gid(int _gid) const
{
    return std::count_if(components_.begin(), components_.end(),
                         [_gid](const Component& c) { return c.must_send_to_gid(_gid); });
}

template<class Real, unsigned D>
int FabComponentBlock<Real, D>::are_all_components_done() const
{
    bool debug = false;
    if (debug) fmt::print("are_all_components_done, gid = {}\n", gid);

    for(const Component& c : components_)
        if (not c.is_done())
            return 0;

    if (debug) fmt::print("are_all_components_done, gid = {}, returning 1\n", gid);
    return 1;
}

template<class Real, unsigned D>
bool FabComponentBlock<Real, D>::check_symmetry(int other_gid,
                                                const std::vector<FabComponentBlock::Component>& received_components)
{
    bool debug = false;
    // check that all edges (outer_vertex -> my_vertex) have corresponding (my_vertex->outer_vertex) edge
    for(const Component& rc : received_components)
    {
        for(AmrVertexId received_deepest : rc.current_neighbors())
        {
            if (received_deepest.gid == gid)
            {
                Component& my_component = get_component_by_deepest(received_deepest);
                if (my_component.current_neighbors().count(received_deepest) == 0)
                {
                    throw std::runtime_error("Asymmetry-1");
//                    return false;
                }
            }
        }
    }

    for(const Component& my_component : components_)
    {
        AmrVertexId my_deepest = my_component.original_deepest();
        for(AmrVertexId other_deepest : my_component.current_neighbors())
        {
            if (other_deepest.gid == other_gid)
            {
                auto iter = find_if(received_components.begin(), received_components.end(),
                                    [other_deepest](const Component& other_component) {
                                        return other_component.original_deepest() == other_deepest;
                                    });
                if (iter == received_components.end())
                {
                    fmt::print("Cannot find component {}, this = {}", other_deepest, container_to_string(components_));
                    throw std::runtime_error("Asymmetry-2");
//                    return false;
                }
                if (iter->current_neighbors().count(my_deepest) == 0)
                {
                    fmt::print("Error, my_deepest = {}, other_deepest = {}, my_component = {}", my_deepest, other_deepest, my_component);
                    throw std::runtime_error("Asymmetry-3");
//                    return false;
                }
            }
        }
    }

    if (debug) fmt::print("gid = {}, symmetry OK\n", gid);
    return true;
}

template<class Real, unsigned D>
void FabComponentBlock<Real, D>::compute_local_integral(Real theta)
{
    global_integral_.clear();
    for(const Component& c : components_)
    {
        if (c.global_deepest() != c.original_deepest())
            continue;

        if (not cmp(c.global_deepest_value(), theta))
            continue;

        for(const auto& d : disjoint_sets_.component_of(c.global_deepest()))
        {
            global_integral_[c.global_deepest()] += original_integral_values_.at(d);
        }
    }
}



template<class Real, unsigned D>
void FabComponentBlock<Real, D>::save(const void* b, diy::BinaryBuffer& bb)
{
    const FabComponentBlock* block = static_cast<const FabComponentBlock*>(b);

    diy::save(bb, block->gid);
    diy::save(bb, block->domain_);
    diy::save(bb, block->negate_);
    diy::save(bb, block->round_);
}

template<class Real, unsigned D>
void FabComponentBlock<Real, D>::load(void* b, diy::BinaryBuffer& bb)
{
    FabComponentBlock* block = static_cast<FabComponentBlock*>(b);

    diy::load(bb, block->gid);
    diy::load(bb, block->domain_);
    diy::load(bb, block->negate_);
    diy::load(bb, block->round_);
}

