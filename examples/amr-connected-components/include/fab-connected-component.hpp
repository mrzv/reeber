template<class Real>
FabConnectedComponent<Real>::FabConnectedComponent()
{
}

template<class Real>
FabConnectedComponent<Real>::FabConnectedComponent(bool negate, const AmrVertexId& deepest, Real deepest_value) :
        negate_(negate),
        original_deepest_(deepest),
        global_deepest_value_(deepest_value),
        current_neighbors_({deepest}),
        current_gids_({deepest.gid}),
        processed_gids_({deepest.gid}),
        n_prev_current_neighbors_(1),
        tree_(negate)
{
}


template<class Real>
bool FabConnectedComponent<Real>::cmp(Real x, Real y) const
{
    return negate_ ? x > y : x < y;
}

template<class Real>
bool FabConnectedComponent<Real>::is_done_sending() const
{
    [[maybe_unused]] bool debug = false;

    assert(std::all_of(processed_gids().begin(), processed_gids().end(),
                       [this](const auto& gid) { return this->current_gids_.count(gid) == 1; }));

    return current_gids().size() == processed_gids().size() and not must_send_neighbors();
}


template<class Real>
void FabConnectedComponent<Real>::add_current_neighbor(const AmrVertexId& new_current_neighbor)
{
    bool debug = false;
    current_neighbors_.insert(new_current_neighbor);
    current_gids_.insert(new_current_neighbor.gid);
    if (debug) { fmt::print("in add_current_neighbor for {}, added ncn = {}, processed_gids = {}\n", original_deepest(), new_current_neighbor, container_to_string(processed_gids_)); }
}

template<class Real>
void FabConnectedComponent<Real>::set_current_neighbors(const AmrVertexSet& new_current_neighbors)
{
    [[maybe_unused]] bool debug = false;

    assert(std::all_of(current_neighbors_.begin(), current_neighbors_.end(),
                       [new_current_neighbors](const auto& x) { return new_current_neighbors.count(x) == 1; }));

    n_prev_current_neighbors_ = current_neighbors_.size();

    current_neighbors_ = new_current_neighbors;

    for(const AmrVertexId& v : new_current_neighbors)
        current_gids_.insert(v.gid);
}

template<class Real>
bool FabConnectedComponent<Real>::must_send_to_gid(int gid) const
{
    return current_gids_.count(gid) == 1;
}

template<class Real>
bool FabConnectedComponent<Real>::must_send_tree_to_gid(int gid) const
{
//    bool debug = (gid == 0);
//    if (debug) fmt::print("HERE must_send_tree_to_gid({}), original_deepest = {}, processed_gids.count = {}, current_gids.count = {}", gid, original_deepest(), processed_gids().count(gid), current_gids_.count(gid));
    return (processed_gids_.count(gid) == 0 and current_gids_.count(gid) == 1);
}


template<class Real>
int FabConnectedComponent<Real>::must_send_neighbors() const
{
    assert(n_prev_current_neighbors_ <= current_neighbors_.size());
    return n_prev_current_neighbors_ < current_neighbors_.size();
}

template<class Real>
void FabConnectedComponent<Real>::mark_all_gids_processed()
{
    assert(std::all_of(processed_gids_.begin(), processed_gids_.end(),
                       [this](const auto& v) { return this->current_gids_.count(v) == 1; }));

    processed_gids_ = current_gids_;
}

template<class Real>
void FabConnectedComponent<Real>::add_edge(const AmrEdge& e)
{
    assert(std::find(edges_.begin(), edges_.end(), e) == edges_.end());
    edges_.push_back(e);
}

template<class Real>
std::string FabConnectedComponent<Real>::to_string() const
{
    return fmt::format(
            "Component(negate = {}, original_deepest = {}, current_gids = {}, processed_gids = {}\n",
            negate_, original_deepest_, container_to_string(current_gids()),
            container_to_string(processed_gids()));
}
