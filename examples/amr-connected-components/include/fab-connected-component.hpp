template<class Real>
FabConnectedComponent<Real>::FabConnectedComponent()
{
}

template<class Real>
FabConnectedComponent<Real>::FabConnectedComponent(bool negate, const AmrVertexId& deepest, Real deepest_value) :
        negate_(negate),
        global_deepest_(deepest),
        original_deepest_(deepest),
        global_deepest_value_(deepest_value),
//        original_deepest_value_(deepest_value),
        current_neighbors_({deepest}),
        processed_neighbors_({deepest}),
        current_gids_({deepest.gid}),
        processed_gids_({deepest.gid})
{
}

template<class Real>
bool FabConnectedComponent<Real>::cmp(Real x, Real y) const
{
    return negate_ ? x > y : x < y;
}

template<class Real>
int FabConnectedComponent<Real>::is_done() const
{
    bool debug = false;

//    assert(std::all_of(processed_neighbors_.begin(), processed_neighbors_.end(),
//                       [this](const AmrVertexId& v) { return this->current_neighbors_.count(v) == 1; }));

    assert(std::all_of(processed_gids_.begin(), processed_gids_.end(),
                       [this](const auto& gid) { return this->current_gids_.count(gid) == 1; }));

//    if (debug) fmt::print("is_done, gid = {}, component = {}: # current_neighbors_ = {}, # processed_neighbors = {}\n", original_deepest_.gid, original_deepest_, current_neighbors_.size(), processed_neighbors_.size());
    if (debug) fmt::print("is_done, gid = {}, component = {}: # current_gids_ = {}, # processed_gids = {}\n", original_deepest_.gid, original_deepest_, current_gids_.size(), processed_gids_.size());

    return current_gids_.size() == processed_gids_.size();
}

template<class Real>
void FabConnectedComponent<Real>::set_global_deepest(const VertexValue& vv)
{
    assert(not cmp(global_deepest_value_, vv.value));
    global_deepest_value_ = vv.value;
    global_deepest_ = vv.vertex;
}

template<class Real>
void FabConnectedComponent<Real>::add_current_neighbor(const AmrVertexId& new_current_neighbor)
{
    current_gids_.insert(new_current_neighbor.gid);
    current_neighbors_.insert(new_current_neighbor);
}

template<class Real>
void FabConnectedComponent<Real>::set_current_neighbors(const AmrVertexSet& new_current_neighbhors)
{
    for(AmrVertexId cn : current_neighbors_)
    {
        assert(new_current_neighbhors.count(cn));
        current_gids_.insert(cn.gid);
    }

    current_neighbors_ = new_current_neighbhors;
}

template<class Real>
int FabConnectedComponent<Real>::must_send_to_gid(int gid) const
{
    return (processed_gids_.count(gid) == 0 and current_gids_.count(gid) == 1);
}

template<class Real>
void FabConnectedComponent<Real>::mark_gid_processed(int _gid)
{
    processed_gids_.insert(_gid);
}

template<class Real>
void FabConnectedComponent<Real>::mark_neighbor_processed(AmrVertexId v)
{
    processed_neighbors_.insert(v);
}

template<class Real>
std::string FabConnectedComponent<Real>::to_string() const
{

//    return fmt::format("(negate = {}, global_deepest = {}, original_deepest = {}, global_integral_value_ = {}, original_integral_value_ = {}, global_deepest_value_ = {}, original_deepest_value = {}, current_neighbors = {}, processed_neighbors = {})",
//    return fmt::format("Component(global_deepest = {}, original_deepest = {},  global_deepest_value_ = {}, original_deepest_value = {}",
//                        global_deepest_, original_deepest_, global_deepest_value_, original_deepest_value_);
    return fmt::format("Component(global_deepest = {}, global_deepest_value_ = {}",
                       global_deepest_, global_deepest_value_);
//    return fmt::format("(original_deepest = {}, original_deepest_value = {}, current_neighbors = {})", original_deepest_, original_deepest_value_, container_to_string(current_neighbors_));
}
