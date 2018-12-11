template<class Real>
FabConnectedComponent<Real>::FabConnectedComponent()
{
}

template<class Real>
FabConnectedComponent<Real>::FabConnectedComponent(bool negate, const AmrVertexId& deepest, Real total_value,
                                                   Real deepest_value) :
        negate_(negate),
        global_deepest_(deepest),
        original_deepest_(deepest),
        global_integral_value_(total_value),
        original_integral_value_(total_value),
        global_deepest_value_(deepest_value),
        original_deepest_value_(deepest_value),
        current_neighbors_({deepest}),
        processed_neighbors_({deepest})
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
    assert(std::all_of(processed_neighbors_.begin(), processed_neighbors_.end(),
                       [this](const AmrVertexId& v) { return this->current_neighbors_.count(v) == 1; }));

    if (debug) fmt::print("is_done, gid = {}, component = {}: # current_neighbors_ = {}, # processed_neighbors = {}\n", original_deepest_.gid, original_deepest_, current_neighbors_.size(), processed_neighbors_.size());

    return current_neighbors_.size() == processed_neighbors_.size();
}

template<class Real>
void FabConnectedComponent<Real>::set_global_deepest(const VertexValue& vv)
{
    // TODO: add assert with cmp. Add negate to component?
    global_deepest_value_ = vv.value;
    global_deepest_ = vv.vertex;
}

template<class Real>
void FabConnectedComponent<Real>::set_current_neighbors(const AmrVertexSet& new_current_neighbhors)
{
    for(AmrVertexId cn : current_neighbors_)
        assert(new_current_neighbhors.count(cn));
    current_neighbors_ = new_current_neighbhors;
}

template<class Real>
int FabConnectedComponent<Real>::must_send_to_gid(int gid) const
{
    for(const AmrVertexId& cn : current_neighbors_)
    {
        if (processed_neighbors_.count(cn))
            continue;
        if (cn.gid == gid)
            return 1;
    }
    return 0;
}

template<class Real>
void FabConnectedComponent<Real>::mark_gid_processed(int _gid)
{
    for(const auto& cn : current_neighbors_)
    {
        if (cn.gid == _gid)
        {
            processed_neighbors_.insert(cn);
        }
    }
}

template<class Real>
std::string FabConnectedComponent<Real>::to_string() const
{

    //fmt::format(
    //       "(negate = {}, global_deepest = {}, original_deepest = {}, global_integral_value_ = {}, original_integral_value_ = {}, global_deepest_value_ = {}, original_deepest_value = {}, current_neighbors = {}, processed_neighbors = {})",
    return fmt::format("(original_deepest = {}, original_deepest_value = {}, current_neighbors = {})",
                       original_deepest_,
                       original_deepest_value_, container_to_string(current_neighbors_));
}
