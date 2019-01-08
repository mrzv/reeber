template<class Real>
FabConnectedComponent<Real>::FabConnectedComponent()
{
}

template<class Real>
FabConnectedComponent<Real>::FabConnectedComponent(bool negate, const AmrVertexId& deepest, Real deepest_value,
                                                   VertexValueMap *_integral_values) :
        negate_(negate),
        global_deepest_(deepest),
        original_deepest_(deepest),
        global_deepest_value_(deepest_value),
        current_neighbors_({deepest.gid}),
        processed_neighbors_({deepest.gid}),
        block_integral_values_(_integral_values),
        tree_(negate)
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
                       [this](const auto& gid) { return this->current_neighbors_.count(gid) == 1; }));

    bool result = current_neighbors_.size() == processed_neighbors_.size();

//    if (result)
//    {
//        for(const auto& v : block_disjoint_sets_->component_of(original_deepest()))
//        {
//            if (block_integral_values_->count(v) == 0)
//                assert(false);
//        }
//    }
    return result;
}

template<class Real>
void FabConnectedComponent<Real>::set_global_deepest(const VertexValue& vv)
{
    assert(not cmp(global_deepest_value_, vv.value));
    global_deepest_value_ = vv.value;
    global_deepest_ = vv.vertex;
}

template<class Real>
void FabConnectedComponent<Real>::add_current_neighbor(int new_current_neighbor)
{
    bool debug = false;
    current_neighbors_.insert(new_current_neighbor);
//    if (debug) { fmt::print("in add_current_neighbor for {}, added ncn = {}, processed_gids = {}\n", original_deepest(), new_current_neighbor, container_to_string(processed_gids_)); }
}

template<class Real>
void FabConnectedComponent<Real>::set_current_neighbors(const GidSet& new_current_neighbhors)
{
//    AmrVertexId debug_v {6, 15811};
//    AmrVertexId debug_v_1 {2, 13462};
//    bool debug = (original_deepest() == debug_v or original_deepest() == debug_v_1);
    bool debug = false;

//    if (debug) { fmt::print("in set_current_neighbors, old cn = {}, new cn = {}\n", container_to_string(current_neighbors_), container_to_string(new_current_neighbhors)); }
//    if (debug) { fmt::print("in set_current_neighbors, current_gids = {}, processed_gids = {}\n", container_to_string(current_gids()), container_to_string(processed_gids())); }

    current_neighbors_ = new_current_neighbhors;
}

template<class Real>
int FabConnectedComponent<Real>::must_send_tree_to_gid(int gid) const
{
    return (processed_neighbors_.count(gid) == 0 and current_neighbors_.count(gid) == 1);
}

//template<class Real>
//int FabConnectedComponent<Real>::must_send_neighbors_to_gid(int gid) const
//{
//    for(const auto& cn : current_neighbors_)
//    {
//        if (cn.gid == gid and processed_gids_.count(cn) == 0)
//        {
//            return true;
//        }
//    }
//    return false;
//}

template<class Real>
void FabConnectedComponent<Real>::mark_all_processed()
{
    assert(std::all_of(processed_neighbors_.begin(), processed_neighbors_.end(),
                       [this](const auto& v) { return this->current_neighbors_.count(v) == 1; }));

    processed_neighbors_ = current_neighbors_;
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
            "Component(negate = {}, original_deepest = {}, current_neighbors = {}, processed_neighbors = {}\n",
            negate_, original_deepest_, container_to_string(current_neighbors()),
            container_to_string(processed_neighbors()));
}
