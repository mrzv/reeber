#ifndef REEBER_SMALL_SET_H
#define REEBER_SMALL_SET_H

#include <vector>
#include <algorithm>
#include <diy/serialization.hpp>

template<class T>
struct SmallSet
{
    using ContainerType = typename std::vector<T>;
    using value_type = typename ContainerType::value_type;
    using iterator = typename ContainerType::iterator;
    using const_iterator = typename ContainerType::const_iterator;
    using size_type = typename ContainerType::size_type;

    SmallSet() {};
    SmallSet(const SmallSet& s) : data_(s.data_) {};
    SmallSet(std::initializer_list<T> init) : data_(init) {};

    ContainerType data_;

    iterator begin() { return  data_.begin(); }
    const_iterator begin() const { return  data_.begin(); }

    iterator end() { return  data_.end(); }
    const_iterator end() const { return  data_.end(); }

    size_type size() const { return data_.size(); }

    auto count(const T& val) const { return std::count(data_.begin(), data_.end(), val); }

    void insert(const T& val) { if (0 == count(val)) data_.push_back(val); }

    template<class IterType>
    void insert(IterType first, IterType last)
    {
        for(IterType iter = first; iter != last; ++iter)
        {
            insert(*iter);
        }
    }

    void clear() { data_.clear(); }

    friend diy::Serialization<SmallSet<T>>;
};


namespace diy {

    template<class R>
    struct Serialization<SmallSet<R>> {
        using SmallSetR = SmallSet<R>;

        static void save(BinaryBuffer& bb, const SmallSetR& x)
        {
            diy::save(bb, x.data_);
        }

        static void load(BinaryBuffer& bb, SmallSetR& x)
        {
            diy::load(bb, x.data_);
        }
    };
}


#endif //REEBER_SMALL_SET_H
