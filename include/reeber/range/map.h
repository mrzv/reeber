#pragma once

#include <iterator>
#include <type_traits>

#include "utility.h"

namespace reeber
{
namespace range
{

template<class Range>
struct map_keys_types
{
    using RangeIterator     = decltype(std::declval<Range>().begin());
    using value_type        = typename std::remove_reference<decltype(std::declval<RangeIterator>()->first)>::type;
    using IteratorParent    = std::iterator<std::forward_iterator_tag, value_type>;
};

template<class Range>
struct map_keys_iterator: public map_keys_types<Range>::IteratorParent
{
    using RangeIterator = typename map_keys_types<Range>::RangeIterator;
    using value_type    = typename map_keys_types<Range>::value_type;
    using iterator      = map_keys_iterator;

                        map_keys_iterator(RangeIterator it):
                            it_(it)                             {}

    const value_type&   operator*() const                       { return it_->first; }

    iterator&           operator++()                            { ++it_; return *this; }
    iterator            operator++(int)                         { iterator it = *this; ++it_; return it; }

    friend bool         operator==(const iterator& x, const iterator& y)    { return x.it_ == y.it_; }
    friend bool         operator!=(const iterator& x, const iterator& y)    { return x.it_ != y.it_; }

    private:
        RangeIterator   it_;
};

template<class Range>
struct map_keys_range: public iterator_range<map_keys_iterator<Range>>
{
    using iterator = map_keys_iterator<Range>;
    using Parent   = iterator_range<iterator>;

                        map_keys_range(Range& r):
                            Parent(iterator(r.begin()), iterator(r.end()))      {}
};

template<class Range>
struct map_values_types
{
    using RangeIterator     = decltype(std::declval<Range>().begin());
    using value_type        = typename std::remove_reference<decltype(std::declval<RangeIterator>()->second)>::type;
    using IteratorParent    = std::iterator<std::forward_iterator_tag, value_type>;
};

template<class Range>
struct map_values_iterator: public map_values_types<Range>::IteratorParent
{
    using RangeIterator = typename map_values_types<Range>::RangeIterator;
    using value_type    = typename map_values_types<Range>::value_type;
    using iterator      = map_values_iterator;

                        map_values_iterator(RangeIterator it):
                            it_(it)                             {}

    const value_type&   operator*() const                       { return it_->second; }

    iterator&           operator++()                            { ++it_; return *this; }
    iterator            operator++(int)                         { iterator it = *this; ++it_; return it; }

    friend bool         operator==(const iterator& x, const iterator& y)    { return x.it_ == y.it_; }
    friend bool         operator!=(const iterator& x, const iterator& y)    { return x.it_ != y.it_; }

    private:
        RangeIterator   it_;
};

template<class Range>
struct map_values_range: public iterator_range<map_values_iterator<Range>>
{
    using iterator = map_values_iterator<Range>;
    using Parent   = iterator_range<iterator>;

                        map_values_range(Range& r):
                            Parent(iterator(r.begin()), iterator(r.end()))      {}
};

template<class Range>
map_values_range<Range>
make_map_values(Range& r)
{
    return map_values_range<Range>(r);
}

template<class Range>
map_keys_range<Range>
make_map_keys(Range& r)
{
    return map_keys_range<Range>(r);
}

// adaptor
struct map_keys_t {};
constexpr map_keys_t map_keys {};

struct map_values_t {};
constexpr map_values_t map_values {};

template<class Range>
map_keys_range<Range>
operator|(Range&& r, const map_keys_t&)
{
    return map_keys_range<Range>(r);
}

template<class Range>
map_values_range<Range>
operator|(Range&& r, const map_values_t&)
{
    return map_values_range<Range>(r);
}

}
}
