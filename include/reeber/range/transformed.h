#pragma once

#include <iterator>
#include <type_traits>

#include "utility.h"

namespace reeber
{
namespace range
{


template<class Range, class Transform>
struct transformed_types
{
    using RangeIterator     = decltype(std::declval<Range>().begin());
    using RangeValue        = typename std::remove_reference<decltype(*std::declval<RangeIterator>())>::type;
    using value_type        = typename std::remove_reference<decltype(std::declval<Transform>()(std::declval<RangeValue>()))>::type;
    using IteratorParent    = std::iterator<std::forward_iterator_tag, value_type>;
};

template<class Range, class Transform>
struct transformed_iterator: public transformed_types<Range, Transform>::IteratorParent
{
    using iterator      = transformed_iterator;
    using value_type    = typename transformed_types<Range, Transform>::value_type;
    using RangeIterator = typename transformed_types<Range, Transform>::RangeIterator;

                        transformed_iterator(RangeIterator it, Transform t):
                            it_(it), t_(t)                      {}

    value_type          operator*() const                       { return t_(*it_); }

    iterator&           operator++()                            { ++it_; return *this; }
    iterator            operator++(int)                         { iterator it = *this; ++it_; return it; }

    friend bool         operator==(const iterator& x, const iterator& y)    { return x.it_ == y.it_; }
    friend bool         operator!=(const iterator& x, const iterator& y)    { return x.it_ != y.it_; }

    private:
        RangeIterator       it_;
        Transform           t_;
};


template<class Range, class Transform>
struct transformed_range: public iterator_range<transformed_iterator<Range, Transform>>
{
    using iterator = transformed_iterator<Range, Transform>;
    using Parent   = iterator_range<iterator>;

                        transformed_range(Range& r, Transform t):
                            Parent(iterator(r.begin(), t), iterator(r.end(), t))    {}
};

template<class Range, class Transform>
transformed_range<Range, Transform>
make_transformed(Range& r, Transform t)
{
    return transformed_range<Range,Transform>(r,t);
}

// adaptor
template<class Transform>
struct transform_t
{
                        transform_t(Transform t):
                            t_(t)                                   {}

    Transform           t_;
};

template<class Transform>
transform_t<Transform>
transformed(Transform t)
{
    return transform_t<Transform>(t);
}

template<class Range, class Transform>
transformed_range<Range, Transform>
operator|(Range&& r, const transform_t<Transform>& t)
{
    return transformed_range<Range,Transform>(r, t.t_);
}

}
}
