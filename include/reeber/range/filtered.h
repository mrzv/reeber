#pragma once

#include <iterator>
#include <type_traits>

#include "utility.h"

namespace reeber
{
namespace range
{

template<class Range, class Filter>
struct filtered_types
{
    using RangeIterator     = decltype(std::declval<Range>().begin());
    using RangeValue        = typename std::remove_reference<decltype(*std::declval<RangeIterator>())>::type;
    using value_type        = RangeValue;
    using reference         = decltype(*std::declval<RangeIterator>());
    using IteratorParent    = std::iterator<std::forward_iterator_tag, value_type, std::ptrdiff_t, value_type*, reference>;
};

template<class Range, class Filter>
struct filtered_iterator: public filtered_types<Range, Filter>::IteratorParent
{
    using iterator      = filtered_iterator;
    using reference     = typename filtered_types<Range, Filter>::reference;
    using RangeIterator = typename filtered_types<Range, Filter>::RangeIterator;

                        filtered_iterator(RangeIterator it, RangeIterator end, Filter f):
                            it_(it), end_(end), f_(f)           { if (!f_(*it_)) increment(); }

    reference           operator*() const                       { return *it_; }

    iterator&           operator++()                            { increment(); return *this; }
    iterator            operator++(int)                         { iterator it = *this; increment(); return it; }

    friend bool         operator==(const iterator& x, const iterator& y)    { return x.it_ == y.it_; }
    friend bool         operator!=(const iterator& x, const iterator& y)    { return x.it_ != y.it_; }

    private:
        void            increment()                             { while (it_ != end_ && !f_(*(++it_))); }

    private:
        RangeIterator   it_;
        RangeIterator   end_;
        Filter          f_;
};

template<class Range, class Filter>
struct filtered_range: public iterator_range<filtered_iterator<Range, Filter>>
{
    using iterator = filtered_iterator<Range, Filter>;
    using Parent   = iterator_range<iterator>;

                        filtered_range(Range& r, Filter f):
                            Parent(iterator(r.begin(), r.end(), f), iterator(r.end(), r.end(), f))    {}
};

template<class Range, class Filter>
filtered_range<Range, Filter>
make_filtered(Range& r, Filter f)
{
    return filtered_range<Range,Filter>(r,f);
}

// adaptor
template<class Filter>
struct filter_t
{
                        filter_t(Filter f):
                            f_(f)                                   {}

    Filter          f_;
};

template<class Filter>
filter_t<Filter>
filtered(const Filter& f)
{
    return filter_t<Filter>(f);
}

template<class Range, class Filter>
filtered_range<Range, Filter>
operator|(Range&& r, const filter_t<Filter>& f)
{
    return filtered_range<Range,Filter>(r, f.f_);
}


}
}
