#pragma once

namespace reeber
{
namespace range
{

template<class Iterator>
struct iterator_range
{
    using iterator = Iterator;

                    iterator_range(Iterator b, Iterator e):
                        begin_(b), end_(e)
    {}

    iterator        begin() const           { return begin_; }
    iterator        end() const             { return end_; }

    iterator        begin_, end_;
};

template<class Iterator>
iterator_range<Iterator>
make_iterator_range(Iterator b, Iterator e)
{
    return iterator_range<Iterator>(b,e);
}

}
}
