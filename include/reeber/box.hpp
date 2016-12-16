#include <dlog/log.h>

template<unsigned D>
reeber::Box<D>
reeber::Box<D>::
intersect(const Box& other) const
{
    Position from, to;
    for (unsigned i = 0; i < D; ++i)
    {
        from[i] = std::max(from_[i], other.from_[i]);
        to[i]   = std::min(to_[i],   other.to_[i]);
    }

    return Box(g_, from, to);
}

template<unsigned D>
bool
reeber::Box<D>::
intersects(const Box& other) const
{
    for (unsigned i = 0; i < D; ++i)
        if (std::max(from_[i], other.from_[i]) > std::min(to_[i], other.to_[i]))
            return false;

    return true;
}

template<unsigned D>
bool
reeber::Box<D>::
contains(const Position& p) const
{
    for (unsigned i = 0; i < D; ++i)
        if (p[i] > to_[i] || p[i] < from_[i])
            return false;
    return true;
}

template<unsigned D>
bool
reeber::Box<D>::
boundary(const Position& p, bool degenerate) const
{
    for (unsigned i = 0; i < D; ++i)
    {
        if (degenerate && from_[i] == to_[i]) continue;
        if (p[i] == from_[i] || p[i] == to_[i])
            return true;
    }

    return false;
}

template<unsigned D>
reeber::Box<D>
reeber::Box<D>::
side(unsigned axis, bool upper) const
{
    Box res(*this);

    if (upper)
        res.from()[axis] = res.to()[axis];
    else
        res.to()[axis] = res.from()[axis];

    return res;
}

template<unsigned D>
void
reeber::Box<D>::
merge(const Box& other)
{
#if DEBUG
    unsigned from_count = 0, to_count = 0;
    for (unsigned i = 0; i < D; ++i)
    {
        if (from_[i] != other.from_[i]) ++from_count;
        if (to_[i]   != other.to_[i])   ++to_count;
    }
    AssertMsg(from_count <= 1 && to_count <= 1, "Can only merge boxes that agree on all but one coordinate");
#endif

    LOG_SEV(trace) << "Merging partitions: " << *this << " and " << other;

    for (unsigned i = 0; i < D; ++i)
    {
        from_[i] = std::min(from_[i], other.from_[i]);
        to_[i]   = std::max(to_[i],   other.to_[i]);
    }
}

/* Box::FreudenthalLinkIterator */
template<unsigned D>
class reeber::Box<D>::FreudenthalLinkIterator:
    public std::iterator<std::forward_iterator_tag, Position>
{
    using Parent = std::iterator<std::forward_iterator_tag, Position>;

    public:
        typedef     typename Parent::value_type                     value_type;
        typedef     typename Parent::difference_type                difference_type;
        typedef     typename Parent::reference                      reference;

                    FreudenthalLinkIterator(): loc_(0), dir_(0)     {}
                    FreudenthalLinkIterator(const Position& p, int loc = 0, int dir = 1):
                        p_(p), v_(p), loc_(loc), dir_(dir)          {}

        static FreudenthalLinkIterator
                    begin(const Position& p)                        { FreudenthalLinkIterator it(p); ++it; return it; }
        static FreudenthalLinkIterator
                    end(const Position& p)                          { return FreudenthalLinkIterator(p, 0, -1); }

        const Position&             operator*() const               { return v_; }
        const Position*             operator->() const              { return &v_; }

        FreudenthalLinkIterator&   operator++()                     { increment(); return *this; }
        FreudenthalLinkIterator    operator++(int)                  { FreudenthalLinkIterator it = *this; increment(); return it; }

        friend bool operator==(const FreudenthalLinkIterator& x, const FreudenthalLinkIterator& y)    { return x.v_ == y.v_; }
        friend bool operator!=(const FreudenthalLinkIterator& x, const FreudenthalLinkIterator& y)    { return x.v_ != y.v_; }

    private:
        void        increment();

    private:
        Position    p_, v_;
        int         loc_;
        int         dir_;
};

template<unsigned D>
void
reeber::Box<D>::FreudenthalLinkIterator::
increment()
{
    loc_ += dir_;
    if (loc_ == (1 << D))
    {
        dir_ = -1;
        loc_ += dir_;
    }

    for (unsigned i = 0; i < D; ++i)
        if (loc_ & (1 << i))
            v_[i] = p_[i] + dir_;
        else
            v_[i] = p_[i];
}
