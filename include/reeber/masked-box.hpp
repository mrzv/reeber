#include <dlog/log.h>

// TODO: think about methods inherited from Box

template<unsigned D>
bool reeber::MaskedBox<D>::core_contains_global(const Position& p) const
{
    for (unsigned i = 0; i < D; ++i)
        if (p[i] > core_to_[i] || p[i] < core_from_[i]) {
            return false;
        }
    return true;
}

template<unsigned D>
bool reeber::MaskedBox<D>::is_on_boundary(const Position& p_global) const
{
    bool result = false;
    for (unsigned i = 0; i < D; ++i)
    {
        if (p_global[i] > bounds_to_[i] || p_global[i] < bounds_from_[i])
        {
            return false;
        }
        if (p_global[i] == bounds_to_[i] or p_global[i] == bounds_from_[i])
        {
            result = true;
        }
    }
    return result;
}

template<unsigned D>
bool reeber::MaskedBox<D>::is_strictly_inside(const Position& p_global) const
{
    for (unsigned i = 0; i < D; ++i)
    {
        if (p_global[i] >= core_to_[i] || p_global[i] <= core_from_[i])
        {
            return false;
        }
    }
    return true;
}

template<unsigned D>
bool reeber::MaskedBox<D>::bounds_contains_global(const Position& p) const
{
    for (unsigned i = 0; i < D; ++i)
        if (p[i] > bounds_to_[i] || p[i] < bounds_from_[i]) {
            return false;
        }
    return true;
}
template<unsigned D>
bool reeber::MaskedBox<D>::contains_local(const Position& p_local) const
{
    Position p_core = p_local - ghost_adjustment_;
    for (unsigned i = 0; i < D; ++i)
        if (p_core[i] < 0 || p_core[i] >= core_shape_[i]) {
            return false;
        }
    return true;
}

template<unsigned D>
void reeber::MaskedBox<D>::save(const void* mb, diy::BinaryBuffer& bb)
{
}

template<unsigned D>
void reeber::MaskedBox<D>::load(void* mb, diy::BinaryBuffer& bb)
{
}


/* MaskedBox::FreudenthalLinkIterator */
template<unsigned D>
class reeber::MaskedBox<D>::FreudenthalLinkIterator :
        public std::iterator<std::forward_iterator_tag, Position>
{
    using Parent = std::iterator<std::forward_iterator_tag, Position>;

public:
    typedef typename Parent::value_type value_type;
    typedef typename Parent::difference_type difference_type;
    typedef typename Parent::reference reference;

    FreudenthalLinkIterator() :
            loc_(0), dir_(0)
    {
    }

    FreudenthalLinkIterator(const Position& p, int loc = 0, int dir = 1) :
            p_(p), v_(p), loc_(loc), dir_(dir)
    {
    }

    static FreudenthalLinkIterator
    begin(const Position& p)
    {
        FreudenthalLinkIterator it(p);
        ++it;
        return it;
    }

    static FreudenthalLinkIterator
    end(const Position& p)
    {
        return FreudenthalLinkIterator(p, 0, -1);
    }

    const Position& operator*() const
    {
        return v_;
    }

    const Position* operator->() const
    {
        return &v_;
    }

    FreudenthalLinkIterator& operator++()
    {
        increment();
        return *this;
    }

    FreudenthalLinkIterator operator++(int)
    {
        FreudenthalLinkIterator it = *this;
        increment();
        return it;
    }

    friend bool operator==(const FreudenthalLinkIterator& x, const FreudenthalLinkIterator& y)
    {
        return x.v_ == y.v_;
    }

    friend bool operator!=(const FreudenthalLinkIterator& x, const FreudenthalLinkIterator& y)
    {
        return x.v_ != y.v_;
    }

private:
    void increment();

private:
    Position p_, v_;
    int loc_;
    int dir_;
};

template<unsigned D>
void
reeber::MaskedBox<D>::FreudenthalLinkIterator::
increment()
{
    loc_ += dir_;
    if (loc_ == (1 << D)) {
        dir_ = -1;
        loc_ += dir_;
    }

    for (unsigned i = 0; i < D; ++i)
        if (loc_ & (1 << i)) {
            v_[i] = p_[i] + dir_;
        } else {
            v_[i] = p_[i];
        }
}
