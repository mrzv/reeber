#pragma once

#ifdef REEBER_USE_TBB

#include <atomic>
#include <tbb/tbb.h>
#include <tbb/concurrent_vector.h>
#include <tbb/concurrent_unordered_map.h>
#include <tbb/scalable_allocator.h>

namespace reeber
{
    using tbb::task_scheduler_init;

    // atomic
    template<class T>
    using atomic = std::atomic<T>;

    template<class T>
    bool compare_exchange(atomic<T>& x, T& expected, T desired)     { return x.compare_exchange_weak(expected, desired); }

    // vector
    template<class T>
    using vector = tbb::concurrent_vector<T>;

    // foreach
    template<class Iterator, class F>
    void                do_foreach(Iterator begin, Iterator end, const F& f)        { tbb::parallel_do(begin, end, f); }

    template<class Range, class F>
    void                for_each_range_(const Range& r, const F& f)
    {
        for (auto& x : r)
            f(x);
    }

    template<class Container, class F>
    void                for_each_range(Container& c, const F& f)
    {
        tbb::parallel_for(c.range(), [&](const typename Container::range_type& r) { for_each_range_(r, f); });
    }

    template<class F>
    void                for_each(size_t from, size_t to, const F& f)                { tbb::parallel_for(from, to, f); }

    // map
    template<class Key, class T,
             class Hash = std::hash<Key>,
             class KeyEqual = std::equal_to<Key>,
             class Allocator = tbb::tbb_allocator<std::pair<Key, T>>>
    using map = tbb::concurrent_unordered_map<Key, T, Hash, KeyEqual, Allocator>;

    template<class... Args>
    using map_range = typename map<Args...>::range_type;

    template<class Key, class... Args>
    void map_erase(map<Key, Args...>& m, const Key& k)                              { m.unsafe_erase(k); }
}

#else

#include <vector>
#include <unordered_map>

namespace reeber
{
    struct task_scheduler_init
    {
                        task_scheduler_init(unsigned)   {}
        static const unsigned automatic = 0;
    };

    // atomic
    template<class T>
    using atomic = T;

    template<class T>
    bool compare_exchange(atomic<T>& x, T& expected, T desired)                     { x = desired; return true; }

    // vector
    template<class T>
    using vector = std::vector<T>;

    // foreach
    template<class Iterator, class F>
    void                do_foreach(Iterator begin, Iterator end, const F& f)        { std::for_each(begin, end, f); }

    template<class Container, class F>
    void                for_each_range(Container& c, const F& f)
    {
        for(auto& x : c)
            f(x);
    }

    template<class F>
    void                for_each(size_t from, size_t to, const F& f)                { for (size_t x = from; x < to; ++x) f(x); }

    // map
    template<class Key, class T,
             class Hash = std::hash<Key>,
             class KeyEqual = std::equal_to<Key>,
             class Allocator = std::allocator<std::pair<const Key, T>>>
    using map = std::unordered_map<Key, T, Hash, KeyEqual, Allocator>;

    template<class... Args>
    using map_range = map<Args...>;

    template<class Key, class... Args>
    void map_erase(map<Key, Args...>& m, const Key& k)                              { m.erase(k); }
}

#endif
