/**
 * Author: Dmitriy Morozov <dmitriy@mrzv.org>
 */

#ifndef DLOG_COUNTERS_H
#define DLOG_COUNTERS_H

#include <iostream>

namespace dlog
{

// Event is expected to provide a `static const char* name() { ... }` function,
// if an instance of the Counter class is to be output to std::ostream.
template<class Event>
struct Counter
{
#ifdef COUNTERS
    Counter&        operator++()                        { ++count; return *this; }
    Counter&        operator++(int)                     { count++; return *this; }
    Counter&        operator+=(int i)                   { count += i; return *this; }

                    operator size_t()                   { return count; }

    friend
    std::ostream&
    operator<<(std::ostream& out, const Counter&)       { out << Event::name() << ": " << count; return out; }

    static size_t   count;
#else
    Counter&        operator++()                        { return *this; }
    Counter&        operator++(int)                     { return *this; }
    Counter&        operator+=(int)                     { return *this; }

                    operator size_t()                   { return 0; }

    friend
    std::ostream&
    operator<<(std::ostream& out, const Counter&)       { out << Event::name() << ": " << "[no counters]"; return out; }
#endif
};

#ifdef COUNTERS
template<class Event>
size_t
Counter<Event>::count = 0;
#endif

#define COUNTER(Event)      dlog::Counter<Event>()

}

#endif
