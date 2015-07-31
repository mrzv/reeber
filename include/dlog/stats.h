/**
 * Author: Dmitriy Morozov <dmitriy@mrzv.org>
 */

#ifndef DLOG_STATS_H
#define DLOG_STATS_H

// Definitions that control log streams at compile time: PROFILE, DEBUG
// Setup:
//      prof.add_stream(out);
//      stats.open("stats.out");    // TODO: add ability to multiplex multiple streams (including log streams)
// Output macros:
//      prof << "event-name";       // enter
//      prof >> "event-name";       // exit
//      stats << "any info";

#include <fstream>
#include <iostream>
#include <vector>
#include <string>

#include <sstream>
#include <sys/time.h>
#include <cstdio>

#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#endif

namespace dlog
{

// Stats
inline
std::fstream&
get_stats_out()
{
    static std::fstream stats_out;
    return stats_out;
}

static std::fstream& stats = get_stats_out();

// Timing
typedef             unsigned long                       time_type;
inline std::string  clock_to_string(time_type time);
inline time_type    get_time();

struct Timer
{
    typedef         time_type                           time;
    typedef         time_type                           duration;

                    Timer()                             { restart(); }
    void            restart()                           { start = get_time(); }
    duration        elapsed() const                     { return get_time() - start; }

    time            start;
};

// Profiling
#ifdef PROFILE
struct Profiler
{
    typedef     time_type                               time;

    struct Event
    {
            Event(const std::string& name_, bool begin_):
                name(name_),
                begin(begin_),
                stamp(get_time())
                                                        {}

        std::string     name;
        bool            begin;
        time            stamp;
    };

    typedef     std::vector<Event>                      EventsVector;

            Profiler()                                  { reset_time(); }
           ~Profiler()                                  { output_events(); }

    void    add_stream(std::ostream& out)               { outs.push_back(&out); }
    void    reset_time()                                { start = get_time(); }

    void    operator<<(std::string name)                { enter(name); }
    void    operator>>(std::string name)                { exit(name); }

    inline
    void    enter(std::string name);
    inline
    void    exit(std::string name);
    void    flush()
    {
        output_events();
        for (size_t i = 0; i < outs.size(); ++i)
            outs[i]->flush();
    }

    private:
    inline
    void    output_events();

    private:
        time            start;
        std::vector<std::ostream*>   outs;
        EventsVector    events;
};
#else
struct Profiler
{
    void    add_stream(std::ostream& out)               {}
    void    reset_time()                                {}

    void    operator<<(std::string name)                {}
    void    operator>>(std::string name)                {}

    void    enter(const std::string& name)              {}
    void    exit(const std::string& name)               {}
    void    flush()                                     {}
};
#endif

inline
Profiler&
get_profiler()
{
    static Profiler profiler;
    return profiler;
}

static Profiler& prof = get_profiler();

struct scoped
{
    scoped(const std::string& name_):
        name(name_)                         { prof.enter(name); }
    ~scoped()                               { prof.exit(name); }

    std::string name;
};

#define PROF        dlog::prof

} // dlog


std::string
dlog::
clock_to_string(time_type time)
{
    char buf[13];       // +1 for the trailing null
    sprintf(buf, "%02d:%02d:%02d.%03d",
            (unsigned) time/1000/60/60,
            (unsigned) time/1000/60 % 60,
            (unsigned) time/1000 % 60,
            (unsigned) time % 1000);
    return buf;
}

dlog::time_type
dlog::
get_time()
{
#ifdef __MACH__ // OS X does not have clock_gettime, use clock_get_time
    clock_serv_t cclock;
    mach_timespec_t ts;
    host_get_clock_service(mach_host_self(), CALENDAR_CLOCK, &cclock);
    clock_get_time(cclock, &ts);
    mach_port_deallocate(mach_task_self(), cclock);
#else
    timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
#endif
    return ts.tv_sec*1000 + ts.tv_nsec/1000000;
}

#ifdef PROFILE
void
dlog::Profiler::
enter(std::string name)
{
    events.push_back(Event(name, true));
#ifdef DEBUG
    flush();
#endif
}

void
dlog::Profiler::
exit(std::string name)
{
    events.push_back(Event(name, false));
#ifdef DEBUG
    flush();
#endif
}

void
dlog::Profiler::
output_events()
{
    for (size_t j = 0; j < outs.size(); ++j)
    {
        std::ostream& out = *outs[j];
        for (size_t i = 0; i < events.size(); ++i)
        {
            const Event& e = events[i];
            out << clock_to_string(e.stamp - start) << (e.begin ? " <" : " >") << e.name << '\n';
        }
    }
    events.clear();
}
#endif // PROFILE

#endif // DLOG_STATS_H
