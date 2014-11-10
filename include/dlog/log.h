/**
 * Author: Dmitriy Morozov <dmitriy@mrzv.org>
 */

#ifndef DLOG_LOG_H
#define DLOG_LOG_H

// Definitions that control log streams at compile time: DEBUG, TRACE
// Setup:
//      add_stream(out, level);
// Output macros:
//      LOG_SEV(...) << ...;
//      LOG_SEV_IF(condition, ...) << ...;
//      AssertMsg(condition, ... << ...);
//      BREAK_IF(condition, ... << ...);


#include <iostream>
#include <vector>
#include <string>
#include <iomanip>
#include <unistd.h>

namespace dlog
{

enum severity_level
{
    trace,
    debug,
    info,
    warning,
    error,
    fatal,
    num_severities
};

inline
severity_level severity(const std::string&);

class StreamWrapper;

typedef unsigned                        Sinks[num_severities];
typedef std::vector<StreamWrapper*>     Streams;
typedef std::vector<std::string>        Scope;

struct ScopeFilterImpl
{
    virtual         ~ScopeFilterImpl()              {}
    virtual bool    operator()(const Scope& s) const =0;
    virtual ScopeFilterImpl*    clone() const =0;
};

template<class Object>
struct ScopeFilterImplConcrete: public ScopeFilterImpl
{
    Object o;
                    ScopeFilterImplConcrete(Object o_): o(o_)   {}
    virtual bool    operator()(const Scope& s) const    { return o(s); }
    virtual ScopeFilterImplConcrete<Object>*
                    clone() const                       { return new ScopeFilterImplConcrete(*this); }
};

struct ScopeFilter
{
    template<class Object>
                    ScopeFilter(Object o):
                        fo(new ScopeFilterImplConcrete<Object>(o))      {}

                    ScopeFilter(const ScopeFilter& other):
                        fo(other.fo->clone())                           {}
    ScopeFilter&    operator=(const ScopeFilter& other)                 { delete fo; fo = other.fo->clone(); return *this; }
                    ~ScopeFilter()                                      { delete fo; }

    bool            operator()(const Scope& s) const                    { return (*fo)(s); }
    ScopeFilterImpl*    fo;
};

struct FilterRecord
{
    FilterRecord(bool i, ScopeFilter f, severity_level s):
        include(i), filter(f), severity(s)      {}

    bool            include;
    ScopeFilter     filter;
    severity_level  severity;
};
typedef std::vector<FilterRecord>       ScopeFilters;

struct ResourceManager;

inline  Sinks&                                get_sinks();
inline  Streams&                              get_streams();
inline  ResourceManager&                      get_manager();
inline  Scope&                                get_scope();
inline  ScopeFilters&                         get_scope_filters();

static  Sinks&              sinks           = get_sinks();
static  Streams&            streams         = get_streams();
static  ResourceManager&    manager         = get_manager();
static  Scope&              scope           = get_scope();
static  ScopeFilters&       scope_filters   = get_scope_filters();

template<typename C, typename T>
std::basic_ostream<C,T>&
operator<< (std::basic_ostream<C,T>& strm, severity_level lvl)
{
    switch (lvl)
    {
        case trace:
            strm << "trace"; break;
        case debug:
            strm << "debug"; break;
        case info:
            strm << "info"; break;
        case warning:
            strm << "warning"; break;
        case error:
            strm << "error"; break;
        case fatal:
            strm << "fatal"; break;
        default:
            strm << static_cast<int>(lvl); break;
    }

    return strm;
}

// Annotators
struct Annotator
{
        virtual         ~Annotator()    {}
        virtual
        std::ostream&   operator()(std::ostream& out, dlog::severity_level) const =0;
};

struct color_pre:   public Annotator    { virtual inline std::ostream& operator()(std::ostream& out, dlog::severity_level) const; };
struct color_post:  public Annotator    { virtual inline std::ostream& operator()(std::ostream& out, dlog::severity_level) const; };
struct stamp:       public Annotator    { virtual inline std::ostream& operator()(std::ostream& out, dlog::severity_level) const; };
struct level:       public Annotator    { virtual inline std::ostream& operator()(std::ostream& out, dlog::severity_level) const; };
struct emphasis:    public Annotator    { virtual inline std::ostream& operator()(std::ostream& out, dlog::severity_level) const; };
struct flush:       public Annotator    { virtual inline std::ostream& operator()(std::ostream& out, dlog::severity_level) const; };

template<class T>
struct aux_reporter_: public dlog::Annotator
{
                        aux_reporter_(const T& aux_): aux(aux_)              {}
        virtual
        std::ostream&   operator()(std::ostream& out, dlog::severity_level) const
        { return out << " [" << aux << "] "; }

        T               aux;
};

template<class T>
aux_reporter_<T>
aux_reporter(const T& aux)      { return aux_reporter_<T>(aux); }

// A wrapper around std::ostream& that records the functions to be executed
// before and after each message.
class StreamWrapper
{
    public:
        typedef         std::vector<Annotator*>                             AnnotatorVector;

    public:
                        StreamWrapper(std::ostream& out_):
                           out(&out_)                                       {}
                        ~StreamWrapper()                                    { for (unsigned i = 0; i < annot_pre.size(); ++i) delete annot_pre[i]; for (unsigned i = 0; i < annot_post.size(); ++i) delete annot_post[i]; }

        void            append_pre(Annotator* a)                            { annot_pre.push_back(a); }
        void            append_post(Annotator* a)                           { annot_post.push_back(a); }

        void            trigger_pre(severity_level level) const             { for (unsigned i = 0; i < annot_pre.size(); ++i)  (*annot_pre[i])(*out, level); }
        void            trigger_post(severity_level level) const            { for (unsigned i = 0; i < annot_post.size(); ++i) (*annot_post[i])(*out, level); }

        template<class T>
        StreamWrapper&  operator<<(const T& t)                              { (*out) << t; return *this; }

    private:
        std::ostream*       out;
        AnnotatorVector     annot_pre;
        AnnotatorVector     annot_post;
};

// One global instance of this class is responsible for cleaning up the streams.
struct ResourceManager
{ ~ResourceManager()                                                        { for (unsigned i = 0; i < streams.size(); ++i) delete streams[i]; } };

// add_stream(...) returns a temporary instance of StreamAnnotator. Its << and >>
// let us supply functions that get triggered before and after each message.
struct StreamAnnotator
{
                        StreamAnnotator(StreamWrapper& out_):
                            out(&out_)                                          {}

    template<class A>
    StreamAnnotator&    operator<<(A a)                                         { out->append_pre(new A(a)); return *this; }

    template<class A>
    StreamAnnotator&    operator>>(A a)                                         { out->append_post(new A(a)); return *this; }

    StreamWrapper*      out;
};

inline
StreamAnnotator       add_stream(std::ostream& out, severity_level from_lvl);

// LOG_SEV(...) returns a temporary SeverityWrapper instance. Its constructor
// and destructor trigger the pre- and post-message events for each stream.
class SeverityWrapper
{
    public:
                            SeverityWrapper(severity_level lvl_):
                                lvl(lvl_)                                   { for (unsigned i = sinks[lvl]; i < streams.size(); ++i) (*streams[i]).trigger_pre(lvl); }
                            ~SeverityWrapper()                              { for (unsigned i = sinks[lvl]; i < streams.size(); ++i) { (*streams[i]).trigger_post(lvl); (*streams[i]) << '\n'; } }

        template<class T>
        SeverityWrapper&    operator<<(const T& t)                          { for (unsigned i = sinks[lvl]; i < streams.size(); ++i) (*streams[i]) << t; return *this; }

    private:
        severity_level      lvl;
};

class LogScope
{
    public:
#ifndef LOG_NO_SCOPE
                            LogScope(std::string s)                         { scope.push_back(s); }
                            ~LogScope()                                     { scope.pop_back(); }
#else
                            LogScope(std::string)                           {}
                            ~LogScope()                                     {}
#endif
};
inline
void add_scope_filter(ScopeFilter filter, severity_level severity, bool include = true);     // "include = false" means exclude
inline
bool check_scope(severity_level severity);

inline
bool    filter_scope(const dlog::Scope& scope, const std::string& target)
{
    for(unsigned i = 0; i < scope.size(); ++i)
        if (scope[i] == target)
            return true;

    return false;
}

struct FilterScope
{
            FilterScope(const std::string& target_):
                target(target_)                             {}
    bool    operator()(const dlog::Scope& scope) const      { return filter_scope(scope, target); }
    std::string target;
};

#ifdef DEBUG
#ifdef TRACE
static const severity_level min_severity = trace;
#else // !TRACE
static const severity_level min_severity = debug;
#endif
#else // !DEBUG
static const severity_level min_severity = info;
#endif


#define __LOG_AVAILABLE_SINKS(severity) (::dlog::sinks[::dlog::severity] < ::dlog::streams.size())
#define LOG_SEV(severity)               if (::dlog::severity >= ::dlog::min_severity           && __LOG_AVAILABLE_SINKS(severity) && ::dlog::check_scope(::dlog::severity))     ::dlog::SeverityWrapper(::dlog::severity)
#define LOG_SEV_IF(cond, severity)      if (::dlog::severity >= ::dlog::min_severity && (cond) && __LOG_AVAILABLE_SINKS(severity) && ::dlog::check_scope(::dlog::severity))     ::dlog::SeverityWrapper(::dlog::severity)
#define LOG_SCOPE(name)                 ::dlog::LogScope __log_scope__(name)

#ifdef DEBUG
#define AssertMsg(cond, msg)            { if(!(cond)) { LOG_SEV(fatal) << msg; std::abort(); } }
#else
#define AssertMsg(cond, msg)            {}
#endif

#ifdef BREAK_ASSERTS
#define BREAK_IF(cond, msg)             AssertMsg(!cond, msg)
#else
#define BREAK_IF(cond, msg)             if (cond) { LOG_SEV(warning) << msg; std::cin.ignore(); }
#endif

#define UNUSED(x)                       (static_cast<void>(x))
}

dlog::Sinks&
dlog::get_sinks()
{
    static Sinks sinks = { 0, 0, 0, 0, 0, 0 };
    return sinks;
}

dlog::Streams&
dlog::get_streams()
{
    static Streams streams;
    return streams;
}

dlog::ResourceManager&
dlog::get_manager()
{
    static ResourceManager manager;
    return manager;
}

dlog::Scope&
dlog::get_scope()
{
    static Scope scope;
    return scope;
}

dlog::ScopeFilters&
dlog::get_scope_filters()
{
    static ScopeFilters scope_filters;
    return scope_filters;
}

dlog::StreamAnnotator
dlog::add_stream(std::ostream& out, severity_level from_lvl)
{
    StreamWrapper* ps = new StreamWrapper(out);
    if (from_lvl == trace)
        streams.push_back(ps);
    else
        streams.insert(streams.begin() + sinks[from_lvl-1], ps);

    for (unsigned i = 0; i < from_lvl; ++i)
        ++sinks[i];

    return StreamAnnotator(*ps);
}

dlog::severity_level
dlog::severity(const std::string& sev)
{
    if (sev == "trace")    return trace;
    if (sev == "debug")    return debug;
    if (sev == "info")     return info;
    if (sev == "warning")  return warning;
    if (sev == "error")    return error;
    if (sev == "fatal")    return fatal;

    return fatal;
}

// Annotators
std::ostream&
dlog::color_pre::operator()(std::ostream& out, dlog::severity_level level) const
{ if (level > dlog::info && isatty(1)) out << "\033[1;31m"; return out; }

std::ostream&
dlog::color_post::operator()(std::ostream& out, dlog::severity_level level) const
{ if (level > dlog::info && isatty(1)) out << "\033[0m"; return out; }

std::ostream&
dlog::stamp::operator()(std::ostream& out, dlog::severity_level) const
{
    time_t t;
    time(&t);

    struct tm*  timeinfo = localtime(&t);
    char        buffer[80];

    strftime (buffer, 80, "%Y-%m-%d %H:%M:%S", timeinfo);
    return out << buffer << " ";
}

std::ostream&
dlog::level::operator()(std::ostream& out, dlog::severity_level level) const
{ return out << "<" << std::setw(7) << level << "> "; }     // 7 is the length of "warning", the longest label

std::ostream&
dlog::emphasis::operator()(std::ostream& out, dlog::severity_level level) const
{ if (level > dlog::info) out << " --- "; else out << "     "; return out; }

std::ostream&
dlog::flush::operator()(std::ostream& out, dlog::severity_level level) const
{ out.flush(); return out; }

void
dlog::add_scope_filter(ScopeFilter filter, severity_level severity, bool include)
{ scope_filters.push_back(FilterRecord(include, filter, severity)); }

bool
dlog::check_scope(severity_level severity)
{
    for (int i = scope_filters.size() - 1; i >= 0; --i)
    {
        bool            include = scope_filters[i].include;
        ScopeFilter     filter  = scope_filters[i].filter;
        severity_level  sev     = scope_filters[i].severity;

        if (filter(scope))
        {
            if (include && severity >= sev)
                return true;

            if (!include && severity <= sev)
                return false;
        }
    }

    return true;    // default is "include all"
}

#endif  // DLOG_LOG_H
