#include <iostream>
#include <numeric>
#include <string>
#include <algorithm>
#include <unordered_set>

#include <dlog/stats.h>
#include <dlog/log.h>
#include <opts/opts.h>

#include <reeber/triplet-merge-tree.h>
namespace r = reeber;

class Topology
{
    public:
                Topology(int size, std::unordered_map<int, std::vector<int>>&& link):
                    size_(size), link_(std::move(link))    {}

        int     size() const                    { return size_; }
        std::vector<int>
                vertices() const                { std::vector<int> v(size_); std::iota(v.begin(), v.end(), 0); return v; }
        const std::vector<int>&
                link(int u) const               { auto it = link_.find(u); if (it == link_.end()) return vector_; else return it->second; }

    private:
        int                                         size_;
        std::unordered_map<int, std::vector<int>>   link_;
        std::vector<int>                            vector_;
};

int main(int argc, char** argv)
{
    using namespace opts;
    int threads;
    Options ops(argc, argv);
    ops >> Option('t', "threads",   threads,      "number of threads to use (with TBB)");
    r::task_scheduler_init init(threads);

    std::ifstream graph(argv[1]);
    if (graph.is_open())
    {
        int n;
        graph >> n;
        int count = n;

        std::unordered_map<int, float>                  f;
        std::unordered_map<int, std::vector<int>>       link;
        std::unordered_map<int, std::tuple<int, int>>   edges;

        int u, v;
        float val;
        while (graph >> u >> v >> val)
        {
            edges[count] = std::make_tuple(u, v);

            if (link.find(u) != link.end()) link[u].push_back(count);
            else link[u] = {count};

            if (link.find(v) != link.end()) link[v].push_back(count);
            else link[v] = {count};

            link[count] = {u, v};

            f[count] = val;

            count++;
        }

        dlog::Timer t;

        r::TripletMergeTree<int, float> mt;
        Topology topology(count, std::move(link));
        r::compute_merge_tree(mt, topology, [&f](int u) { return f[u]; });

        fmt::print("Time: {}\n", t.elapsed());

        std::unordered_set<int> mst;
        const auto& mt_ = mt;
        for (auto& kv : mt_.nodes())
        {
            auto s = std::get<0>(kv.second->parent());
            if (s != kv.second) mst.insert(s->vertex);
        }
        
        std::ofstream out(argv[2]);

        for (int i : mst)
        {
            int u, v;
            std::tie(u, v) = edges[i];
            out << u << " " << v << "\n";
        }

        return 0;
    }
    else return -1;
}
