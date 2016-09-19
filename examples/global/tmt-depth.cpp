#include <iostream>

#include <opts/opts.h>

#include <reeber/grid.h>
#include <reeber/triplet-merge-tree.h>
#include <reeber/triplet-merge-tree-serialization.h>
namespace r = reeber;

typedef     REEBER_REAL                       Real;
typedef     r::Grid<Real, 3>                  Grid;
typedef     Grid::Index                       Index;
typedef     Grid::Value                       Value;
typedef     r::TripletMergeTree<Index, Value> TripletMergeTree;
typedef     TripletMergeTree::Neighbor        Neighbor;

int main(int argc, char** argv)
{
    using namespace opts;
    Options ops(argc, argv);

    std::string infn;
    if (  ops >> Present('h', "help", "show help message") ||
        !(ops >> PosOption(infn)))
    {
        fmt::print("Usage: {} IN.tmt\n{}", argv[0], ops);
        return 1;
    }

    diy::MemoryBuffer bb;
    bb.read(infn);

    TripletMergeTree mt;
    diy::load(bb, mt);

    std::unordered_map<Neighbor,size_t>  depth;
    for (auto& x : static_cast<const TripletMergeTree&>(mt).nodes())
    {
        Neighbor n = x.second;

        std::vector<Neighbor> path { n };
        Neighbor next = std::get<1>(n->parent());
        while (n != next)
        {
            auto it = depth.find(n);
            if (it != depth.end())
            {
                size_t cur_depth = it->second;
                for (size_t i = 0; i < path.size(); ++i)
                {
                    Neighbor y = path[i];
                    depth[y] = cur_depth + (path.size() - i);
                }
                path.clear();
                break;
            }

            n = next;
            path.push_back(n);
            next = std::get<1>(n->parent());
        }
        for (size_t i = 0; i < path.size(); ++i)
        {
            Neighbor y = path[i];
            depth[y] = (path.size() - i);
        }
    }

    size_t max_depth = 0;
    for (auto& x : depth)
        if (x.second > max_depth)
            max_depth = x.second;
    fmt::print("Max depth: {}\n", max_depth);
}
