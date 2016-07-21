#include <iostream>
#include <string>

#include <dlog/stats.h>
#include <dlog/log.h>
#include <opts/opts.h>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>

#include "format.h"

int main(int argc, char** argv)
{
    using namespace boost;
    typedef adjacency_list < vecS, vecS, undirectedS, no_property, property < edge_weight_t, int > > Graph;
    typedef graph_traits < Graph >::edge_descriptor Edge;
    typedef std::pair<int, int> E;

    std::ifstream graph(argv[1]);
    if (graph.is_open())
    {
        int num_nodes;
        graph >> num_nodes;

        std::vector<E> edges;
        std::vector<double> weights;

        int u, v;
        double val;
        int i = 0;
        while (graph >> u >> v >> val)
        {
            Edge e;
            edges.push_back(E(u, v));
            weights.push_back(val);
            i++;
        }

        Graph g(edges.begin(), edges.end(), weights.begin(), num_nodes);
        
        dlog::Timer t;

        std::vector<Edge> mst;
        kruskal_minimum_spanning_tree(g, std::back_inserter(mst));

        fmt::print("Time: {}\n", t.elapsed());

        std::ofstream out(argv[2]);

        for (Edge e : mst) out << source(e, g) << " " << target(e, g) << "\n";

        return 0;
    }
    else return -1;
}
