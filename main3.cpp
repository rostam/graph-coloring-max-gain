#include <iostream>
#include "Mtx2Graph.hpp"
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include "graph.h"
#include <boost/graph/sequential_vertex_coloring.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <fstream>
#include <numeric>
#include <omp.h>

using std::cout;
using std::endl;

int compute_misses(int num_colors, const std::vector<int> &color_vec, const graph& g, int color) {

}

int main(int argc, const char *argv[]) {
    using boost::numeric::ublas::matrix;
    matrix_market mm("mats/bcsstk08.mtx");
    std::map<int, std::vector<int>> my = mm.to_mymat();
    graph g;
    for (const auto &p : my) {
        const std::vector<int> &v = p.second;
        for(int i=0;i<v.size();i++) {
            for(int j=i+1;j<v.size();j++) {
                g.add_edge(i,j,0);
            }
        }
    }

    cout << g.num_e() << endl;
    for(int i=0;i < 50;i++) {
        g.remove_e(g.first_edge());
        auto[num_colors_natural, color_vec_natural] = g.greedy_color(1000000);
        std::vector<int> ord = g.optimum_order();
        auto[num_colors_newIdea, color_vec_newIdea] = g.greedy_color_limited(ord, 100000);
        std::vector<int> lfo_ord = g.largest_first_order();
        auto[num_colors_lfo, color_vec_lfo] = g.greedy_color_limited(lfo_ord, 100000);
        auto[num_colors_sat, color_vec_sat] = g.saturation_degree_ordering_coloring(100000);
        for (int color = 0; color <= 10; color += 1) {
            std::cout << color << endl;
            int all_misses_natural = compute_misses(num_colors_natural, color_vec_natural, g, color);
            int all_misses_newIdea = compute_misses(num_colors_newIdea, color_vec_newIdea, g, color);
            int all_misses_lfo = compute_misses(num_colors_lfo, color_vec_lfo, g, color);
            int all_misses_sat = compute_misses(num_colors_sat, color_vec_sat, g, color);
            cout << color << "," << all_misses_natural << "," << all_misses_newIdea << ","
                << all_misses_lfo << "," << all_misses_sat << endl;
//            if(color == sample) {
//                out2 << k << "," << all_misses_natural << "," << all_misses_newIdea << ","
//                     << all_misses_lfo << "," << all_misses_sat << endl;
//                break;
//            }
        }
//        cout << "num of colors of full coloring: " << num_colors_natural_full << endl;
//        cout << i << endl;
    }


    return 0;
}