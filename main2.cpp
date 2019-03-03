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

int
compute_misses(int num_colors, const std::vector<int> &color_vec, boost::numeric::ublas::matrix<int> &m, int color);

using std::cout;
using std::endl;

// k is the number of edges that we want to remove
graph matrix2graph_limited(const boost::numeric::ublas::matrix<int> &m, int k) {
    graph g(m.size1());
    std::vector<std::tuple<int, int, int>> edges;
    for (int i = 0; i < m.size1(); i++) {
        for (int j = i + 1; j < m.size1(); j++) {
            int E_discovered = 0, E_missed = 0;
            for (int k = 0; k < m.size2(); k++) {
                if (m(k, i) != 0 && m(k, j) == 0) {
                    E_discovered++;
                } else if (m(k, i) == 0 && m(k, j) != 0) {
                    E_discovered++;
                } else if (m(k, i) != 0 && m(k, j) != 0) {
                    E_missed += 2;
                }
            }
            int weight = E_discovered - E_missed;
            if (E_missed != 0) edges.push_back({i, j, weight});
        }
    }
    sort(begin(edges), end(edges), [&](std::tuple<int, int, int> t1, std::tuple<int, int, int> t2) {
        return get<2>(t1) > get<2>(t2);
    });

    for (int i = k; i < edges.size(); i++) {
        auto[v1, v2, w] = edges[i];
        g.add_edge(v1, v2, w);
    }
    return g;
}

std::tuple<int, int, int> get_bounds(int num_colors_natural_full) {
    if(num_colors_natural_full > 220) {
        return {num_colors_natural_full - 200, num_colors_natural_full, 20};
    } else if(num_colors_natural_full > 120) {
        return {num_colors_natural_full - 100, num_colors_natural_full, 10};
    } else if (num_colors_natural_full > 35) {
        return {num_colors_natural_full - 30, num_colors_natural_full, 3};
    } else if (num_colors_natural_full > 25) {
        return {num_colors_natural_full - 20, num_colors_natural_full, 2};
    } else {
        return {num_colors_natural_full - 10, num_colors_natural_full, 1};
    }
}

int main(int argc, const char *argv[]) {
    using boost::numeric::ublas::matrix;
    using boost::numeric::ublas::matrix_column;
//    auto matrix_arr = {"nos3.mtx", "plbuckle.mtx", "bcsstk08.mtx", "1138_bus.mtx", "G51.mtx", "bcsstm13.mtx", "gemat11.mtx"};
    auto matrix_arr = {"G51.mtx"};
    int sample = 75;
    std::ofstream out2(std::string("G51.mtx") + "k.csv");
    out2 << "k,mnat,mnew,mlfo,msat" << endl;
    for (int k = 0; k < 10; k+=100) {
        for (auto matrix_name : matrix_arr) {
            std::cout << matrix_name << " " << std::endl;
            matrix_market mm(matrix_name);
            matrix<int> m = mm.to_ublas_matrix();
            graph g = matrix2graph_limited(m, k);
            auto[num_colors_natural_full, color_vec_natural_full] = g.greedy_color(1000000);
            cout << "num of colors of full coloring: " << num_colors_natural_full << endl;
            if (num_colors_natural_full <= 10) return 10;
            auto[from, to, step] = get_bounds(num_colors_natural_full);
            std::cout << from << " " << to << " " << step << std::endl;
            std::ofstream out(std::string(matrix_name) + std::to_string(k) + ".csv");
            out << "numOfColor,mnat,mnew,mlfo,msat" << endl;
            auto[num_colors_natural, color_vec_natural] = g.greedy_color(100000);
            std::vector<int> ord = g.optimum_order();
            auto[num_colors_newIdea, color_vec_newIdea] = g.greedy_color_limited(ord, 100000);
            std::vector<int> lfo_ord = g.largest_first_order();
            auto[num_colors_lfo, color_vec_lfo] = g.greedy_color_limited(lfo_ord, 100000);
            auto[num_colors_sat, color_vec_sat] = g.saturation_degree_ordering_coloring(100000);
            for (int color = 75; color <= to; color += step) {
                std::cerr << color << endl;
                int all_misses_natural = compute_misses(num_colors_natural, color_vec_natural, m, color);
                int all_misses_newIdea = compute_misses(num_colors_newIdea, color_vec_newIdea, m, color);
                int all_misses_lfo = compute_misses(num_colors_lfo, color_vec_lfo, m, color);
                int all_misses_sat = compute_misses(num_colors_sat, color_vec_sat, m, color);
                out << color << "," << all_misses_natural << "," << all_misses_newIdea << ","
                    << all_misses_lfo << "," << all_misses_sat << endl;
                if(color == sample) {
                    out2 << k << "," << all_misses_natural << "," << all_misses_newIdea << ","
                         << all_misses_lfo << "," << all_misses_sat << endl;
                    break;
                }
            }
            out.flush();
            out.close();
        }
    }
    out2.flush();
    out2.close();
    return 0;
}

int compute_misses(int num_colors, const std::vector<int> &color_vec, boost::numeric::ublas::matrix<int> &m, int color) {
    std::vector<boost::numeric::ublas::vector<int>> misses(color);
//    boost::numeric::ublas::vector<int> zero_color = boost::numeric::ublas::zero_vector<int>(m.size2());
    for (int i = 0; i < color; i++) {
        misses[i] = boost::numeric::ublas::zero_vector<int>(m.size2());
    }
    int cnt = 0;
    for (int i = 0; i < color_vec.size(); i++) {
//        if(color_vec[i] >= color) {
//            zero_color += column(m,i);
//        } else {
//            misses[color_vec[i]] += column(m, i);
//        }
        if(color_vec[i] < color) {
            misses[color_vec[i]] += column(m, i);
        } else {
            cnt++;
        }
    }
    std::cerr << "The numbber of vertices with color zero: " <<  cnt << std::endl;
    int all_sum = 0;
    for (auto &misse : misses) {
        for (int j : misse) {
            if (j == 1)
                all_sum += j;
        }
    }
    return all_sum;
}