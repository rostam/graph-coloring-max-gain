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
compute_misses(int num_colors, const std::vector<int> &color_vec, boost::numeric::ublas::matrix<int> &m);

using std::cout;
using std::endl;

graph matrix2graph_limited(const boost::numeric::ublas::matrix<int> &m, int index) {
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
//            if (weight > index) {
//                g.add_edge(i, j, weight);
//            }
//            if(weight < index) g.add_edge(i, j, weight);
        }
    }
    sort(begin(edges), end(edges), [&](std::tuple<int, int, int> t1, std::tuple<int, int, int> t2) {
        return get<2>(t1) > get<2>(t2);
    });

    for (int i = index; i < edges.size(); i++) {
        auto[v1, v2, w] = edges[i];
        g.add_edge(v1, v2, w);
    }
    return g;
}

graph matrix2graph(const boost::numeric::ublas::matrix<int> &m) {
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

    for (int i = 0; i < edges.size(); i++) {
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
//    auto matrix_arr = {"nos3.mtx", "plbuckle.mtx", "bcsstk08.mtx", "1138_bus.mtx", "G51.mtx", "bcsstm13.mtx"};
    auto matrix_arr = {"gemat11.mtx"};
    for(auto matrix_name : matrix_arr) {
        std::cout << matrix_name << " " << std::endl;
        matrix_market mm(matrix_name);
        matrix<int> m = mm.to_ublas_matrix();
//    omp_set_num_threads(4);
//#pragma omp parallel for
        graph g = matrix2graph(m);
        auto[num_colors_natural_full, color_vec_natural_full] = g.greedy_color(1000000);
        cout << "num of colors of full coloring: " << num_colors_natural_full << endl;
        if (num_colors_natural_full <= 10) return 10;
        auto [from, to, step] = get_bounds(num_colors_natural_full);
        std::cout << from << " " << to << " " << step << std::endl;
        std::ofstream out(std::string(matrix_name) + ".csv");
        out << "numOfColor,mnat,mnew,mlfo,msat" << endl;
        for (int color = from; color <= to; color += step) {
            std::cerr << color << endl;
            auto[num_colors_natural, color_vec_natural] = g.greedy_color(color);
            int all_misses_natural = compute_misses(num_colors_natural, color_vec_natural, m);

            std::vector<int> ord = g.optimum_order();
            auto[num_colors_newIdea, color_vec_newIdea] = g.greedy_color_limited(ord, color);
            int all_misses_newIdea = compute_misses(num_colors_newIdea, color_vec_newIdea, m);

            std::vector<int> lfo_ord = g.largest_first_order();
            auto[num_colors_lfo, color_vec_lfo] = g.greedy_color_limited(lfo_ord, color);
            int all_misses_lfo = compute_misses(num_colors_lfo, color_vec_lfo, m);

            auto[num_colors_sat, color_vec_sat] = g.saturation_degree_ordering_coloring(color);
            int all_misses_sat = compute_misses(num_colors_sat, color_vec_sat, m);

//#pragma omp critical
            out << color << ","
                << all_misses_natural << "," << all_misses_newIdea << ","
                << all_misses_lfo << "," << all_misses_sat << endl;

//        }
        }
        out.flush();
        out.close();
    }
    return 0;
}

int compute_misses(int num_colors, const std::vector<int> &color_vec, boost::numeric::ublas::matrix<int> &m) {
    std::vector<boost::numeric::ublas::vector<int>> misses(num_colors);
    for (int i = 0; i < num_colors; i++) {
        misses[i] = boost::numeric::ublas::zero_vector<int>(m.size2());
    }
    for (int i = 0; i < color_vec.size(); i++) {
        if(color_vec[i] >= misses.size()) continue;
        misses[color_vec[i]] += column(m, i);
    }
//    std::cerr << misses << std::endl;
    int all_sum = 0;
    for (auto &misse : misses) {
        for (int j : misse) {
            if (j == 1)
                all_sum += j;
//        misses[i] = 2;
//        all_sum += sum(misses[i]);
        }
    }
    return all_sum;
}