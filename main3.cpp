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
#include <chrono>

std::tuple<int, int, int, int>
compute_discovered_misses(const std::vector<int> &color_vec, boost::numeric::ublas::matrix<int> &m, int color);

using std::cout;
using std::cerr;
using std::endl;

// k is the number of edges that we want to remove
graph matrix2graph_limited(const boost::numeric::ublas::matrix<int> &m, int kk) {
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
            if (E_missed != 0) edges.emplace_back(i, j, weight);
        }
    }

    sort(begin(edges), end(edges),[&](std::tuple<int, int, int> t1, std::tuple<int, int, int> t2) { return get<2>(t1) > get<2>(t2); });
    for (int i = kk; i < edges.size(); i++) {
        auto[v1, v2, w] = edges[i];
        g.add_edge(v1, v2, w);
    }
//    std::cerr << "num of edges: " << g.num_e() << std::endl;
    return g;
}

std::tuple<int, int, int> get_bounds(int num_colors_natural_full) {
    if (num_colors_natural_full > 220) {
        return {num_colors_natural_full - 200, num_colors_natural_full, 10};
    } else if (num_colors_natural_full > 120) {
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
//    auto matrix_arr = {"mats/nos3.mtx", "mats/plbuckle.mtx", "mats/bcsstk08.mtx", "mats/bcsstk09.mtx", "mats/G51.mtx", "mats/bcsstm13.mtx", "mats/gemat11.mtx"};
    auto matrix_arr = {"mats/plbuckle.mtx"};
    std::ofstream out(std::string("plbuckle_res.csv"));
    out << "numOfColor,nat,ago,lfo,sat,single,multiple,k,ck" << endl;
    auto start = std::chrono::steady_clock::now();
    for (auto matrix_name : matrix_arr) {
        matrix_market mm(matrix_name);
        matrix<int> m = mm.to_ublas_matrix();
        graph g = matrix2graph_limited(m, 0);
        auto[num_colors_natural_full, color_vec_natural_full, after_discovered] = g.greedy_color_martin_idea(g.natural_order(), m, 300);
        auto[from, to, step] = get_bounds(num_colors_natural_full);
//#pragma omp parallel for
        for (int k = 0; k < 1000; k += 100) {
            graph g = matrix2graph_limited(m, k);
            auto[num_colors_natural_full, color_vec_natural_full] = g.greedy_color(1000000);
            auto[num_colors_natural, color_vec_natural] = g.greedy_color(100000);
            std::vector<int> lfo_ord = g.largest_first_order();
            auto[num_colors_lfo, color_vec_lfo] = g.greedy_color_limited(lfo_ord, 100000);
            auto[num_colors_sat, color_vec_sat] = g.saturation_degree_ordering_coloring(100000);
            std::vector<int> ord = g.optimum_order();
            auto[num_colors_newIdea, color_vec_newIdea] = g.greedy_color_limited(ord, 100000);
            for (int color = from; color <= to; color += step) {
                auto[all_sum_nat, all_discovered_color_zero_sum_nat, all_misses_nat, all_misses_color_zero_sum_nat] = compute_discovered_misses(color_vec_natural, m, color);
                auto[all_sum_new, all_discovered_color_zero_sum_new, all_misses_new, all_misses_color_zero_sum_new] = compute_discovered_misses(color_vec_newIdea, m, color);
                auto[all_sum_lfo, all_discovered_color_zero_sum_lfo, all_misses_lfo, all_misses_color_zero_sum_lfo] = compute_discovered_misses(color_vec_lfo, m, color);
                auto[all_sum_sat, all_discovered_color_zero_sum_sat, all_misses_sat, all_misses_color_zero_sum_sat] = compute_discovered_misses(color_vec_sat, m, color);
                int all_discovered_natural = all_sum_nat + all_discovered_color_zero_sum_nat;//compute_discovered(color_vec_natural, m, color);
                int all_discovered_newIdea = all_sum_new + all_discovered_color_zero_sum_new;//compute_discovered(color_vec_newIdea, m, color);
                int all_discovered_lfo = all_sum_lfo + all_discovered_color_zero_sum_lfo;//compute_discovered(color_vec_lfo, m, color);
                int all_discovered_sat = all_sum_sat + all_discovered_color_zero_sum_sat;//compute_discovered(color_vec_sat, m, color);
                auto[num_colors_natural_full, color_vec_natural_full, after_discovered] = g.greedy_color_martin_idea(g.natural_order(), m, color);
                auto[num_colors_multiple, color_vec_multiple] = g.greedy_color_limited(g.natural_order(), color);
                auto[multiple_disc, multiple_zero_disc, multiple_misses, multiple_misses_zero] = compute_discovered_misses(color_vec_multiple, m, color);
                out << color << "," << all_sum_nat << "," << all_sum_new << "," << all_sum_lfo << "," << all_sum_sat << "," <<  after_discovered << "," << multiple_disc << "," << k << "," << num_colors_natural_full<< endl;
            }
        }
    }
    auto end = std::chrono::steady_clock::now();
    cerr << "Elapsed time in milliseconds for the main loop: "
         << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
         << " ms" << endl;
    out.flush();
    out.close();
    return 0;
}

std::tuple<int, int, int, int>
compute_discovered_misses(const std::vector<int> &color_vec, boost::numeric::ublas::matrix<int> &m, int color) {
    auto start = std::chrono::steady_clock::now();
    std::vector<boost::numeric::ublas::vector<int>> discovered(color);
    boost::numeric::ublas::vector<int> zeros = boost::numeric::ublas::zero_vector<int>(m.size2());
    for (int i = 0; i < color; i++) {
        discovered[i] = boost::numeric::ublas::zero_vector<int>(m.size2());
    }
    int cnt = 0;
    for (int i = 0; i < color_vec.size(); i++) {
        if (color_vec[i] < color) {
            discovered[color_vec[i]] += column(m, i);
        } else {
            cnt++;
            zeros += column(m, i);
        }
    }
//    std::cerr << "The number of vertices with color zero: " << cnt << std::endl;
    int all_sum = 0;
    int all_misses = 0;
    for (auto &misse : discovered) {
        for (int j : misse) {
            if (j == 1)
                all_sum += 1;
            else
                all_misses += j;
        }
    }

    int all_discovered_color_zero_sum = 0;
    int all_misses_color_zero_sum = 0;
    for (auto &zero_color_miss : zeros)
        if (zero_color_miss == 1)
            all_discovered_color_zero_sum += 1;
        else
            all_misses_color_zero_sum += zero_color_miss;

    auto end = std::chrono::steady_clock::now();
//    cerr << "Elapsed time in milliseconds : "
//         << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count()
//         << " ms" << endl;
    return {all_sum, all_discovered_color_zero_sum, all_misses, all_misses_color_zero_sum};
}