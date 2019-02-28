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

graph matrix2graph(const boost::numeric::ublas::matrix<int> &m, int index) {
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

int main(int argc, const char *argv[]) {
    using boost::numeric::ublas::matrix;
    using boost::numeric::ublas::matrix_column;
    auto matrix_name = argv[1];
    std::cout << matrix_name << " ";
//    int min_possible_color = std::stoi(argv[2]);
//    int max_possible_color = std::stoi(argv[3]);
//    matrix_market mm("/home/rostam/b1_ss.mtx");
    matrix_market mm(matrix_name);
//    matrix_market mm("/home/rostam/nos3.mtx");
    matrix<int> m = mm.to_ublas_matrix();
    int min_index = std::stoi(argv[2]);
    int max_index = std::stoi(argv[3]);
    std::cout << min_index << " " << max_index << " ";
    int color = 100000;

    std::ofstream out("results_" + std::string(matrix_name) + "_" + std::string(argv[2]) + "_" + std::string(argv[3]) + ".csv");
    out << "num_edges,cnat,cnew,clfo,csat,mnat,mnew,mlfo,msat" << endl;
//    omp_set_num_threads(4);
//#pragma omp parallel for
    for (int index = min_index; index <= max_index; index+=1) {
        cout << index << endl;
        graph g = matrix2graph(m, index);
//        typedef property_map<Graph, boost::vertex_index_t>::const_type vertex_index_map;
//        boost::iterator_property_map<int *, vertex_index_map> color(&color_vec.front(), boost::get(boost::vertex_index, g));
//        int num_colors = boost::sequential_vertex_coloring(g, color);

//        out << "threshold_for_edge_weights " << index << endl;
//        for (int color = min_possible_color; color < max_possible_color; color++) {
        auto [num_colors_natural, color_vec_natural] = g.greedy_color(color);
        int all_misses_natural = compute_misses(num_colors_natural, color_vec_natural, m);

        std::vector<int> ord = g.optimum_order();
        auto[num_colors_newIdea, color_vec_newIdea] = g.greedy_color_order(ord, color);
        int all_misses_newIdea = compute_misses(num_colors_newIdea, color_vec_newIdea, m);

        std::vector<int> lfo_ord = g.largest_first_order();
        auto[num_colors_lfo, color_vec_lfo] = g.greedy_color_order(lfo_ord, color);
        int all_misses_lfo = compute_misses(num_colors_lfo, color_vec_lfo, m);

        auto[num_colors_sat, color_vec_sat] = g.saturation_degree_ordering_coloring(color);
        int all_misses_sat = compute_misses(num_colors_sat, color_vec_sat, m);

//#pragma omp critical
        out << index << "," << num_colors_natural << "," << num_colors_newIdea << ","
            << num_colors_lfo << "," << num_colors_sat << ","
            << all_misses_natural << "," << all_misses_newIdea << ","
            << all_misses_lfo << "," << all_misses_sat << endl;

//        }
    }
    out.flush();
    out.close();
    return 0;
}

int compute_misses(int num_colors, const std::vector<int> &color_vec, boost::numeric::ublas::matrix<int> &m) {
    int all_sum = 0;
    std::vector<boost::numeric::ublas::vector<int>> misses(num_colors);
    for (int i = 0; i < num_colors; i++) {
        misses[i] = boost::numeric::ublas::zero_vector<int>(m.size2());
    }
    for (int i = 0; i < color_vec.size(); i++) {
        misses[color_vec[i]] += column(m, i);
    }

    for (int i = 0; i < misses.size(); i++) {
        for (int j = 0; j < misses[i].size(); j++) {
            if (misses[i][j] > 1) all_sum+=misses[i][j];
//        misses[i] = 2;
//        all_sum += sum(misses[i]);
        }
    }

    return all_sum;
}