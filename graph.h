//
// Created by rostam on 14.06.18.
//
#ifndef MY_GCOL_GRAPH_H
#define MY_GCOL_GRAPH_H

#include "prettyprint.hpp"
#include "boost/config.hpp"
#include "boost/graph/adjacency_list.hpp"
#include "boost/cstdlib.hpp"
#include <memory>
#include <iostream>
#include <fstream>
#include <boost/dynamic_bitset.hpp>
#include <tuple>

using boost::property;
using boost::edge_name;
using boost::edge_weight;
using boost::vertex_color_t;
using boost::edge_weight_t;
using boost::edge_name_t;
using boost::property_map;
using boost::graph_traits;
using boost::vertex_color;
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
        property<boost::vertex_color_t, int>,
        property<boost::edge_weight_t, int,
                property<boost::edge_name_t, std::string>>> Graph;

typedef boost::graph_traits<Graph>::vertex_iterator V_iter;
typedef boost::graph_traits<Graph>::edge_iterator E_iter;
typedef Graph::vertex_descriptor Ver;
typedef Graph::edge_descriptor Edge;
typedef boost::property_map<Graph, boost::edge_weight_t>::type edge_weight_type;
typedef boost::dynamic_bitset<> dynbit;

class graph {
    Graph g;
    property_map<Graph, vertex_color_t>::type color;
public:
    /**
     * A constructor which creates a graph with the given number of vertices
     * @param num_vertices the given number of vertices
     */
    graph(unsigned long num_vertices) : g(Graph(num_vertices)) {}

    graph()=default;

    /**
     * Adds an edge with the given vertices as endpoints and its weight
     * @param i
     * @param j
     * @param weight
     */
    void add_edge(int i, int j, int weight) { boost::add_edge(i, j, weight, g); }

    void init_colors() {
        color = get(vertex_color, g);
        for_each_v([&](Ver v) { put_color_v(v, 0); });
    }

    /**
     * Number of vertices
     * @returns number of vertices
     */
    unsigned long num_v() { return boost::num_vertices(g); }

    /**
     * Number of edges
     * @returns number of edges
     */
    unsigned long num_e() { return boost::num_edges(g); }

    /**
     * Returns the color of a given vertex
     * @param v the given matrix
     * @returns the color of the vertex v
     */
    int get_color_v(int v) { return boost::get(vertex_color, g, v); }

    /**
     * Set the color of the given vertex with the given color
     * @param v the given vertex
     * @param c the given color
     */
    void put_color_v(int v, int c) { put(color, v, c); }


    /**
     * for each vertex
     * @tparam Lambda
     * @param g
     * @param func
     */
    template<typename Lambda>
    void for_each_v(Lambda func) {
        V_iter vi, vi_end;
        std::tie(vi, vi_end) = boost::vertices(g);
        std::for_each(vi, vi_end, func);
    }

    /**
     * for each edge
     * @tparam Lambda
     * @param g
     * @param func
     */
    template<typename Lambda>
    void for_each_e(Lambda func) {
        E_iter ei, ei_end;
        tie(ei, ei_end) = edges(g);
        std::for_each(ei, ei_end, func);
    }

    E_iter first_edge() {
        E_iter ei, ei_end;
        tie(ei, ei_end) = edges(g);
        return ei;
    }

    void remove_e(E_iter e) {
        boost::remove_edge(boost::source(*e, g), boost::target(*e, g), g);
    }

    /**
     * for each neighbor of v
     * @tparam Lambda
     * @param g
     * @param v
     * @param func
     */
    template<typename Lambda>
    void for_each_n(const Ver &v, Lambda func) {
        auto av = adjacent_vertices(v, g);
        std::for_each(av.first, av.second, func);
    }

    /**
     * for each distance-2 neighbor of v
     * @tparam Lambda
     * @param g
     * @param v
     * @param func
     */
    template<typename Lambda>
    void for_each_2n(int v, Lambda func) {
        std::set<int> tmp;
        for_each_n(v, [&](Ver a) {
            for_each_n(v, [&](Ver a2) {
                if (v != a2)
                    tmp.insert(a2);
            });
        });

        for (int i : tmp) func(i);
    }

    std::tuple<int, std::vector<int>> greedy_color_limited(const std::vector<int>& order, int max_color) {
        init_colors();
        for (int v : order) {
            std::vector<unsigned int> forbiddenColors(num_v(), -1);
            forbiddenColors[0] = v;
            for_each_n(v, [&](int n) {
                int c = get_color_v(n);
                if (c > 0)forbiddenColors[c] = v;
            });
            //Find first color which can be assigned to v
            auto result = find_if(forbiddenColors.begin(), forbiddenColors.end(), [&](int i) { return i != v; });
            auto res_color = distance(forbiddenColors.begin(), result);
            int c = get_suitable_color(res_color, max_color, v);
            put_color_v(v, c);
        }
        return tuple_numOfColor_Colors();
    }

    std::tuple<int, std::vector<int>> greedy_color_order(const std::vector<int>& order, int max_color) {
        init_colors();
        for (int v : order) {
            // Since vertex numbering starts at zero, we initialize forbiddenColors with the (non-existing) vertex -1
            std::vector<unsigned int> forbiddenColors(num_v(), -1);
            forbiddenColors[0] = v;
            for_each_n(v, [&](int n) {
                int c = get_color_v(n);
                if (c > 0) forbiddenColors[c] = v;
            });
            //Find first color which can be assigned to v
            auto result = find_if(forbiddenColors.begin(), forbiddenColors.end(), [&](int i) { return i != v; });
            auto res_color = distance(forbiddenColors.begin(), result);
            int c = get_suitable_color(res_color, max_color, v);
            put_color_v(v, c);
        }
        return tuple_numOfColor_Colors();
    }

    // In this idea, we will first color greedily and then select between the color groups
    // those with the best discovered numbers.
    // Martin Idea: Max Discovered
    std::tuple<int, std::vector<int>, int> greedy_color_max_discovered(const std::vector<int>& order, const boost::numeric::ublas::matrix<int> m, int max_color) {
        init_colors();
        for (int v : order) {
            std::vector<unsigned int> forbiddenColors(num_v(), -1);
            forbiddenColors[0] = v;
            for_each_n(v, [&](int n) {
                int c = get_color_v(n);
                if (c > 0) forbiddenColors[c] = v;
            });
            auto result = find_if(forbiddenColors.begin(), forbiddenColors.end(), [&](int i) { return i != v; });
            auto res_color = distance(forbiddenColors.begin(), result);
            put_color_v(v, res_color);
        }
        std::vector<int> colors;
        std::set<int> unique_colors;
        for_each_v([&](int v) {
            int the_color = get_color_v(v) - 1;
            colors.push_back(the_color);
            unique_colors.insert(the_color);
        });
        std::vector<boost::numeric::ublas::vector<int>> discovered(unique_colors.size());
        boost::numeric::ublas::vector<int> zeros = boost::numeric::ublas::zero_vector<int>(m.size1());
        for (int i = 0; i < unique_colors.size(); i++) {
            discovered[i] = boost::numeric::ublas::zero_vector<int>(m.size1());
        }
        int cnt = 0;
        for (int i = 0; i < colors.size(); i++) {
            discovered[colors[i]] += column(m, i);
        }
        std::vector<std::pair<int,int>> discover_index_color(unique_colors.size());
        int col = 0;
        for (auto &d : discovered) {
            int all_sum = 0;
            for (int j : d) {
                if (j == 1)
                    all_sum += 1;
            }
            discover_index_color.emplace_back(std::make_pair(col, all_sum));
            col++;
        }

        sort(begin(discover_index_color), end(discover_index_color), [&](std::pair<int,int> t1, std::pair<int,int> t2) { return t1.second > t2.second;});
        std::vector<std::pair<int,int>> max_color_colors;
        int after_discovered = 0;
        int color_to = max_color;
        if(max_color > unique_colors.size()) color_to = unique_colors.size();
        for(int i=0;i<color_to;i++) {
            max_color_colors.emplace_back(discover_index_color[i]);
            after_discovered += discover_index_color[i].second;
        }
        return {unique_colors.size(), colors, after_discovered};
    }

    /**
     * Computes the number of distinct colors in the neighbors of the given vertex
     * @param v The given vertex
     * @return the number of distinct colors in the neighbors of v
     */
    int num_colors_of_neighbors(int v) {
        std::set<int> unique_colors;
        for_each_n(v, [&](int n) {
            int c = get_color_v(n);
            unique_colors.insert(c);
        });
        return unique_colors.size();
    }

    /**
     * Coloring with the consideration that we have a maxnimum number of colors
     * @param max_color the maximum number of colors that can be used for coloring.
     * @return the coloring and number of colors
     */
    std::tuple<int, std::vector<int>> greedy_color(int max_color) {
        init_colors();
        std::vector<int> order(num_v(), 0);
        for_each_v([&](int v) {
            order[v] = v;
        });
        return greedy_color_order(order, max_color);
    }

    /*
     * Computes the natural ordering which is nothing than a list of numbers from zero
     * the number of vertices minus 1
     * @return the natural ordering
     */
    std::vector<int> natural_order() {
        init_colors();
        std::vector<int> order(num_v(), 0);
        for_each_v([&](int v) {
            order[v] = v;
        });
        return order;
    }

    /**
     *
     * @return
     */
    std::vector<int> optimum_order() {
        std::vector<std::tuple<int, int>> vertex_weight;
        for_each_v([&](int v) {
            int sum = 0;
            int max = -1000;
            for_each_n(v, [&](int n) {
                int w = get(boost::edge_weight_t(), g, edge(v, n, g).first);
                if (w > max)
                    max = w;
                sum += w;
            });
            vertex_weight.emplace_back(v, sum);
        });

        std::sort(begin(vertex_weight), end(vertex_weight), [](auto const &t1, auto const &t2) {
            return get<1>(t1) > get<1>(t2);
        });

        std::vector<int> ret;
        for (auto &t : vertex_weight)
            ret.push_back(get<0>(t));
        return ret;
    }

    /**
     *
     * @return
     */
    std::vector<int> largest_first_order() {
        std::vector<std::tuple<int, int>> vertex_weight;
        for_each_v([&](int v) {
            vertex_weight.emplace_back(v, degree(v, g));
        });

        std::sort(begin(vertex_weight), end(vertex_weight), [](auto const &t1, auto const &t2) {
            return get<1>(t1) > get<1>(t2); // or use a custom compare function
        });

        std::vector<int> ret;
        for (std::tuple<int, int> t : vertex_weight)
            ret.push_back(get<0>(t));
        return ret;
    }

    /**
     *
     * @param max_color
     * @return
     */
    std::tuple<int, std::vector<int>> saturation_degree_ordering_coloring(int max_color) {
        init_colors();
        V_iter vi, vi_end;
        std::tie(vi, vi_end) = vertices(g);
        std::list<int> vs(vi, vi_end);
        while (!vs.empty()) {
            int sat_v;
            std::vector<std::tuple<int, int>> vertex_num_color_neighbors;
            for (int v : vs) {
                vertex_num_color_neighbors.emplace_back(v, num_colors_of_neighbors(v));
            }

            std::tuple<int, int> vc =
                    *std::max_element(std::begin(vertex_num_color_neighbors), std::end(vertex_num_color_neighbors),
                                      [&](std::tuple<int, int> vc1, std::tuple<int, int> vc2) {
                                          return get<1>(vc1) < get<1>(vc2);
                                      });
            sat_v = std::get<0>(vc);

            std::vector<unsigned int> forbiddenColors(num_v(), -1);
            forbiddenColors[0] = sat_v;
            for_each_n(sat_v, [&](int n) {
                int c = get_color_v(n);
                if (c > 0)forbiddenColors[c] = sat_v;
            });
            //Find first color which can be assigned to v
            auto result = find_if(forbiddenColors.begin(), forbiddenColors.end(), [&](int i) { return i != sat_v; });
            auto res_color = distance(forbiddenColors.begin(), result);
            int c = get_suitable_color(res_color, max_color, sat_v);
            put_color_v(sat_v, c);
            vs.remove(sat_v);
        }
        return tuple_numOfColor_Colors();
    }

    /**
     *
     * @return
     */
    std::tuple<int, std::vector<int>> tuple_numOfColor_Colors() {
        std::vector<int> colors;
        std::set<int> unique_colors;
        for_each_v([&](int v) {
            colors.push_back(get_color_v(v));
            int the_color = get_color_v(v) - 1;
            if (the_color != -1)
                unique_colors.insert(the_color);
            else
                unique_colors.insert(0);
        });
        return {unique_colors.size(), colors};
    }

    /**
     *
     * @param res_color
     * @param max_color
     * @param v
     * @return
     */
    int get_suitable_color(int res_color, int max_color, int v) {
        if (res_color < max_color) {
            return res_color;
        } else {
            typename graph_traits<Graph>::out_edge_iterator ei, ei_end;
            int nv = 0;
            double max_w = -1000;
            for (boost::tie(ei, ei_end) = out_edges(v, g); ei != ei_end; ++ei) {
                double w = get(boost::edge_weight_t(), g, *ei);
                auto source = boost::source(*ei, g);
                auto target = boost::target(*ei, g);
                if (w > max_w) {
                    if (boost::get(vertex_color, g, nv) != -1) {
                        max_w = w;
                        nv = target;
                    }
                }
            }
//            if(max_w < -1) return 0;
//            std::cerr<<max_w<<std::endl;
            if (boost::get(vertex_color, g, nv) != -1) {
                return boost::get(vertex_color, g, nv);
            } else {
                return 1;
            }
        }
    }

//
//    std::pair<int,int> d2color() {
////        std::vector<unsigned int> V = V_c;
//        property_map<Graph, vertex_color_t>::type color = get(vertex_color, G_b);
//        std::vector<unsigned int> N_2;
//        vector<unsigned int> forbiddenColors(num_vertices(G_b), -1);
//        //All edges in E_S have edge_weight=1; otherwise edge_weight=0
//        //Initialize colors
//        for_each(V.begin(), V.end(), [&](Ver v) { put(color, v, 0); });
//        //Iterate over all vertices which should be colored
//        for_each(V.begin(), V.end(), [&](unsigned int v) {
//            forbiddenColors[0] = v;
//            if (IncidentToReqEdge(G_b, v)) {
//                //Get the distance-2 neighbors
//                if (restricted)
//                    N_2 = neighbors::N_2restricted(G_b, v);
//                else
//                    N_2 = neighbors::N_2(G_b, v);
//                //Iterate over distance-2 neighbors
//                for_each(N_2.begin(), N_2.end(), [&](unsigned int n_2) {
//                    //Mark colors which are used by distance-2 neighbors in forbiddenColors
//                    if (get(vertex_color, G_b, n_2) > 0) {
//                        forbiddenColors[get(vertex_color, G_b, n_2)] = v;
//                    }
//                });
//
//                //Find first color which can be assigned to v
//                vector<unsigned int>::iterator result = find_if(forbiddenColors.begin(), forbiddenColors.end(),
//                                                                bind1st(not_equal_to<int>(), v));
//
//                //Color v
//                put(color, v, distance(forbiddenColors.begin(), result));
//            } else {
//                put(color, v, 0);
//            }
//        });
//        return num_colors_d2(G_b);
//    }

};




#endif //MY_GCOL_GRAPH_H
