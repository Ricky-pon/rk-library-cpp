#pragma once
// https://github.com/Ricky-pon/rk-library-cpp/blob/main/rklib/graph/graph.hpp

#include <cassert>
#include <vector>

namespace rklib {

template<class T = int>
struct Edge {
    int from, to;
    T cost;
    int index;

    Edge() = default;
    Edge(const int from, const int to, const T cost = 1, const int index = 0)
        : from(from),
          to(to),
          cost(cost),
          index(index)
    {
    }

    // v でない方の端点
    int opposite(int v) const
    {
        assert(from == v || to == v);
        return from == v ? to : from;
    }
};

template<class T = int>
struct Graph {
    explicit Graph(int n = 0) : adj_(n) {}

    // 有向辺 from->to を追加し、id を返す
    int add_edge(const int from, const int to, const T cost = 1)
    {
        assert(0 <= from && from < size());
        assert(0 <= to && to < size());
        const int id = edge_list_.size();
        edge_list_.emplace_back(from, to, cost, id);
        adj_[from].emplace_back(from, to, cost, id);
        return id;
    }

    // 無向辺 u-v を追加し、id を返す
    int add_undirected_edge(const int u, const int v, const T cost = 1)
    {
        const int id = add_edge(u, v, cost);
        adj_[v].emplace_back(v, u, cost, id);
        return id;
    }

    // 頂点数
    int size() const { return adj_.size(); }

    // 辺数
    int edge_size() const { return edge_list_.size(); }

    // v から出る辺のリスト
    const std::vector<Edge<T>>& operator[](int v) const
    {
        assert(0 <= v && v < size());
        return adj_[v];
    }

    // id に対応する辺を返す。無向グラフの場合、向きは追加時の引数順
    const Edge<T>& edge(int id) const
    {
        assert(0 <= id && id < edge_size());
        return edge_list_[id];
    }

    // 辺のリスト（id 順）
    const std::vector<Edge<T>>& edges() const { return edge_list_; }

private:
    std::vector<std::vector<Edge<T>>> adj_;
    std::vector<Edge<T>> edge_list_;
};

}  // namespace rklib
