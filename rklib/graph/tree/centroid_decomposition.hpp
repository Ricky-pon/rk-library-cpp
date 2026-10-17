#pragma once
// https://github.com/Ricky-pon/rk-library-cpp/blob/main/rklib/graph/tree/centroid_decomposition.hpp

#include <algorithm>
#include <cassert>
#include <utility>
#include <vector>

#include "rklib/graph/graph.hpp"

namespace rklib {

// 重心分解
// 重心木の構造だけを持ち、問題ごとのデータは利用者が重心を添字とする配列で持つ
template<class T = int>
struct CentroidDecomposition {
public:
    struct Level {
        int cent;  // この段の重心
        T dist;    // dist(v, c)
        int idx;   // vertices(c) における v の位置
    };

    // key: vertices(c) の並び順（省略時は頂点番号順）
    //   例：行きがけ順の番号を渡すと、各部分木が連続した区間になる
    CentroidDecomposition(const Graph<T>& gr, std::vector<int> key = {})
        : n_(gr.size()),
          key_(std::move(key)),
          par_(n_, -1),
          vs_(n_),
          levels_(n_)
    {
        if (key_.empty()) {
            key_.resize(n_);
            for (int v = 0; v < n_; ++v) {
                key_[v] = v;
            }
        }
        assert((int)key_.size() == n_);
        if (n_ == 0) {
            return;
        }
        assert(gr.edge_size() == n_ - 1);
        auto by_key = [&](int u, int v) { return key_[u] < key_[v]; };

        std::vector<bool> removed(n_, false);
        std::vector<int> sz(n_), pv(n_), order;
        std::vector<T> dist(n_);
        // (成分内の頂点, 親の重心)
        std::vector<std::pair<int, int>> stk = {{0, -1}};

        auto bfs = [&](int s) {
            order.assign(1, s);
            pv[s] = -1;
            dist[s] = 0;
            for (int i = 0; i < (int)order.size(); ++i) {
                int v = order[i];
                for (auto& e : gr[v]) {
                    if (e.to == pv[v] || removed[e.to]) {
                        continue;
                    }
                    pv[e.to] = v;
                    dist[e.to] = dist[v] + e.cost;
                    order.push_back(e.to);
                }
            }
        };

        while (!stk.empty()) {
            auto [s, p] = stk.back();
            stk.pop_back();

            bfs(s);
            int total = order.size();
            for (int i = total - 1; i >= 0; --i) {
                int v = order[i];
                sz[v] = 1;
                for (auto& e : gr[v]) {
                    if (e.to != pv[v] && !removed[e.to]) {
                        sz[v] += sz[e.to];
                    }
                }
            }
            int c = s;
            for (bool moved = true; moved;) {
                moved = false;
                for (auto& e : gr[c]) {
                    if (e.to != pv[c] && !removed[e.to] &&
                        sz[e.to] > total / 2) {
                        c = e.to;
                        moved = true;
                        break;
                    }
                }
            }

            bfs(c);
            par_[c] = p;
            vs_[c] = order;
            std::sort(vs_[c].begin(), vs_[c].end(), by_key);
            for (int i = 0; i < total; ++i) {
                int v = vs_[c][i];
                levels_[v].push_back({c, dist[v], i});
            }

            removed[c] = true;
            for (auto& e : gr[c]) {
                if (!removed[e.to]) {
                    stk.emplace_back(e.to, c);
                }
            }
        }
    }

    int size() const { return n_; }

    // 重心木における c の親（c が根なら -1）
    int parent(int c) const
    {
        assert(0 <= c && c < n_);
        return par_[c];
    }

    // 重心が c である成分の頂点（key の昇順）
    const std::vector<int>& vertices(int c) const
    {
        assert(0 <= c && c < n_);
        return vs_[c];
    }

    // v を含む成分の重心を、根側から順に並べたもの
    // （最後は v 自身）
    const std::vector<Level>& levels(int v) const
    {
        assert(0 <= v && v < n_);
        return levels_[v];
    }

    // vertices(c) のうち key が [l, r) に入るものの添字の範囲
    std::pair<int, int> index_range(int c, int l, int r) const
    {
        assert(0 <= c && c < n_);
        auto& a = vs_[c];
        auto lt = [&](int v, int x) { return key_[v] < x; };
        int li = std::lower_bound(a.begin(), a.end(), l, lt) - a.begin();
        int ri = std::lower_bound(a.begin(), a.end(), r, lt) - a.begin();
        return {li, std::max(li, ri)};
    }

private:
    int n_;
    std::vector<int> key_;
    std::vector<int> par_;
    std::vector<std::vector<int>> vs_;
    std::vector<std::vector<Level>> levels_;
};

}  // namespace rklib
