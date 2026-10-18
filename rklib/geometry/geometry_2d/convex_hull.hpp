#pragma once
// https://github.com/Ricky-pon/rk-library-cpp/blob/main/rklib/geometry/geometry_2d/convex_hull.hpp

#include <algorithm>
#include <cassert>
#include <numeric>
#include <utility>
#include <vector>

#include "rklib/geometry/geometry_2d/base.hpp"

namespace rklib {

namespace geo2d {

namespace internal {

// convex_hull と同じ順番で、点の添字を返す。重複する点は、添字が最も小さいものを残す
template<class T, class W>
std::vector<int> convex_hull_index(const std::vector<BasicPoint<T, W>>& p,
                                   bool keep_collinear = false)
{
    const int n = p.size();
    std::vector<int> ord(n);
    std::iota(ord.begin(), ord.end(), 0);
    std::stable_sort(ord.begin(), ord.end(),
                     [&](int i, int j) { return lex_less(p[i], p[j]); });
    ord.erase(std::unique(ord.begin(), ord.end(),
                          [&](int i, int j) { return eq(p[i], p[j]); }),
              ord.end());
    const int m = ord.size();
    if (m <= 2) {
        return ord;
    }

    // p[a] -> p[b] -> p[c] と進むとき、p[b] を取り除くか。
    // 右に曲がるときは常に取り除き、まっすぐ進むときは keep_collinear でなければ取り除く
    auto bad = [&](int a, int b, int c) {
        int s = sgn(cross(p[b] - p[a], p[c] - p[b]));
        return keep_collinear ? s < 0 : s <= 0;
    };

    // 上側は ord[upper_start] から逆順にたどる。ord[m - 1] は下側の終点なので含めない
    int upper_start = m - 2;
    if (keep_collinear) {
        // すべての点が一直線上にあると、下側と上側で同じ点を往復してしまうので、先に返す
        bool collinear = true;
        for (int i = 2; i < m; ++i) {
            if (sgn(cross(p[ord[1]] - p[ord[0]], p[ord[i]] - p[ord[0]])) != 0) {
                collinear = false;
                break;
            }
        }
        if (collinear) {
            return ord;
        }
        // 右端の縦の辺上の点は、下側で既に通っているので飛ばす
        while (upper_start >= 0 &&
               sgn(p[ord[upper_start]].x - p[ord[m - 1]].x) == 0) {
            --upper_start;
        }
    }

    std::vector<int> ch(2 * m);
    int k = 0;
    for (int i = 0; i < m; ++i) {
        while (k >= 2 && bad(ch[k - 2], ch[k - 1], ord[i])) {
            --k;
        }
        ch[k++] = ord[i];
    }
    for (int i = upper_start, t = k + 1; i >= 0; --i) {
        while (k >= t && bad(ch[k - 2], ch[k - 1], ord[i])) {
            --k;
        }
        ch[k++] = ord[i];
    }
    // 最後の点は、始点 ord[0] に戻ってきたものなので除く
    ch.resize(k - 1);
    return ch;
}

}  // namespace internal

// 凸包を、(x, y) の辞書順で最小の点から反時計回りに返す。重複する点（eps 込みで等しい点）は 1 つにまとめる。
// keep_collinear = false: 頂点だけを返す。すべての点が一直線上にあるときは両端の 2 点、
//   すべての点が同じときは 1 点になる。
// keep_collinear = true: 辺上の点も含め、境界上の点をすべて返す。
//   すべての点が一直線上にあるときは、辞書順に並べて返す。
template<class T, class W>
std::vector<BasicPoint<T, W>> convex_hull(
    const std::vector<BasicPoint<T, W>>& p, const bool keep_collinear = false)
{
    std::vector<BasicPoint<T, W>> res;
    for (int i : internal::convex_hull_index(p, keep_collinear)) {
        res.push_back(p[i]);
    }
    return res;
}

// 凸多角形 p, q のミンコフスキー和。p, q は反時計回りで、一直線上に並ぶ 3 点を含まないこと
// （keep_collinear = false の convex_hull の結果など）。1 点や 2 点の多角形でもよい。
// 結果は辞書順最小の点から反時計回りに並ぶ
template<class T, class W>
std::vector<BasicPoint<T, W>> minkowski_sum(std::vector<BasicPoint<T, W>> p,
                                            std::vector<BasicPoint<T, W>> q)
{
    using P = BasicPoint<T, W>;
    if (p.empty() || q.empty()) {
        return {};
    }
    // 2 つの多角形の辺を偏角の順に併合するため、どちらも辞書順最小の点から始める
    auto reorder = [](std::vector<P>& poly) {
        std::rotate(poly.begin(),
                    std::min_element(poly.begin(), poly.end(), lex_less<T, W>),
                    poly.end());
    };
    reorder(p);
    reorder(q);
    int n = p.size(), m = q.size();
    // i = n や j = m のときにも次の辺を求められるように、先頭の 2 点を末尾に足す
    p.push_back(p[0]);
    p.push_back(p[1 % n]);
    q.push_back(q[0]);
    q.push_back(q[1 % m]);

    std::vector<P> res;
    // 偏角が小さい方の辺を進む。平行なら両方とも進む
    int i = 0, j = 0;
    while (i < n || j < m) {
        res.push_back(p[i] + q[j]);
        int s = sgn(cross(p[i + 1] - p[i], q[j + 1] - q[j]));
        if (s >= 0 && i < n) {
            ++i;
        }
        if (s <= 0 && j < m) {
            ++j;
        }
    }
    return res;
}

// 最遠点対を (距離の 2 乗, (i, j)) で返す。i, j は p の添字。O(n log n)
// p は凸でなくてよいが、空でないこと。すべての点が同じときは i == j になる
template<class T, class W>
std::pair<W, std::pair<int, int>> farthest_pair(
    const std::vector<BasicPoint<T, W>>& p)
{
    assert(!p.empty());
    std::vector<int> idx = internal::convex_hull_index(p);
    const int n = idx.size();
    if (n == 1) {
        return {W(0), {idx[0], idx[0]}};
    }
    if (n == 2) {
        return {norm(p[idx[0]] - p[idx[1]]), {idx[0], idx[1]}};
    }

    auto at = [&](int k) -> const BasicPoint<T, W>& { return p[idx[k]]; };
    // 辞書順で最小の点と最大の点は対蹠点の組なので、そこから始めてキャリパーを半周回す
    int i = 0, j = 0;
    for (int k = 0; k < n; ++k) {
        if (lex_less(at(k), at(i))) {
            i = k;
        }
        if (lex_less(at(j), at(k))) {
            j = k;
        }
    }
    W d = norm(at(i) - at(j));
    int a = i, b = j, si = i, sj = j;
    while (i != sj || j != si) {
        W cur = norm(at(i) - at(j));
        if (d < cur) {
            d = cur;
            a = i;
            b = j;
        }
        if (sgn(cross(at((i + 1) % n) - at(i), at((j + 1) % n) - at(j))) < 0) {
            i = (i + 1) % n;
        }
        else {
            j = (j + 1) % n;
        }
    }
    return {d, {idx[a], idx[b]}};
}

}  // namespace geo2d

}  // namespace rklib
