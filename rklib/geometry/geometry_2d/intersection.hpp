#pragma once
// https://github.com/Ricky-pon/rk-library-cpp/blob/main/rklib/geometry/geometry_2d/intersection.hpp

#include <cassert>

#include "rklib/geometry/geometry_2d/base.hpp"

namespace rklib {

namespace geo2d {

// 直線 l と線分 s が共有点を持つか。端点が l 上にある場合、線分が直線に重なる場合も含む
template<class T, class W>
bool has_intersection(const BasicLine<T, W>& l, const BasicSegment<T, W>& s)
{
    const auto [p, v] = l;
    assert((v != BasicPoint<T, W>()) && "the line must not be degenerated");
    const auto [a, b] = s;
    return sgn(cross(v, a - p)) * sgn(cross(v, b - p)) <= 0;
}

template<class T, class W>
bool has_intersection(const BasicSegment<T, W>& s1,
                      const BasicSegment<T, W>& s2)
{
    const auto [a, b] = s1;
    const auto [c, d] = s2;
    if (eq(a, b) && eq(c, d)) {
        return eq(a, c);
    }
    if (eq(a, b)) {
        return ccw(c, d, a) == CCW::ON_SEGMENT;
    }
    if (eq(c, d)) {
        return ccw(a, b, c) == CCW::ON_SEGMENT;
    }
    auto f = [](CCW x) { return static_cast<int>(x); };
    return f(ccw(a, b, c)) * f(ccw(a, b, d)) <= 0 &&
           f(ccw(c, d, a)) * f(ccw(c, d, b)) <= 0;
}

template<class T, class W>
    requires ExactDivisionScalar<T>
BasicPoint<T, W> intersection(const BasicLine<T, W>& l1,
                              const BasicLine<T, W>& l2)
{
    const auto [a, v] = l1;
    const auto [b, w] = l2;
    assert(!is_parallel(v, w) && "the lines must not be parallel");
    return a + v * T(cross(b - a, w) / cross(v, w));
}

// 線分の交点。s1 と s2 は交わり、かつ平行でない（重なる線分は扱えない）
template<class T, class W>
    requires ExactDivisionScalar<T>
BasicPoint<T, W> intersection(const BasicSegment<T, W>& s1,
                              const BasicSegment<T, W>& s2)
{
    const auto [a, b] = s1;
    const auto [c, d] = s2;
    assert(has_intersection(s1, s2) && "the segments must intersect");
    return intersection(BasicLine<T, W>(a, b - a), BasicLine<T, W>(c, d - c));
}

}  // namespace geo2d

}  // namespace rklib
