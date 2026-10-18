#pragma once
// https://github.com/Ricky-pon/rk-library-cpp/blob/main/rklib/geometry/geometry_2d/distance.hpp

#include <algorithm>
#include <cassert>
#include <cmath>

#include "rklib/geometry/geometry_2d/base.hpp"
#include "rklib/geometry/geometry_2d/intersection.hpp"

namespace rklib {

namespace geo2d {

template<class T, class W>
W distance_squared(const BasicPoint<T, W>& p, const BasicPoint<T, W>& q)
{
    return norm(p - q);
}

template<class T, class W>
    requires ExactDivisionScalar<T>
W distance_squared(const BasicLine<T, W>& l, const BasicPoint<T, W>& p)
{
    const auto [a, v] = l;
    assert((v != BasicVec<T, W>()) && "the direction of l must not be zero");
    W c = cross(v, p - a);
    return c * c / norm(v);
}

template<class T, class W>
    requires ExactDivisionScalar<T>
W distance_squared(const BasicSegment<T, W>& s, const BasicPoint<T, W>& p)
{
    const auto [a, b] = s;
    if (eq(a, b) || sgn(dot(b - a, p - a)) < 0) {
        return norm(p - a);
    }
    if (sgn(dot(a - b, p - b)) < 0) {
        return norm(p - b);
    }
    return distance_squared(BasicLine<T, W>(a, b - a), p);
}

template<class T, class W>
    requires ExactDivisionScalar<T>
W distance_squared(const BasicLine<T, W>& l1, const BasicLine<T, W>& l2)
{
    if (is_parallel(l1.v, l2.v)) {
        return distance_squared(l1, l2.a);
    }
    return W(0);
}

template<class T, class W>
    requires ExactDivisionScalar<T>
W distance_squared(const BasicLine<T, W>& l, const BasicSegment<T, W>& s)
{
    if (has_intersection(l, s)) {
        return W(0);
    }
    return std::min(distance_squared(l, s.a), distance_squared(l, s.b));
}

template<class T, class W>
    requires ExactDivisionScalar<T>
W distance_squared(const BasicSegment<T, W>& s1, const BasicSegment<T, W>& s2)
{
    if (has_intersection(s1, s2)) {
        return W(0);
    }
    return std::min({distance_squared(s1, s2.a), distance_squared(s1, s2.b),
                     distance_squared(s2, s1.a), distance_squared(s2, s1.b)});
}

template<class T, class W>
    requires FloatScalar<T>
T distance(const BasicPoint<T, W>& p, const BasicPoint<T, W>& q)
{
    return abs(p - q);
}

template<class T, class W>
    requires FloatScalar<T>
T distance(const BasicLine<T, W>& l, const BasicPoint<T, W>& p)
{
    const auto [a, v] = l;
    assert((v != BasicVec<T, W>()) && "the direction of l must not be zero");
    return std::abs(T(cross(v, p - a))) / abs(v);
}

template<class T, class W>
    requires FloatScalar<T>
T distance(const BasicSegment<T, W>& s, const BasicPoint<T, W>& p)
{
    const auto [a, b] = s;
    if (eq(a, b) || sgn(dot(b - a, p - a)) < 0) {
        return abs(p - a);
    }
    if (sgn(dot(a - b, p - b)) < 0) {
        return abs(p - b);
    }
    return distance(BasicLine<T, W>(a, b - a), p);
}

template<class T, class W>
    requires FloatScalar<T>
T distance(const BasicLine<T, W>& l1, const BasicLine<T, W>& l2)
{
    if (is_parallel(l1.v, l2.v)) {
        return distance(l1, l2.a);
    }
    return T(0);
}

template<class T, class W>
    requires FloatScalar<T>
T distance(const BasicLine<T, W>& l, const BasicSegment<T, W>& s)
{
    if (has_intersection(l, s)) {
        return T(0);
    }
    return std::min(distance(l, s.a), distance(l, s.b));
}

template<class T, class W>
    requires FloatScalar<T>
T distance(const BasicSegment<T, W>& s1, const BasicSegment<T, W>& s2)
{
    if (has_intersection(s1, s2)) {
        return T(0);
    }
    return std::min({distance(s1, s2.a), distance(s1, s2.b), distance(s2, s1.a),
                     distance(s2, s1.b)});
}

}  // namespace geo2d

}  // namespace rklib