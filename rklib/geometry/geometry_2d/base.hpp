#pragma once
// https://github.com/Ricky-pon/rk-library-cpp/blob/main/rklib/geometry/geometry_2d/base.hpp

#include <algorithm>
#include <cassert>
#include <cmath>
#include <concepts>
#include <vector>

#include "rklib/geometry/scalar.hpp"

namespace rklib {

namespace geo2d {

// rklib::geo の共通部分を取り込む。
// 点用の eq などを定義するとスカラー用が隠れるので、using namespace ではなく名前ごとに取り込む。
using geo::compare;
using geo::eps;
using geo::eq;
using geo::ExactDivisionScalar;
using geo::FloatScalar;
using geo::geq;
using geo::gt;
using geo::leq;
using geo::lt;
using geo::Scalar;
using geo::scalar_abs;
using geo::ScalarTraits;
using geo::set_eps;
using geo::sgn;

template<Scalar T, class W = typename ScalarTraits<T>::wide>
struct BasicPoint {
    // 浮動小数点数で W を変えると、eps や sqrt の型が二重になるため
    static_assert(FloatScalar<T> ? std::same_as<T, W> : !FloatScalar<W>,
                  "W must be T for floating point T, and exact for exact T");

    T x, y;

    BasicPoint() : x(0), y(0) {}
    BasicPoint(T x, T y) : x(x), y(y) {}

    BasicPoint& operator+=(const BasicPoint& rhs)
    {
        x += rhs.x;
        y += rhs.y;
        return *this;
    }
    BasicPoint& operator-=(const BasicPoint& rhs)
    {
        x -= rhs.x;
        y -= rhs.y;
        return *this;
    }
    BasicPoint& operator*=(const T& rhs)
    {
        x *= rhs;
        y *= rhs;
        return *this;
    }
    BasicPoint& operator/=(const T& rhs)
        requires ExactDivisionScalar<T>
    {
        x /= rhs;
        y /= rhs;
        return *this;
    }

    BasicPoint operator+() const { return *this; }
    BasicPoint operator-() const { return BasicPoint(-x, -y); }

    friend BasicPoint operator+(const BasicPoint& lhs, const BasicPoint& rhs)
    {
        return BasicPoint(lhs) += rhs;
    }
    friend BasicPoint operator-(const BasicPoint& lhs, const BasicPoint& rhs)
    {
        return BasicPoint(lhs) -= rhs;
    }
    friend BasicPoint operator*(const BasicPoint& lhs, const T& rhs)
    {
        return BasicPoint(lhs) *= rhs;
    }
    friend BasicPoint operator*(const T& lhs, const BasicPoint& rhs)
    {
        return BasicPoint(rhs) *= lhs;
    }
    friend BasicPoint operator/(const BasicPoint& lhs, const T& rhs)
        requires ExactDivisionScalar<T>
    {
        return BasicPoint(lhs) /= rhs;
    }
    // 整数の点に小数のスカラー倍を行うと、黙って切り捨てられるため禁止
    template<std::floating_point U>
        requires(!FloatScalar<T>)
    friend void operator*(const BasicPoint&, const U&) = delete;
    template<std::floating_point U>
        requires(!FloatScalar<T>)
    friend void operator*(const U&, const BasicPoint&) = delete;
    template<std::floating_point U>
        requires(!FloatScalar<T>)
    friend void operator/(const BasicPoint&, const U&) = delete;

    // 誤差を考えずに比べる（ソートや map のキー用）。誤差を考えるときは eq を使う
    auto operator<=>(const BasicPoint&) const = default;
};

template<Scalar T, class W = typename ScalarTraits<T>::wide>
using BasicVec = BasicPoint<T, W>;

template<Scalar T, class W = typename ScalarTraits<T>::wide>
struct BasicLine {
    BasicPoint<T, W> a;
    BasicVec<T, W> v;

    BasicLine() {}
    BasicLine(BasicPoint<T, W> a, BasicVec<T, W> v) : a(a), v(v) {}
};

template<Scalar T, class W = typename ScalarTraits<T>::wide>
struct BasicSegment {
    BasicPoint<T, W> a, b;

    BasicSegment() {}
    BasicSegment(BasicPoint<T, W> a, BasicPoint<T, W> b) : a(a), b(b) {}
};

template<Scalar T, class W = typename ScalarTraits<T>::wide>
struct BasicCircle {
    BasicPoint<T, W> c;
    T r;

    BasicCircle() : r(0) {}
    BasicCircle(BasicPoint<T, W> c, T r) : c(c), r(r) {}
};

// よく使う型の別名をまとめたもの。使う側では次のように書く。
//   using G2 = rklib::geo2d::Geo2D<long double>;
//   using Point = G2::Point;
//   rklib::geo::set_eps<long double>(1e-9);
// 別名を名前空間に直接置くと、使う側でつけた Point などと衝突するため、構造体の中に置いている。
template<Scalar T, class W = typename ScalarTraits<T>::wide>
struct Geo2D {
    using Point = BasicPoint<T, W>;
    using Vec = Point;
    using Points = std::vector<Point>;
    using Polygon = Points;
    using Line = BasicLine<T, W>;
    using Segment = BasicSegment<T, W>;
    using Circle = BasicCircle<T, W>;
};

template<class T, class W>
W dot(const BasicVec<T, W>& v, const BasicVec<T, W>& w)
{
    return W(v.x) * W(w.x) + W(v.y) * W(w.y);
}

template<class T, class W>
W cross(const BasicVec<T, W>& v, const BasicVec<T, W>& w)
{
    return W(v.x) * W(w.y) - W(v.y) * W(w.x);
}

template<class T, class W>
bool eq(const BasicPoint<T, W>& p, const BasicPoint<T, W>& q)
{
    return eq(p.x, q.x) && eq(p.y, q.y);
}

// |v|^2
template<class T, class W>
W norm(const BasicVec<T, W>& v)
{
    return dot(v, v);
}

// |v|
template<class T, class W>
    requires FloatScalar<T>
T abs(const BasicVec<T, W>& v)
{
    return std::sqrt(norm(v));
}

// ベクトルの偏角。戻り値の範囲は [-π, π]
template<class T, class W>
    requires FloatScalar<T>
T arg(const BasicVec<T, W>& v)
{
    return std::atan2(v.y, v.x);
}

// v から w へのなす角。反時計回りを正とする
template<class T, class W>
    requires FloatScalar<T>
T angle(const BasicVec<T, W>& v, const BasicVec<T, W>& w)
{
    return arg(BasicPoint<T, W>(dot(v, w), cross(v, w)));
}

template<class T, class W>
    requires FloatScalar<T>
BasicVec<T, W> rotation(const BasicVec<T, W>& v, const T theta)
{
    return BasicVec<T, W>(v.x * std::cos(theta) - v.y * std::sin(theta),
                          v.x * std::sin(theta) + v.y * std::cos(theta));
}

template<class T, class W>
    requires FloatScalar<T>
BasicPoint<T, W> rotation(const BasicPoint<T, W>& p, const T theta,
                          const BasicPoint<T, W>& center)
{
    return center + rotation(p - center, theta);
}

template<class T, class W>
    requires ExactDivisionScalar<T>
BasicPoint<T, W> projection(const BasicPoint<T, W>& p, const BasicLine<T, W>& l)
{
    const auto& [a, v] = l;
    assert((v != BasicVec<T, W>()) && "the direction of l must not be zero");
    return a + v * dot(v, p - a) / norm(v);
}

// 直線 l を挟んで、p と対称な位置にある点
template<class T, class W>
    requires ExactDivisionScalar<T>
BasicPoint<T, W> reflection(const BasicPoint<T, W>& p, const BasicLine<T, W>& l)
{
    return projection(p, l) * 2 - p;
}

// 点 c が、有向線分 ab に対してどのような位置関係にあるか
enum class CCW : int {
    COUNTER_CLOCKWISE = 1,
    CLOCKWISE = -1,
    ON_LINE_BACK = -2,  // c -> a -> b
    ON_LINE_FRONT = 2,  // a -> b -> c
    ON_SEGMENT = 0,     // a -> c -> b (c は a, b と一致しうる)
};

template<class T, class W>
CCW ccw(const BasicPoint<T, W>& a, const BasicPoint<T, W>& b,
        const BasicPoint<T, W>& c)
{
    const BasicVec<T, W> v = b - a, w = c - a;
    const int sign = sgn(cross(v, w));
    if (sign > 0) {
        return CCW::COUNTER_CLOCKWISE;
    }
    if (sign < 0) {
        return CCW::CLOCKWISE;
    }
    if (sgn(dot(v, w)) < 0) {
        return CCW::ON_LINE_BACK;
    }
    if (sgn(dot(a - b, c - b)) < 0) {
        return CCW::ON_LINE_FRONT;
    }
    return CCW::ON_SEGMENT;
}

template<class T, class W>
bool is_parallel(const BasicVec<T, W>& v, const BasicVec<T, W>& w)
{
    return sgn(cross(v, w)) == 0;
}

template<class T, class W>
bool is_orthogonal(const BasicVec<T, W>& v, const BasicVec<T, W>& w)
{
    return sgn(dot(v, w)) == 0;
}

// 偏角を (-π, π] の範囲で比較する。
// 原点は偏角 0 として扱い、正の x 軸の手前に置く
template<class T, class W>
bool arg_less(const BasicVec<T, W>& v, const BasicVec<T, W>& w)
{
    // 0: (-π, 0), 1: origin, 2: [0, π), 3: π
    auto region = [](const BasicVec<T, W>& u) {
        const int sx = sgn(u.x), sy = sgn(u.y);
        if (sy < 0) {
            return 0;
        }
        if (sy == 0 && sx == 0) {
            return 1;
        }
        if (sy > 0 || sx > 0) {
            return 2;
        }
        return 3;
    };

    const int rv = region(v), rw = region(w);
    if (rv != rw) {
        return rv < rw;
    }
    return sgn(cross(v, w)) > 0;
}

template<class T, class W>
void arg_sort(std::vector<BasicVec<T, W>>& vs)
{
    std::sort(vs.begin(), vs.end(), arg_less<T, W>);
}

}  // namespace geo2d

}  // namespace rklib
