#pragma once
// https://github.com/Ricky-pon/rk-library-cpp/blob/main/rklib/geometry/scalar.hpp

#include <concepts>
#include <type_traits>

namespace rklib {

template<class T>
struct Rational;

namespace geo {

// 座標に使える数の型。int, long long, __int128, Rational, double, long double など
template<class T>
concept Scalar = std::totally_ordered<T> && requires(T a, T b) {
    T(0);
    { a + b } -> std::convertible_to<T>;
    { a - b } -> std::convertible_to<T>;
    { a * b } -> std::convertible_to<T>;
    { a += b } -> std::convertible_to<T>;
    { a -= b } -> std::convertible_to<T>;
    { a *= b } -> std::convertible_to<T>;
    { -a } -> std::convertible_to<T>;
};

// is_float           : 浮動小数点数か（比較に eps を使うか）
// has_exact_division : 割り算が正確か（交点などを求められるか）
// wide               : 内積と外積を計算する型（オーバーフローを避けるため）
template<Scalar T>
struct ScalarTraits {
    static constexpr bool is_float = std::is_floating_point_v<T>;
    static constexpr bool has_exact_division = std::is_floating_point_v<T>;
    using wide = T;
};

template<class T>
    requires Scalar<T> && std::signed_integral<T> && (sizeof(T) <= 4)
struct ScalarTraits<T> {
    static constexpr bool is_float = false;
    static constexpr bool has_exact_division = false;
    using wide = long long;
};

template<class T>
    requires Scalar<T> && std::signed_integral<T> && (sizeof(T) == 8)
struct ScalarTraits<T> {
    static constexpr bool is_float = false;
    static constexpr bool has_exact_division = false;
    using wide = __int128;
};

template<class U>
struct ScalarTraits<Rational<U>> {
    static constexpr bool is_float = false;
    static constexpr bool has_exact_division = true;
    using wide = Rational<U>;
};

template<class T>
concept ExactDivisionScalar = Scalar<T> && ScalarTraits<T>::has_exact_division;

template<class T>
concept FloatScalar = ExactDivisionScalar<T> && ScalarTraits<T>::is_float;

// 絶対誤差。長さにも内積や外積にも同じ値を使う
template<FloatScalar T>
inline T eps = T(1e-9);

// 型の指定を必須にしている。
// 引数から型を決めた場合、eps<long double> を書き換えたい場面で set_eps(1e-7) と書くと、
// eps<double> が変わってしまう。
template<FloatScalar T>
void set_eps(std::type_identity_t<T> new_eps)
{
    eps<T> = new_eps;
}

template<class T>
int sgn(const T& x)
{
    if constexpr (FloatScalar<T>) {
        if (x < -eps<T>) {
            return -1;
        }
        if (x > eps<T>) {
            return 1;
        }
        return 0;
    }
    else {
        if (x < T(0)) {
            return -1;
        }
        if (x > T(0)) {
            return 1;
        }
        return 0;
    }
}

// geo2d などの中では点用の abs が見つかり、std::abs が使われないため、
// スカラーにはこちらを使う。
template<class T>
T scalar_abs(const T& x)
{
    return x < T(0) ? -x : x;
}

// x < y なら -1, x == y なら 0, x > y なら 1
template<Scalar T, Scalar U>
    requires std::common_with<T, U>
int compare(const T& x, const U& y)
{
    using C = std::common_type_t<T, U>;
    if constexpr (FloatScalar<C>) {
        return sgn(C(x) - C(y));
    }
    else {
        // x - y がオーバーフローしうるため、引き算せずに比べる
        return int(C(y) < C(x)) - int(C(x) < C(y));
    }
}

// x == y
template<Scalar T, Scalar U>
    requires std::common_with<T, U>
bool eq(const T& x, const U& y)
{
    return compare(x, y) == 0;
}

// x < y
template<Scalar T, Scalar U>
    requires std::common_with<T, U>
bool lt(const T& x, const U& y)
{
    return compare(x, y) < 0;
}

// x > y
template<Scalar T, Scalar U>
    requires std::common_with<T, U>
bool gt(const T& x, const U& y)
{
    return compare(x, y) > 0;
}

// x <= y
template<Scalar T, Scalar U>
    requires std::common_with<T, U>
bool leq(const T& x, const U& y)
{
    return compare(x, y) <= 0;
}

// x >= y
template<Scalar T, Scalar U>
    requires std::common_with<T, U>
bool geq(const T& x, const U& y)
{
    return compare(x, y) >= 0;
}

}  // namespace geo

}  // namespace rklib
