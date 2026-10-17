#pragma once
// https://github.com/Ricky-pon/rk-library-cpp/blob/main/rklib/number/rational.hpp

#include <cassert>
#include <compare>
#include <concepts>
#include <numeric>
#include <ostream>

#include "rklib/number/wide.hpp"

namespace rklib {

// 途中の積が T に収まる範囲で正しい。比較だけは wide_t<T> で計算する
template<class T = long long>
struct Rational {
    Rational() : num(0), den(1) {}
    Rational(T n) : num(n), den(1) {}
    Rational(T n, T d) : num(n), den(d)
    {
        assert(den != T(0));
        if (den < T(0)) {
            num = -num, den = -den;
        }
        reduce();
    }
    // 小数は黙って切り捨てられてしまうため禁止
    template<std::floating_point U>
    Rational(U) = delete;

    T numerator() const { return num; }
    T denominator() const { return den; }

    Rational& operator+=(const Rational& rhs)
    {
        if (den == rhs.den) {
            num += rhs.num;
            reduce();
            return *this;
        }
        T r_den = rhs.den;
        if (den == T(1) || r_den == T(1)) {
            // 分母の一方が 1 なら、結果はすでに既約
            num = num * r_den + rhs.num * den;
            den *= r_den;
            return *this;
        }
        T g = gcd_of(den, r_den);
        den /= g;
        num = num * (r_den / g) + rhs.num * den;
        g = gcd_of(num, g);
        num /= g;
        den *= r_den / g;
        return *this;
    }
    Rational& operator-=(const Rational& rhs)
    {
        if (den == rhs.den) {
            num -= rhs.num;
            reduce();
            return *this;
        }
        T r_den = rhs.den;
        if (den == T(1) || r_den == T(1)) {
            // 分母の一方が 1 なら、結果はすでに既約
            num = num * r_den - rhs.num * den;
            den *= r_den;
            return *this;
        }
        T g = gcd_of(den, r_den);
        den /= g;
        num = num * (r_den / g) - rhs.num * den;
        g = gcd_of(num, g);
        num /= g;
        den *= r_den / g;
        return *this;
    }
    Rational& operator*=(const Rational& rhs)
    {
        T r_num = rhs.num;
        T r_den = rhs.den;

        if (den == T(1) && r_den == T(1)) {
            num *= r_num;
            return *this;
        }
        if (r_den == T(1)) {
            T g = gcd_of(r_num, den);
            num *= r_num / g;
            den /= g;
            return *this;
        }
        if (den == T(1)) {
            T g = gcd_of(num, r_den);
            num = (num / g) * r_num;
            den = r_den / g;
            return *this;
        }

        T g1 = gcd_of(num, r_den);
        T g2 = gcd_of(r_num, den);
        num = (num / g1) * (r_num / g2);
        den = (den / g2) * (r_den / g1);
        return *this;
    }
    Rational& operator/=(const Rational& rhs)
    {
        T r_num = rhs.num;
        T r_den = rhs.den;

        assert(r_num != T(0));
        if (num == T(0)) {
            return *this;
        }

        T gcd1 = gcd_of(num, r_num);
        T gcd2 = (den == T(1) || r_den == T(1)) ? T(1) : gcd_of(r_den, den);
        num = (num / gcd1) * (r_den / gcd2);
        den = (den / gcd2) * (r_num / gcd1);

        if (den < T(0)) {
            num = -num;
            den = -den;
        }
        return *this;
    }

    friend Rational operator+(const Rational& lhs, const Rational& rhs)
    {
        return Rational(lhs) += rhs;
    }
    friend Rational operator-(const Rational& lhs, const Rational& rhs)
    {
        return Rational(lhs) -= rhs;
    }
    friend Rational operator*(const Rational& lhs, const Rational& rhs)
    {
        return Rational(lhs) *= rhs;
    }
    friend Rational operator/(const Rational& lhs, const Rational& rhs)
    {
        return Rational(lhs) /= rhs;
    }

    Rational operator+() const { return *this; }
    Rational operator-() const { return Rational(-num, den, Reduced{}); }

    Rational& operator++()
    {
        num += den;
        return *this;
    }
    Rational& operator--()
    {
        num -= den;
        return *this;
    }

    // 常に約分されているので、分子と分母がそれぞれ等しければ等しい
    friend bool operator==(const Rational&, const Rational&) = default;
    friend std::strong_ordering operator<=>(const Rational& lhs,
                                            const Rational& rhs)
    {
        using W = wide_t<T>;
        const W l = W(lhs.num) * W(rhs.den), r = W(rhs.num) * W(lhs.den);
        if (l < r) {
            return std::strong_ordering::less;
        }
        if (r < l) {
            return std::strong_ordering::greater;
        }
        return std::strong_ordering::equal;
    }

    friend std::ostream& operator<<(std::ostream& os, const Rational& r)
    {
        if (r.den == T(1)) {
            return os << r.num;
        }
        return os << r.num << "/" << r.den;
    }

    explicit operator double() const { return double(num) / double(den); }
    explicit operator long double() const
    {
        return (long double)num / (long double)den;
    }

private:
    T num, den;

    // 約分済みの値から作る（gcd の計算を省く）
    struct Reduced {};
    Rational(T n, T d, Reduced) : num(n), den(d) {}

    // 分母が正であることを前提に約分する
    void reduce()
    {
        if (den == T(1)) {
            return;
        }
        T g = gcd_of(num, den);
        num /= g;
        den /= g;
    }

    static T gcd_of(const T& a, const T& b)
    {
        if constexpr (std::integral<T>) {
            return std::gcd(a, b);
        }
        else {
            return gcd(a, b);
        }
    }
};

}  // namespace rklib
