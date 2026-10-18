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
    Rational() : num_(0), den_(1) {}
    Rational(T n) : num_(n), den_(1) {}
    // T が BigInteger などの場合、int -> T -> Rational の 2 段の変換は暗黙には
    // 行われないため、整数から直接作れるようにする
    template<std::integral U>
        requires(!std::same_as<U, T> && std::constructible_from<T, U>)
    Rational(U n) : Rational(T(n))
    {
    }
    Rational(T n, T d) : num_(n), den_(d)
    {
        assert(den_ != T(0));
        if (den_ < T(0)) {
            num_ = -num_, den_ = -den_;
        }
        reduce();
    }
    // 小数は黙って切り捨てられてしまうため禁止
    template<std::floating_point U>
    Rational(U) = delete;

    T numerator() const { return num_; }
    T denominator() const { return den_; }

    Rational& operator+=(const Rational& rhs)
    {
        if (den_ == rhs.den_) {
            num_ += rhs.num_;
            reduce();
            return *this;
        }
        T r_den = rhs.den_;
        if (den_ == T(1) || r_den == T(1)) {
            // 分母の一方が 1 なら、結果はすでに既約
            num_ = num_ * r_den + rhs.num_ * den_;
            den_ *= r_den;
            return *this;
        }
        T g = gcd_of(den_, r_den);
        den_ /= g;
        num_ = num_ * (r_den / g) + rhs.num_ * den_;
        g = gcd_of(num_, g);
        num_ /= g;
        den_ *= r_den / g;
        return *this;
    }
    Rational& operator-=(const Rational& rhs)
    {
        if (den_ == rhs.den_) {
            num_ -= rhs.num_;
            reduce();
            return *this;
        }
        T r_den = rhs.den_;
        if (den_ == T(1) || r_den == T(1)) {
            // 分母の一方が 1 なら、結果はすでに既約
            num_ = num_ * r_den - rhs.num_ * den_;
            den_ *= r_den;
            return *this;
        }
        T g = gcd_of(den_, r_den);
        den_ /= g;
        num_ = num_ * (r_den / g) - rhs.num_ * den_;
        g = gcd_of(num_, g);
        num_ /= g;
        den_ *= r_den / g;
        return *this;
    }
    Rational& operator*=(const Rational& rhs)
    {
        const T r_num = rhs.num_;
        const T r_den = rhs.den_;

        if (den_ == T(1) && r_den == T(1)) {
            num_ *= r_num;
            return *this;
        }
        if (r_den == T(1)) {
            T g = gcd_of(r_num, den_);
            num_ *= r_num / g;
            den_ /= g;
            return *this;
        }
        if (den_ == T(1)) {
            T g = gcd_of(num_, r_den);
            num_ = (num_ / g) * r_num;
            den_ = r_den / g;
            return *this;
        }

        T g1 = gcd_of(num_, r_den);
        T g2 = gcd_of(r_num, den_);
        num_ = (num_ / g1) * (r_num / g2);
        den_ = (den_ / g2) * (r_den / g1);
        return *this;
    }
    Rational& operator/=(const Rational& rhs)
    {
        const T r_num = rhs.num_;
        const T r_den = rhs.den_;

        assert(r_num != T(0));
        if (num_ == T(0)) {
            return *this;
        }

        T gcd1 = gcd_of(num_, r_num);
        T gcd2 = (den_ == T(1) || r_den == T(1)) ? T(1) : gcd_of(r_den, den_);
        num_ = (num_ / gcd1) * (r_den / gcd2);
        den_ = (den_ / gcd2) * (r_num / gcd1);

        if (den_ < T(0)) {
            num_ = -num_;
            den_ = -den_;
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
    Rational operator-() const { return Rational(-num_, den_, Reduced{}); }

    Rational& operator++()
    {
        num_ += den_;
        return *this;
    }
    Rational& operator--()
    {
        num_ -= den_;
        return *this;
    }

    // 常に約分されているので、分子と分母がそれぞれ等しければ等しい
    friend bool operator==(const Rational&, const Rational&) = default;
    friend std::strong_ordering operator<=>(const Rational& lhs,
                                            const Rational& rhs)
    {
        using W = wide_t<T>;
        const W l = W(lhs.num_) * W(rhs.den_), r = W(rhs.num_) * W(lhs.den_);
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
        if (r.den_ == T(1)) {
            return os << r.num_;
        }
        return os << r.num_ << "/" << r.den_;
    }

    explicit operator double() const { return double(num_) / double(den_); }
    explicit operator long double() const
    {
        return (long double)num_ / (long double)den_;
    }

private:
    T num_, den_;

    // 約分済みの値から作る（gcd の計算を省く）
    struct Reduced {};
    Rational(T n, T d, Reduced) : num_(n), den_(d) {}

    // 分母が正であることを前提に約分する
    void reduce()
    {
        if (den_ == T(1)) {
            return;
        }
        T g = gcd_of(num_, den_);
        num_ /= g;
        den_ /= g;
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
