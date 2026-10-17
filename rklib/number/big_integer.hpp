#pragma once
// https://github.com/Ricky-pon/rk-library-cpp/blob/main/rklib/number/big_integer.hpp

#include <algorithm>
#include <cassert>
#include <compare>
#include <concepts>
#include <cstddef>
#include <iostream>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

#include <atcoder/convolution>

#include "rklib/number/wide.hpp"
#include "rklib/utility/utility.hpp"

namespace rklib {

struct BigInteger {
public:
    BigInteger() : BigInteger(0) {}
    template<std::integral T>
    BigInteger(T x) : BigInteger(from_digits(std::vector<wide_t<T>>{x}))
    {
    }
    explicit BigInteger(std::string_view s) : BigInteger(from_string(s)) {}

    BigInteger& operator+=(const BigInteger& rhs) { return add(rhs, rhs.sgn_); }
    BigInteger& operator-=(const BigInteger& rhs)
    {
        return add(rhs, -rhs.sgn_);
    }
    BigInteger& operator*=(const BigInteger& rhs)
    {
        int rhs_sign = rhs.sgn_;
        if (sgn_ == 0) {
            return *this;
        }
        if (rhs_sign == 0) {
            return *this = rhs;
        }

        int res_sgn = sgn_ * rhs_sign;

        if (int(std::min(v_.size(), rhs.v_.size())) <=
            naive_multiply_threshold) {
            v_ = multiply_naive(v_, rhs.v_);
        }
        else if (std::max(v_.size(), rhs.v_.size()) > ntt_multiply_threshold &&
                 v_.size() + rhs.v_.size() - 1 <= (1 << 24)) {
            v_ = multiply_ntt(v_, rhs.v_);
        }
        else {
            std::vector<long long> a(v_.size()), b(rhs.v_.size());
            for (size_t i = 0; i < a.size(); i++) {
                a[i] = (long long)v_[i];
            }
            for (size_t i = 0; i < b.size(); i++) {
                b[i] = (long long)rhs.v_[i];
            }
            int sz = 1;
            while (sz < int(a.size()) || sz < int(b.size())) {
                sz *= 2;
            }
            a.resize(sz, 0);
            b.resize(sz, 0);

            this->v_ = carry_digits(multiply_karatsuba(a, b));
        }

        this->sgn_ = res_sgn;

        return *this;
    }
    BigInteger& operator/=(const BigInteger& rhs)
    {
        int rhs_sign = rhs.sgn_;
        assert(rhs_sign != 0);
        if (sgn_ == 0) {
            return *this;
        }
        if (rhs == BigInteger(1) || rhs == BigInteger(-1)) {
            sgn_ *= rhs_sign;
            return *this;
        }

        int res_sgn = sgn_ * rhs_sign;

        (*this) = divide_newton_raphson(*this, rhs);

        this->sgn_ = (this->v_.back() == 0 ? 0 : res_sgn);

        return *this;
    }
    BigInteger& operator%=(const BigInteger& rhs)
    {
        assert(rhs.sgn_ != 0);

        return *this = (*this) - (*this) / rhs * rhs;
    }

    friend BigInteger operator+(const BigInteger& lhs, const BigInteger& rhs)
    {
        return BigInteger(lhs) += rhs;
    }
    friend BigInteger operator-(const BigInteger& lhs, const BigInteger& rhs)
    {
        return BigInteger(lhs) -= rhs;
    }
    friend BigInteger operator*(const BigInteger& lhs, const BigInteger& rhs)
    {
        return BigInteger(lhs) *= rhs;
    }
    friend BigInteger operator/(const BigInteger& lhs, const BigInteger& rhs)
    {
        return BigInteger(lhs) /= rhs;
    }
    friend BigInteger operator%(const BigInteger& lhs, const BigInteger& rhs)
    {
        return BigInteger(lhs) %= rhs;
    }

    BigInteger operator+() const { return *this; }
    BigInteger operator-() const { return {-sgn_, v_}; }

    BigInteger& operator++() { return *this += 1; }
    BigInteger& operator--() { return *this -= 1; }

    friend bool operator==(const BigInteger& lhs, const BigInteger& rhs)
    {
        if (lhs.sgn_ != rhs.sgn_ || lhs.v_.size() != rhs.v_.size()) {
            return false;
        }
        for (size_t i = 0; i < lhs.v_.size(); i++) {
            if (lhs.v_[i] != rhs.v_[i]) {
                return false;
            }
        }
        return true;
    }
    friend std::strong_ordering operator<=>(const BigInteger& lhs,
                                            const BigInteger& rhs)
    {
        int lhs_sign = lhs.sgn_, rhs_sign = rhs.sgn_;
        if (lhs_sign < rhs_sign) {
            return std::strong_ordering::less;
        }
        if (lhs_sign > rhs_sign) {
            return std::strong_ordering::greater;
        }

        int nl = lhs.v_.size(), nr = rhs.v_.size();
        if (lhs_sign * nl < rhs_sign * nr) {
            return std::strong_ordering::less;
        }
        if (lhs_sign * nl > rhs_sign * nr) {
            return std::strong_ordering::greater;
        }

        for (int i = nl - 1; i >= 0; i--) {
            if (lhs_sign * lhs.v_[i] < rhs_sign * rhs.v_[i]) {
                return std::strong_ordering::less;
            }
            if (lhs_sign * lhs.v_[i] > rhs_sign * rhs.v_[i]) {
                return std::strong_ordering::greater;
            }
        }

        return std::strong_ordering::equal;
    }

    friend std::istream& operator>>(std::istream& is, BigInteger& rhs)
    {
        std::string s;
        is >> s;
        rhs = BigInteger(s);
        return is;
    }
    friend std::ostream& operator<<(std::ostream& os, const BigInteger& rhs)
    {
        if (rhs.sgn_ < 0) {
            os << "-";
        }
        for (int i = int(rhs.v_.size()) - 1; i >= 0; i--) {
            if (i == int(rhs.v_.size()) - 1) {
                os << rhs.v_[i];
            }
            else {
                os << std::to_string(rhs.v_[i] + base).substr(1, log10_base);
            }
        }
        return os;
    }

    explicit operator double() const
    {
        double res = 0;
        for (int i = v_.size() - 1; i >= 0; i--) {
            res = res * base + v_[i];
        }
        res *= sgn_;
        return res;
    }
    explicit operator long double() const
    {
        long double res = 0;
        for (int i = v_.size() - 1; i >= 0; i--) {
            res = res * base + v_[i];
        }
        res *= sgn_;
        return res;
    }

    friend BigInteger abs(BigInteger x)
    {
        if (x.sgn_ < 0) {
            x.sgn_ = 1;
        }
        return x;
    }

    friend BigInteger gcd(BigInteger a, BigInteger b)
    {
        a = abs(std::move(a));
        b = abs(std::move(b));
        while (b.sgn_ != 0) {
            a %= b;
            std::swap(a, b);
        }
        return a;
    }

private:
    static constexpr int base = 10000, log10_base = 4;
    static constexpr int naive_multiply_threshold = 75,
                         ntt_multiply_threshold = 4000;
    int sgn_;
    std::vector<int> v_;

    BigInteger(int sgn, std::vector<int> v) : sgn_(sgn), v_(std::move(v)) {}

    // 各桁の値が a (下の位から)である数を作る。
    // ただし、0 でない最上位の桁の符号が、数全体の符号と一致すること
    template<std::integral T>
    static BigInteger from_digits(std::vector<T> a)
    {
        while (a.size() >= 2 && a.back() == 0) {
            a.pop_back();
        }
        const int sign = int(a.back() > 0) - int(a.back() < 0);
        if (sign < 0) {
            std::transform(a.begin(), a.end(), a.begin(),
                           [](T x) { return -x; });
        }
        return BigInteger(sign, carry_digits(a));
    }

    // 繰り上がりを処理して、各桁を [0, base) にする。
    template<std::integral T>
    static std::vector<int> carry_digits(const std::vector<T>& a)
    {
        std::vector<int> res(a.size());
        T carry = 0;
        for (int i = 0; i < int(a.size()); ++i) {
            const T tmp = a[i] + carry;
            carry = div_floor(tmp, T(base));
            res[i] = int(tmp - carry * base);
        }
        while (carry > 0) {
            res.push_back(int(carry % base));
            carry /= base;
        }
        while (res.size() >= 2 && res.back() == 0) {
            res.pop_back();
        }
        return res;
    }

    static BigInteger from_string(std::string_view s)
    {
        int sign = 1;
        if (!s.empty() && (s[0] == '+' || s[0] == '-')) {
            sign = (s[0] == '-' ? -1 : 1);
            s.remove_prefix(1);
        }
        assert(!s.empty());

        // 下の位から log10_base 文字ずつ区切る
        std::vector<int> digits;
        for (int r = int(s.size()); r > 0; r -= log10_base) {
            int d = 0;
            for (int i = std::max(0, r - log10_base); i < r; ++i) {
                assert('0' <= s[i] && s[i] <= '9');
                d = d * 10 + (s[i] - '0');
            }
            digits.push_back(sign * d);
        }
        return from_digits(std::move(digits));
    }

    BigInteger& add(const BigInteger& rhs, const int rhs_sign)
    {
        if (rhs_sign == 0) {
            return *this;
        }
        if (sgn_ == 0) {
            *this = rhs;
            sgn_ = rhs_sign;
            return *this;
        }
        if (v_.size() < rhs.v_.size()) {
            v_.resize(rhs.v_.size(), 0);
        }
        for (size_t i = 0; i < rhs.v_.size(); i++) {
            if (sgn_ == rhs_sign) {
                v_[i] += rhs.v_[i];
            }
            else {
                v_[i] -= rhs.v_[i];
            }
        }
        const int sign = sgn_;
        *this = from_digits(std::move(v_));
        sgn_ *= sign;
        return *this;
    }

    static std::vector<int> multiply_naive(const std::vector<int>& a,
                                           const std::vector<int>& b)
    {
        std::vector<long long> c(a.size() + b.size() - 1, 0);
        for (size_t i = 0; i < a.size(); i++) {
            for (size_t j = 0; j < b.size(); j++) {
                c[i + j] += a[i] * b[j];
            }
        }

        return carry_digits(c);
    }
    static std::vector<int> multiply_ntt(const std::vector<int>& a,
                                         const std::vector<int>& b)
    {
        std::vector<long long> al(a.size()), bl(b.size());
        for (size_t i = 0; i < a.size(); i++) {
            al[i] = (long long)a[i];
        }
        for (size_t i = 0; i < b.size(); i++) {
            bl[i] = (long long)b[i];
        }

        std::vector<long long> c = atcoder::convolution_ll(al, bl);

        return carry_digits(c);
    }
    static std::vector<long long> multiply_karatsuba(std::vector<long long> a,
                                                     std::vector<long long> b)
    {
        int na = a.size(), nb = b.size();
        if (std::min(na, nb) <= naive_multiply_threshold) {
            std::vector<long long> res(a.size() + b.size() - 1, 0);
            for (size_t i = 0; i < a.size(); i++) {
                for (size_t j = 0; j < b.size(); j++) {
                    res[i + j] += a[i] * b[j];
                }
            }
            return res;
        }

        int n = std::max(na, nb);
        if (na < n) {
            a.resize(n, 0);
        }
        if (nb < n) {
            b.resize(n, 0);
        }

        int k = n / 2;
        std::vector<long long> x(
            std::vector<long long>(a.begin() + k, a.end()));
        std::vector<long long> y(
            std::vector<long long>(a.begin(), a.begin() + k));
        std::vector<long long> z(
            std::vector<long long>(b.begin() + k, b.end()));
        std::vector<long long> w(
            std::vector<long long>(b.begin(), b.begin() + k));

        std::vector<long long> xz = multiply_karatsuba(x, z);
        std::vector<long long> yw = multiply_karatsuba(y, w);

        for (int i = 0; i < k; i++) {
            x[i] += y[i];
            z[i] += w[i];
        }
        std::vector<long long> t = multiply_karatsuba(x, z);
        for (size_t i = 0; i < xz.size(); i++) {
            t[i] -= xz[i] + yw[i];
        }

        a = std::vector<long long>(2 * k, 0);
        a.insert(a.end(), xz.begin(), xz.end());
        b = std::vector<long long>(k, 0);
        b.insert(b.end(), t.begin(), t.end());

        for (size_t i = 0; i < b.size(); i++) {
            a[i] += b[i];
        }
        for (size_t i = 0; i < yw.size(); i++) {
            a[i] += yw[i];
        }

        return a;
    }

    static void shift_right(BigInteger& a, size_t len)
    {
        if (a == BigInteger(0)) {
            return;
        }
        if (a.v_.size() <= len) {
            a = 0;
            return;
        }
        a.v_.erase(a.v_.begin(), a.v_.begin() + len);
    }

    static BigInteger divide_newton_raphson(BigInteger a, BigInteger b)
    {
        a.sgn_ = b.sgn_ = 1;
        if (a < b) {
            return 0;
        }
        size_t len = a.v_.size() + 2;
        BigInteger x = base * base, prv = x;
        while (true) {
            BigInteger tmp = x * x * b;
            shift_right(tmp, len);
            x += x;
            x -= tmp;
            if (x.v_.size() > len) {
                shift_right(x, x.v_.size() - len);
            }

            if ((x - prv).v_.size() <= 1) {
                break;
            }
            prv = x;
        }
        BigInteger res = a * x;
        res.v_.erase(res.v_.begin(), res.v_.begin() + len);
        if (res * b + b <= a) {
            ++res;
        }

        return res;
    }
};

}  // namespace rklib
