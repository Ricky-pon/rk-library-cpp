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
    // bool は組み込みの整数と同じく 0, 1 として扱う。
    // wide_t<bool> は bool のままで繰り上がりを計算できないため、int にしてから作る
    template<std::integral T>
    BigInteger(T x)
        : BigInteger(from_digits(
              std::vector<
                  wide_t<std::conditional_t<std::same_as<T, bool>, int, T>>>{
                  x}))
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
        else if (std::min(v_.size(), rhs.v_.size()) > ntt_multiply_threshold) {
            v_ = multiply_ntt(v_, rhs.v_);
        }
        else {
            v_ = carry_digits(multiply_karatsuba(
                std::vector<long long>(v_.begin(), v_.end()),
                std::vector<long long>(rhs.v_.begin(), rhs.v_.end())));
        }

        this->sgn_ = res_sgn;

        return *this;
    }
    BigInteger& operator/=(const BigInteger& rhs)
    {
        *this = divmod(*this, rhs).first;
        return *this;
    }
    BigInteger& operator%=(const BigInteger& rhs)
    {
        *this = divmod(*this, rhs).second;
        return *this;
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
    // 商と余り。商は 0 方向に丸め、余りの符号は lhs と同じ（組み込みの整数と同じ）
    friend std::pair<BigInteger, BigInteger> divmod(const BigInteger& lhs,
                                                    const BigInteger& rhs)
    {
        assert(rhs.sgn_ != 0);
        if (lhs.sgn_ == 0) {
            return {0, 0};
        }
        auto [q, r] = divmod_abs(abs(lhs), abs(rhs));
        q.sgn_ *= lhs.sgn_ * rhs.sgn_;
        r.sgn_ *= lhs.sgn_;
        return {std::move(q), std::move(r)};
    }

    BigInteger operator+() const { return *this; }
    BigInteger operator-() const { return {-sgn_, v_}; }

    BigInteger& operator++() { return *this += 1; }
    BigInteger& operator--() { return *this -= 1; }

    BigInteger operator++(int)
    {
        BigInteger res = *this;
        ++*this;
        return res;
    }
    BigInteger operator--(int)
    {
        BigInteger res = *this;
        --*this;
        return res;
    }

    friend bool operator==(const BigInteger& lhs,
                           const BigInteger& rhs) = default;
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
        if (is >> s) {
            rhs = BigInteger(s);
        }
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

    template<std::floating_point T>
    explicit operator T() const
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

    template<class T, class U>
        requires(std::same_as<T, BigInteger> || std::same_as<U, BigInteger>) &&
                std::convertible_to<T, BigInteger> &&
                std::convertible_to<U, BigInteger>
    friend BigInteger gcd(T a, U b)
    {
        BigInteger x = abs(BigInteger(std::move(a)));
        BigInteger y = abs(BigInteger(std::move(b)));
        while (y.sgn_ != 0) {
            x %= y;
            std::swap(x, y);
        }
        return x;
    }

private:
    static constexpr int base = 10000, log10_base = 4;
    static constexpr int naive_multiply_threshold = 75,
                         ntt_multiply_threshold = 2000;
    int sgn_;
    std::vector<int> v_;

    BigInteger(int sgn, std::vector<int> v) : sgn_(sgn), v_(std::move(v)) {}

    // ---- 正規化・文字列からの変換 ----

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

    // ---- 足し算・引き算 ----

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

    // ---- 掛け算 ----

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
        // convolution_ll で計算できる長さの上限
        assert(a.size() + b.size() - 1 <= (1 << 24));
        return carry_digits(atcoder::convolution_ll(
            std::vector<long long>(a.begin(), a.end()),
            std::vector<long long>(b.begin(), b.end())));
    }
    // 短い方が ntt_multiply_threshold 桁以下で使う。
    // a0 + a1 のように足してから掛けるので値は段ごとに最大 2 倍になるが、
    // 半分に分ける段数は log2(ntt_multiply_threshold / naive_multiply_threshold)
    // 程度なので long long に収まる
    static std::vector<long long> multiply_karatsuba(std::vector<long long> a,
                                                     std::vector<long long> b)
    {
        const int na = a.size(), nb = b.size();
        if (na < nb) {
            return multiply_karatsuba(b, a);
        }
        if (nb <= naive_multiply_threshold) {
            std::vector<long long> res(na + nb - 1, 0);
            for (int i = 0; i < na; i++) {
                for (int j = 0; j < nb; j++) {
                    res[i + j] += a[i] * b[j];
                }
            }
            return res;
        }

        // a = a1 base^k + a0, b = b1 base^k + b0 (a0, b0 は k 桁)
        const int k = (na + 1) / 2;
        std::vector<long long> res(na + nb - 1, 0);
        auto add_to = [&](const std::vector<long long>& v, int shift) {
            for (size_t i = 0; i < v.size(); i++) {
                res[shift + i] += v[i];
            }
        };
        const std::vector<long long> a0(a.begin(), a.begin() + k),
            a1(a.begin() + k, a.end());
        if (nb <= k) {
            // b が a の半分以下の長さなら、a だけを分ける
            add_to(multiply_karatsuba(a0, b), 0);
            add_to(multiply_karatsuba(a1, b), k);
            return res;
        }

        const std::vector<long long> b0(b.begin(), b.begin() + k),
            b1(b.begin() + k, b.end());
        const std::vector<long long> z0 = multiply_karatsuba(a0, b0);
        const std::vector<long long> z2 = multiply_karatsuba(a1, b1);
        // z1 = (a0 + a1)(b0 + b1) - z0 - z2
        std::vector<long long> s = a0, t = b0;
        for (size_t i = 0; i < a1.size(); i++) {
            s[i] += a1[i];
        }
        for (size_t i = 0; i < b1.size(); i++) {
            t[i] += b1[i];
        }
        std::vector<long long> z1 = multiply_karatsuba(s, t);
        for (size_t i = 0; i < z0.size(); i++) {
            z1[i] -= z0[i];
        }
        for (size_t i = 0; i < z2.size(); i++) {
            z1[i] -= z2[i];
        }

        add_to(z0, 0);
        add_to(z1, k);
        add_to(z2, 2 * k);
        return res;
    }

    // ---- 割り算 ----

    static constexpr int naive_divide_threshold = 128;

    // a * base^d を返す。d < 0 のときは切り捨て
    static BigInteger shift(BigInteger a, int d)
    {
        if (a.sgn_ == 0 || d == 0) {
            return a;
        }
        if (d > 0) {
            a.v_.insert(a.v_.begin(), d, 0);
            return a;
        }
        if (int(a.v_.size()) <= -d) {
            return 0;
        }
        a.v_.erase(a.v_.begin(), a.v_.begin() - d);
        return a;
    }

    // a, b > 0 の商と余り
    static std::pair<BigInteger, BigInteger> divmod_abs(const BigInteger& a,
                                                        const BigInteger& b)
    {
        if (a < b) {
            return {0, a};
        }
        const int n = a.v_.size(), m = b.v_.size();
        if (m == 1) {
            return divmod_small(a, b.v_[0]);
        }
        if (std::min(m, n - m + 1) <= naive_divide_threshold) {
            return divmod_naive(a, b);
        }
        return divmod_newton(a, b);
    }

    // 1 桁の数 d で割る。O(n)
    //   n: a を base 進数表記したときの桁数
    static std::pair<BigInteger, BigInteger> divmod_small(const BigInteger& a,
                                                          const int d)
    {
        if (d == 1) {
            return {a, 0};
        }
        std::vector<int> q(a.v_.size());
        int r = 0;
        for (int i = int(a.v_.size()) - 1; i >= 0; i--) {
            const int cur = r * base + a.v_[i];
            q[i] = cur / d;
            r = cur % d;
        }
        return {from_digits(std::move(q)), r};
    }

    // 筆算。b は 2 桁以上。O((n - m + 1)m)
    //   n: a を base 進数表記したときの桁数
    //   m: b を base 進数表記したときの桁数
    static std::pair<BigInteger, BigInteger> divmod_naive(const BigInteger& a,
                                                          const BigInteger& b)
    {
        const int n = a.v_.size(), m = b.v_.size();
        const long long b_top = (long long)b.v_[m - 1] * base + b.v_[m - 2];

        // r += coef * b を m + 1 桁で計算し、最上位からあふれた分を返す
        std::vector<int> r(a.v_.begin() + (n - m + 1), a.v_.end());
        r.push_back(0);
        auto add_mul = [&](long long coef) {
            long long carry = 0;
            for (int j = 0; j <= m; j++) {
                const long long cur =
                    r[j] + carry + (j < m ? coef * b.v_[j] : 0);
                carry = div_floor(cur, (long long)base);
                r[j] = int(cur - carry * base);
            }
            return carry;
        };
        // r >= b か (r は m + 1 桁)
        auto geq_b = [&]() {
            if (r[m] != 0) {
                return true;
            }
            for (int j = m - 1; j >= 0; j--) {
                if (r[j] != b.v_[j]) {
                    return r[j] > b.v_[j];
                }
            }
            return true;
        };

        std::vector<int> q(n - m + 1);
        for (int i = n - m; i >= 0; i--) {
            // r < b を保ったまま、a の桁を 1 つ下ろす
            r.insert(r.begin(), a.v_[i]);
            // r の上位 3 桁 / b の上位 2 桁 で商の桁を見積もる（若干ずれるので後で補正）
            const long long r_top =
                ((long long)r[m] * base + r[m - 1]) * base + r[m - 2];
            int d = int(std::min<long long>(r_top / b_top, base - 1));
            long long over = add_mul(-d);
            while (over < 0) {
                d--;
                over += add_mul(1);
            }
            while (geq_b()) {
                d++;
                add_mul(-1);
            }
            q[i] = d;
            r.pop_back();
        }
        return {from_digits(std::move(q)), from_digits(std::move(r))};
    }

    // b の上位 k 桁を c としたときの base^(2k) / c の近似値（若干ずれるので後で補正）
    // k が b の桁数より大きいときは、b の下に 0 を補って c とする
    static BigInteger reciprocal(const BigInteger& b, int k)
    {
        const BigInteger c = shift(b, k - int(b.v_.size()));
        if (k <= naive_divide_threshold) {
            return divmod_naive(shift(1, 2 * k), c).first;
        }
        // 半分強の精度で求めた逆数から、Newton 法を 1 回まわして精度を倍にする
        const int h = (k + 1) / 2 + 1;
        const BigInteger x = shift(reciprocal(b, h), k - h);
        return x + x - shift(x * x * c, -2 * k);
    }

    // Newton 法。O(M(n))
    //   n: a を base 進数表記したときの桁数
    //   M(n): n 桁同士の掛け算の計算量
    static std::pair<BigInteger, BigInteger> divmod_newton(const BigInteger& a,
                                                           const BigInteger& b)
    {
        const int n = a.v_.size(), m = b.v_.size(), k = n - m + 2;
        BigInteger q = shift(a * reciprocal(b, k), -(k + m));
        BigInteger r = a - q * b;
        while (r.sgn_ < 0) {
            --q;
            r += b;
        }
        while (r >= b) {
            ++q;
            r -= b;
        }
        return {std::move(q), std::move(r)};
    }
};

}  // namespace rklib
