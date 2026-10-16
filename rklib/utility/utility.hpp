#pragma once
// https://github.com/Ricky-pon/rk-library-cpp/blob/main/rklib/utility/utility.hpp

namespace rklib {

// a <- max(a, b)
template<class T, class U>
    requires requires(U u) { T{u}; }
bool chmax(T& a, const U& b)
{
    if (a < b) {
        a = b;
        return true;
    }
    return false;
}

// a <- min(a, b)
template<class T, class U>
    requires requires(U u) { T{u}; }
bool chmin(T& a, const U& b)
{
    if (a > b) {
        a = b;
        return true;
    }
    return false;
}

// floor(num / den)
template<class T>
T div_floor(T num, T den)
{
    if (den < 0) {
        num = -num, den = -den;
    }
    return num >= 0 ? num / den : (num + 1) / den - 1;
}

// ceil(num / den)
template<class T>
T div_ceil(T num, T den)
{
    if (den < 0) {
        num = -num, den = -den;
    }
    return num <= 0 ? num / den : (num - 1) / den + 1;
}

}  // namespace rklib
