#pragma once
// https://github.com/Ricky-pon/rk-library-cpp/blob/main/rklib/number/wide.hpp

#include <concepts>

namespace rklib {

// 2 つの T の積がオーバーフローしない型。より広い型がなければ T 自身
template<class T>
struct Wide {
    using type = T;
};

template<std::signed_integral T>
    requires(sizeof(T) <= 4)
struct Wide<T> {
    using type = long long;
};

template<std::signed_integral T>
    requires(sizeof(T) == 8)
struct Wide<T> {
    using type = __int128;
};

template<class T>
using wide_t = typename Wide<T>::type;

}  // namespace rklib
