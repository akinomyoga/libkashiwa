// -*- mode: c++ -*-
#ifndef kashiwa_arithmetic_hpp
#define kashiwa_arithmetic_hpp
#include <cmath>
#include <type_traits>
#include <mwg/except.h>
#include "def.hpp"

// This header provides customization points for basic arithmetic operations
// for general algebraic types, including the built-in integral and
// floating-point numbers as well as the program-defined types.  The operations
// associated with the C++ operators are not covered because their
// customization points are natively provided by operator overloading.

namespace kashiwa {

  template<typename K>
  constexpr K inv(K const& x) { return 1 / x; }

  template<typename K>
  constexpr void chinv(K& x) { x = inv(x); }

  // chneg (destructive negate)
  //
  // 型によっては value = -value もしくは value *= -1 を効率よく行うことができる場合がある。
  // 型毎に chneg 関数でその様な操作を提供できる様にするための物。
  // 或る型について chneg を定義する場合は、
  // kashiwa::overloads の下に多重定義を用意する。
  //
  // Note: neg() is not provided because it corresponds to the existing
  // operator-().
  namespace overloads {
    template<typename T> using result_of_chneg = decltype(std::declval<T>().chneg());
    template<typename T> using has_chneg = is_instantiatable<result_of_chneg, T>;

    template<typename K, nullptr_if_t<!has_chneg<K>::value> = nullptr>
    constexpr void chneg(K& value, adl_inducer) { value = -value; }
    template<typename K, nullptr_if_t<has_chneg<K>::value> = nullptr>
    constexpr void chneg(K& value, adl_inducer) { value.chneg(); }
  }

  template<typename K>
  constexpr void chneg(K& value) {
    chneg(value, overloads::adl_inducer());
  }

  template<typename K, typename M>
  constexpr M mod(K const& x, M const& m) {
    M ret = x % m;
    if (ret < 0) ret += m;
    return ret;
  }

  template<typename K, typename M>
  void chmod(K& x, M const& m) {
    x %= m;
    if (x < 0) x += m;
  }

  template<typename K>
  K div(K const& x, K const& y, K& q) {
    q = x / y;
    return x % y;
  }

  int div(int x, int y, int& q) {
    std::div_t d = std::div(x, y);
    q = d.quot;
    return d.rem;
  }
  long div(long x, long y, long& q) {
    std::ldiv_t d = std::ldiv(x, y);
    q = d.quot;
    return d.rem;
  }
  long long div(long long x, long long y, long long& q) {
    std::lldiv_t d = std::lldiv(x, y);
    q = d.quot;
    return d.rem;
  }

  template<typename K>
  constexpr K gcd(K lhs, K rhs) {
    if (lhs < 0) chneg(lhs);
    if (rhs < 0) chneg(rhs);
#ifdef __clang__
# pragma clang diagnostic push
# pragma clang diagnostic ignored "-Wswitch-bool"
#endif
    switch (lhs > rhs)
      for (;;) {
      case true:
        if (rhs == 0) return lhs;
        lhs %= rhs;
      default:
        if (lhs == 0) return rhs;
        rhs %= lhs;
      }
#ifdef __clang__
# pragma clang diagnostic pop
#endif
  }

  template<typename K>
  constexpr K lcm(K const& lhs, K const& rhs) {
    return (lhs / gcd(lhs, rhs)) * rhs;
  }

  template<typename K, typename I, typename OpMulEq>
  constexpr K ipow(K base, I exp, OpMulEq chmul) {
    // mwg_assert(exp >= 0); // This is not "constexpr".
    if (exp <= 0) return K(1);

    while (exp % 2 == 0) {
      chmul(base, base);
      exp /= 2;
    }

    K ret { base };
    while ((exp /= 2) > 0) {
      chmul(base, base);
      if (exp % 2) chmul(ret, base);
    }

    return ret;
  }

  template<typename K, typename I, nullptr_if_t<std::is_integral<I>::value> = nullptr>
  constexpr K ipow(K base, I exp) {
    return ipow(base, exp, [] (K& u, K const& v) { u *= v; });
  }

  template<typename K, typename I>
  constexpr K ipow_mod(K base, I exp, K const& mod) {
    return ipow(base, exp, [&mod] (K& u, K const& v) { u *= v; u %= mod; });
  }

  /// @return This function returns the inverse of the VALUE under modulo MOD.
  /// When the return value is 0, it implies that the inverse of VALUE does not
  /// exist.
  template<typename K, typename M>
  constexpr M inv_mod(K const& value, M const& mod) {
    M y = mod, b = 0;
    if (y < 0) kashiwa::chneg(y);
    M x = kashiwa::mod(value, y), a = 1;

    M q;
    for (;;) {
      if (x == 0) {
        if (y != 1) return M(0);
        kashiwa::chmod(b, mod);
        return b;
      }
      y = kashiwa::div(y, x, q);
      b -= q * a;

      if (y == 0) {
        if (x != 1) return M(0);
        kashiwa::chmod(a, mod);
        return a;
      }
      x = kashiwa::div(x, y, q);
      a -= q * b;
    }
  }

}
#endif
