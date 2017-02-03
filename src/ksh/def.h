// -*- mode: c++ -*-
#ifndef KASHIWA_DEF_H
#define KASHIWA_DEF_H
#include <cstddef>
#include <type_traits>
#include <utility>
namespace kashiwa {
  typedef unsigned char byte;

  struct invalid_type {};

  namespace overloads {
    struct adl_inducer {};
  }

  //
  // destructive_negate
  //
  // 型によっては value = -value もしくは value *= -1 を効率よく行うことができる場合がある。
  // 型毎に destructive_negate 関数でその様な操作を提供できる様にするための物。
  // 或る型について destructive_negate を定義する場合は、
  // kashiwa::overloads の下に多重定義を用意する。
  //
  namespace overloads {
#define kashiwa_define_is_valid_expression(Name,T,X,Expr) \
    struct Name { \
      using invalid_type = ::kashiwa::invalid_type; \
      template<typename X> static constexpr invalid_type check(...); \
      template<typename X> static constexpr auto check(int) -> decltype((Expr)); \
      using return_type = decltype(check<T>(0)); \
      enum {value = !std::is_same<return_type, invalid_type>::value}; \
    };

    template<typename T>
    kashiwa_define_is_valid_expression(has_destructive_negate, T, X, (std::declval<X>().destructive_negate()));

    template<typename K, typename std::enable_if<!has_destructive_negate<K>::value, std::nullptr_t>::type = nullptr>
    constexpr void destructive_negate(K& value, adl_inducer) {value = -value;}
    template<typename K, typename std::enable_if<has_destructive_negate<K>::value, std::nullptr_t>::type = nullptr>
    constexpr void destructive_negate(K& value, adl_inducer) {value.destructive_negate();}
  }

  template<typename K>
  constexpr void destructive_negate(K& value) {
    destructive_negate(value, overloads::adl_inducer());
  }

}
#endif
