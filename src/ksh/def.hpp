// -*- mode: c++ -*-
#ifndef kashiwa_def_hpp
#define kashiwa_def_hpp
#include <cstddef>
#include <type_traits>
#include <utility>

#ifndef ksh_restrict
# if defined(_MSC_VER) && _MSC_VER >= 1400
#  define ksh_restrict __restrict
# elif defined(__GNUC__)
#  define ksh_restrict __restrict__
# else
#  define ksh_restrict
# endif
#endif

#ifndef ksh_unused
# define ksh_unused(var) ((void) var)
#endif

namespace kashiwa {

  typedef unsigned char byte;

  struct invalid_type {};

  namespace overloads {
    struct adl_inducer {};
  }

  template<typename T> struct identity {typedef T type;};
  template<typename T> using identity_t = typename identity<T>::type;

  //
  // For std::nullptr_t trick (http://qiita.com/kazatsuyu/items/203584ef4cb8b9e52462 by kazatsuyu)
  //
  template<bool B> using nullptr_if_t = typename std::enable_if<B, std::nullptr_t>::type;

  // std::void_t (n3911, C++17) + 変種
  //
  //   std::void_t について思うこと。意味的には void_t は以下の様な定義であるべきではないか。
  //   もしくは態々 void_t を使わずに以下の様に記述するべきだったのではないか。
  //   template<typename... Types> using void_t = std::enable_if_t<std::is_valid_v<Types...>>
  //
  //   ※と思ったが、gcc でコンパイルできない。gcc-6.3.0 以下のバグの様だ。
  //
  template<typename...> using void_t = void;
  template<typename...> using is_valid = std::true_type;
  template<typename...> constexpr bool is_valid_v = true;

  //
  // is_instantiatable (n4502 is_detected 変種)
  //
  namespace detail {
    template<typename Accept, template<typename...> class Template, typename... Args>
    struct instantiater: std::false_type {};
    template<template<typename...> class Template, typename... Args>
    struct instantiater<is_valid<Template<Args...>>, Template, Args...>:
      std::true_type, kashiwa::identity<Template<Args...>> {};
  }
  template<template<typename...> class Template, typename... Args>
  using is_instantiatable = detail::instantiater<std::true_type, Template, Args...>;

  //
  // destructive_negate
  //
  // 型によっては value = -value もしくは value *= -1 を効率よく行うことができる場合がある。
  // 型毎に destructive_negate 関数でその様な操作を提供できる様にするための物。
  // 或る型について destructive_negate を定義する場合は、
  // kashiwa::overloads の下に多重定義を用意する。
  //
  namespace overloads {
    template<typename T> using result_of_destructive_negate = decltype(std::declval<T>().destructive_negate());
    template<typename T> using has_destructive_negate = is_instantiatable<result_of_destructive_negate, T>;

    template<typename K, nullptr_if_t<!has_destructive_negate<K>::value> = nullptr>
    constexpr void destructive_negate(K& value, adl_inducer) {value = -value;}
    template<typename K, nullptr_if_t<has_destructive_negate<K>::value> = nullptr>
    constexpr void destructive_negate(K& value, adl_inducer) {value.destructive_negate();}
  }

  template<typename K>
  constexpr void destructive_negate(K& value) {
    destructive_negate(value, overloads::adl_inducer());
  }

}
#endif
