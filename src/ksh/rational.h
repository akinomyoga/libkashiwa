// -*- mode: c++ -*-
#ifndef KASHIWA_RATIONAL_H
#define KASHIWA_RATIONAL_H
#include <cstddef>
#include <type_traits>
#include <functional>
#include <ostream>
#include "def.h"
namespace kashiwa {

  namespace lambda {
    template<typename T = void>
    struct identity {
      constexpr T& operator()(T& value) const {return value;}
      constexpr T const& operator()(T const& value) const {return value;}
    };
    template<>
    struct identity<void> {
      template<typename T>
      constexpr T& operator()(T& value) const {return value;}
      template<typename T>
      constexpr T const& operator()(T const& value) const {return value;}
    };
  }

  template<typename K>
  constexpr K gcd(K lhs, K rhs) {
    if (lhs < 0) destructive_negate(lhs);
    if (rhs < 0) destructive_negate(rhs);
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
  constexpr K lcm(K lhs, K rhs) {return (lhs / gcd(lhs, rhs)) * rhs;}

  struct already_canonicalized_tag {};

  template<typename K>
  struct rational {
    typedef K underlying_type;
    underlying_type m_num;
    underlying_type m_den;

    constexpr rational(K const& num = 0, K const& den = 1): m_num(num), m_den(den) {
      if (m_num == 0 && m_den == 0) return;
      K const _gcd = gcd(m_num, m_den);
      m_num /= _gcd;
      m_den /= _gcd;
      if (m_den < 0) {
        kashiwa::destructive_negate(m_num);
        kashiwa::destructive_negate(m_den);
      }
    }

    constexpr rational(K const& num, K const& den, already_canonicalized_tag): m_num(num), m_den(den) {}

    constexpr K const& numerator() const {return m_num;}
    constexpr K const& denominator() const {return m_den;}
    constexpr rational const& operator+() const {return *this;}
    constexpr rational operator-() const {return {-m_num, m_den};}

    constexpr explicit operator bool() const {return m_num != 0;}

    constexpr void destructive_negate() {kashiwa::destructive_negate(m_num);}
  };

  template<typename K>
  constexpr bool isfinite(rational<K> const& value) {return value.denominator() != 0;}
  template<typename K>
  constexpr bool isnan(rational<K> const& value) {return value.denominator() == 0 && value.numerator() == 0;}
  template<typename K>
  constexpr bool isinf(rational<K> const& value) {return value.denominator() == 0 && value.numerator() != 0;}

  //
  // rational == rational
  //
  template<typename K>
  constexpr bool operator==(rational<K> const& lhs, rational<K> const& rhs) {
    return lhs.numerator() == rhs.numerator() && lhs.denominator() == rhs.denominator();
  }
  template<typename K>
  constexpr bool operator!=(rational<K> const& lhs, rational<K> const& rhs) {return !(lhs == rhs);}

  //
  // rational == scalar
  //
  template<typename K, typename L, typename std::enable_if<std::is_convertible<L, K>::value, std::nullptr_t>::type = nullptr>
  constexpr bool operator==(rational<K> const& lhs, L const& rhs) {
    return lhs.numerator() == rhs && lhs.denominator() == 1;
  }
  template<typename K, typename L, typename std::enable_if<std::is_convertible<L, K>::value, std::nullptr_t>::type = nullptr>
  constexpr bool operator==(L const& lhs, rational<K> const& rhs) {return rhs == lhs;}
  template<typename K, typename L, typename std::enable_if<std::is_convertible<L, K>::value, std::nullptr_t>::type = nullptr>
  constexpr bool operator!=(rational<K> const& lhs, L const& rhs) {return !(lhs == rhs);}
  template<typename K, typename L, typename std::enable_if<std::is_convertible<L, K>::value, std::nullptr_t>::type = nullptr>
  constexpr bool operator!=(L const& lhs, rational<K> const& rhs) {return !(rhs == lhs);}

  //
  // rational + rational
  //
  namespace rational_detail {
    template<typename K, typename Add, typename Neg>
    constexpr rational<K> impl_add(rational<K> const& lhs, rational<K> const& rhs, Add add, Neg neg) {
      K const& a = lhs.numerator();
      K const& b = rhs.numerator();
      K const& c = lhs.denominator();
      K const& d = rhs.denominator();
      if (c == 0 || d == 0) {
        if (c == d)
          return {a == add(0, b)? a: 0, 0, already_canonicalized_tag()};
        else
          return c == 0? lhs: neg(rhs);
      }

      K const _gcd = gcd(c, d);
      if (_gcd == 1)
        return {add(a * d, b * c), c * d};
      else {
        K const reducedC = c / _gcd;
        K const reducedD = d / _gcd;
        return {add(a * reducedD, b * reducedC), c * reducedD};
      }
    }
  }
  template<typename K>
  constexpr rational<K> operator+(rational<K> const& lhs, rational<K> const& rhs) {
    return rational_detail::impl_add(lhs, rhs, std::plus<K>(), kashiwa::lambda::identity<rational<K>>());
  }
  template<typename K>
  constexpr rational<K> operator-(rational<K> const& lhs, rational<K> const& rhs) {
    return rational_detail::impl_add(lhs, rhs, std::minus<K>(), std::negate<rational<K>>());
  }
  template<typename K>
  rational<K> operator+=(rational<K>& lhs, rational<K> const& rhs) {return lhs = lhs + rhs;}
  template<typename K>
  rational<K> operator-=(rational<K>& lhs, rational<K> const& rhs) {return lhs = lhs - rhs;}

  //
  // rational + scalar
  //
  template<typename K, typename L, typename std::enable_if<std::is_convertible<L, K>::value, std::nullptr_t>::type = nullptr>
  constexpr rational<K> operator+(rational<K> const& r1, L const& k2) {
    K const& n = r1.numerator();
    K const& d = r1.denominator();
    if (d == 0) return r1;
    return {n + k2 * d, d};
  }
  template<typename K, typename L, typename std::enable_if<std::is_convertible<L, K>::value, std::nullptr_t>::type = nullptr>
  constexpr rational<K> operator-(rational<K> const& r1, L const& k2) {
    K const& n = r1.numerator();
    K const& d = r1.denominator();
    if (d == 0) return r1;
    return {n - k2 * d, d};
  }
  template<typename K, typename L, typename std::enable_if<std::is_convertible<L, K>::value, std::nullptr_t>::type = nullptr>
  constexpr rational<K> operator+(L const& lhs, rational<K> const& rhs) {return rhs + lhs;}
  template<typename K, typename L, typename std::enable_if<std::is_convertible<L, K>::value, std::nullptr_t>::type = nullptr>
  constexpr rational<K> operator-(L const& k1, rational<K> const& r2) {
    K const& n = r2.numerator();
    K const& d = r2.denominator();
    if (d == 0) return -r2;
    return {k1 * d - n, d};
  }
  template<typename K, typename L, typename std::enable_if<std::is_convertible<L, K>::value, std::nullptr_t>::type = nullptr>
  rational<K> operator+=(rational<K>& lhs, L const& rhs) {return lhs = lhs + rhs;}
  template<typename K, typename L, typename std::enable_if<std::is_convertible<L, K>::value, std::nullptr_t>::type = nullptr>
  rational<K> operator-=(rational<K>& lhs, L const& rhs) {return lhs = lhs - rhs;}

  //
  // rational * rational
  //
  template<typename K>
  constexpr rational<K> operator*(rational<K> const& lhs, rational<K> const& rhs) {
    K a = lhs.numerator();
    K b = rhs.numerator();
    K c = lhs.denominator();
    K d = rhs.denominator();
    if (c == 0 || d == 0) {
      if (c == d)
        return {a * b, 0, already_canonicalized_tag()}; // {nan, inf} * {nan, inf}
      else if (a == 0 || b == 0)
        return {0, 0, already_canonicalized_tag()}; // 0 * {nan, inf} or nonzero * nan -> nan
      else
        return {(a > 0? 1: -1) * (b > 0? 1: -1), 0, already_canonicalized_tag()}; // nonzero * inf -> inf
    }

    K const _gcd1 = gcd(a, d);
    if (_gcd1 > 1) {
      a /= _gcd1;
      d /= _gcd1;
    }
    K const _gcd2 = gcd(b, c);
    if (_gcd2 > 1) {
      b /= _gcd2;
      c /= _gcd2;
    }

    return {a * b, c * d, already_canonicalized_tag()};
  }
  template<typename K>
  constexpr rational<K> operator/(rational<K> const& lhs, rational<K> const& rhs) {
    return lhs * rational<K> {rhs.denominator(), rhs.numerator()};
  }
  template<typename K>
  rational<K> operator*=(rational<K>& lhs, rational<K> const& rhs) {return lhs = lhs * rhs;}
  template<typename K>
  rational<K> operator/=(rational<K>& lhs, rational<K> const& rhs) {return lhs = lhs / rhs;}

  //
  // rational * scalar
  //
  template<typename K, typename L, typename std::enable_if<std::is_convertible<L, K>::value, std::nullptr_t>::type = nullptr>
  constexpr rational<K> operator*(rational<K> const& lhs, L const& rhs) {
    rational<K> const tmp {rhs, lhs.denominator()};
    return {lhs.numerator() * tmp.numerator(), tmp.denominator(), already_canonicalized_tag()};
  }
  template<typename K, typename L, typename std::enable_if<std::is_convertible<L, K>::value, std::nullptr_t>::type = nullptr>
  constexpr rational<K> operator/(rational<K> const& lhs, L const& rhs) {
    rational<K> const tmp {lhs.numerator(), rhs};
    return {tmp.numerator(), lhs.denominator() * tmp.denominator(), already_canonicalized_tag()};
  }
  template<typename K, typename L, typename std::enable_if<std::is_convertible<L, K>::value, std::nullptr_t>::type = nullptr>
  constexpr rational<K> operator*(L const& lhs, rational<K> const& rhs) {return rhs * lhs;}
  template<typename K, typename L, typename std::enable_if<std::is_convertible<L, K>::value, std::nullptr_t>::type = nullptr>
  constexpr rational<K> operator/(L const& lhs, rational<K> const& rhs) {
    rational<K> const tmp {lhs, rhs.numerator()};
    return {tmp.numerator() * rhs.denominator(), tmp.denominator(), already_canonicalized_tag()};
  }
  template<typename K, typename L, typename std::enable_if<std::is_convertible<L, K>::value, std::nullptr_t>::type = nullptr>
  rational<K> operator*=(rational<K>& lhs, L const& rhs) {return lhs = lhs * rhs;}
  template<typename K, typename L, typename std::enable_if<std::is_convertible<L, K>::value, std::nullptr_t>::type = nullptr>
  rational<K> operator/=(rational<K>& lhs, L const& rhs) {return lhs = lhs / rhs;}

  //
  // rational < rational
  //
  template<typename K> constexpr bool operator<(rational<K> const& lhs, rational<K> const& rhs) {return (lhs - rhs).numerator() < 0;}
  template<typename K> constexpr bool operator<=(rational<K> const& lhs, rational<K> const& rhs) {return (lhs - rhs).numerator() < 0;}
  template<typename K> constexpr bool operator>(rational<K> const& lhs, rational<K> const& rhs) {return (lhs - rhs).numerator() > 0;}
  template<typename K> constexpr bool operator>=(rational<K> const& lhs, rational<K> const& rhs) {return (lhs - rhs).numerator() >= 0;}

  //
  // rational < scalar
  //
  namespace rational_detail {
    template<typename K, typename L, typename Compare>
    constexpr bool impl_compare(rational<K> const& lhs, L const& rhs, Compare compare, int inf) {
      if (lhs.denominator() == 0)
        return lhs.numerator() == inf;
      else
        return compare(lhs.numerator(), lhs.denominator() * rhs); // overflow?
    }
  }
  template<typename K, typename L, typename std::enable_if<std::is_convertible<L, K>::value, std::nullptr_t>::type = nullptr>
  constexpr bool operator<(rational<K> const& lhs, L const& rhs) {
    return rational_detail::impl_compare(lhs, rhs, std::less<K>(), -1);
  }
  template<typename K, typename L, typename std::enable_if<std::is_convertible<L, K>::value, std::nullptr_t>::type = nullptr>
  constexpr bool operator<=(rational<K> const& lhs, L const& rhs) {
    return rational_detail::impl_compare(lhs, rhs, std::less_equal<K>(), -1);
  }
  template<typename K, typename L, typename std::enable_if<std::is_convertible<L, K>::value, std::nullptr_t>::type = nullptr>
  constexpr bool operator>(rational<K> const& lhs, L const& rhs) {
    return rational_detail::impl_compare(lhs, rhs, std::greater<K>(), 1);
  }
  template<typename K, typename L, typename std::enable_if<std::is_convertible<L, K>::value, std::nullptr_t>::type = nullptr>
  constexpr bool operator>=(rational<K> const& lhs, L const& rhs) {
    return rational_detail::impl_compare(lhs, rhs, std::greater_equal<K>(), 1);
  }
  template<typename K, typename L, typename std::enable_if<std::is_convertible<L, K>::value, std::nullptr_t>::type = nullptr>
  constexpr bool operator< (L const& lhs, rational<K> const& rhs) {return rhs > lhs;}
  template<typename K, typename L, typename std::enable_if<std::is_convertible<L, K>::value, std::nullptr_t>::type = nullptr>
  constexpr bool operator> (L const& lhs, rational<K> const& rhs) {return rhs < lhs;}
  template<typename K, typename L, typename std::enable_if<std::is_convertible<L, K>::value, std::nullptr_t>::type = nullptr>
  constexpr bool operator<=(L const& lhs, rational<K> const& rhs) {return rhs >= lhs;}
  template<typename K, typename L, typename std::enable_if<std::is_convertible<L, K>::value, std::nullptr_t>::type = nullptr>
  constexpr bool operator>=(L const& lhs, rational<K> const& rhs) {return rhs <= lhs;}

  //
  // ostr << rational
  //
  template<typename K>
  std::ostream& operator<<(std::ostream& ostr, rational<K> const& value) {
    if (value.denominator() == 0)
      return ostr << (value.numerator() == 0? "nan": value.numerator() > 0? "inf": "-inf");
    else if (value.denominator() == 1)
      return ostr << value.numerator();
    else
      return ostr << '(' << value.numerator() << '/' << value.denominator() << ')';
  }

}
#endif
