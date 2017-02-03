// -*- mode: c++ -*-
#ifndef KASHIWA_RATIONAL_H
#define KASHIWA_RATIONAL_H
#include <cstddef>
#include <type_traits>
#include <functional>
#include <ostream>
// #include <vector>
// #include <initializer_list>
// #include <algorithm>
// #include <utility>
namespace kashiwa {

  template<typename K>
  constexpr K gcd(K lhs, K rhs) {
    if (lhs < 0) lhs = -lhs;
    if (rhs < 0) rhs = -rhs;
    switch (lhs > rhs)
      for (;;) {
      case true:
        if (rhs == 0) return lhs;
        lhs %= rhs;
      default:
        if (lhs == 0) return rhs;
        rhs %= lhs;
      }
  }

  struct _canonical_tag {};

  template<typename K>
  struct rational {
    typedef K value_type;
    value_type m_num;
    value_type m_den;

    constexpr rational(K const& num = 0, K const& den = 1): m_num(num), m_den(den) {
      if (m_num == 0 && m_den == 0) return;
      K const _gcd = gcd(m_num, m_den);
      m_num /= _gcd;
      m_den /= _gcd;
      if (m_den < 0) {
        m_num = -m_num;
        m_den = -m_den;
      }
    }

    constexpr rational(K const& num, K const& den, _canonical_tag): m_num(num), m_den(den) {}

    constexpr K const& numerator() const {return m_num;}
    constexpr K const& denominator() const {return m_den;}
    constexpr rational const& operator+() const {return *this;}
    constexpr rational operator-() const {return {-m_num, m_den};}
  };

  template<typename K>
  constexpr bool isfinite(rational<K> const& value) {return value.denominator() != 0;}
  template<typename K>
  constexpr bool isnan(rational<K> const& value) {return value.denominator() == 0 && value.numerator() == 0;}
  template<typename K>
  constexpr bool isinf(rational<K> const& value) {return value.denominator() == 0 && value.numerator() != 0;}

  // rational == rational
  template<typename K>
  constexpr bool operator==(rational<K> const& lhs, rational<K> const& rhs) {
    return lhs.numerator() == rhs.numerator() && lhs.denominator() == rhs.denominator();
  }
  template<typename K>
  constexpr bool operator!=(rational<K> const& lhs, rational<K> const& rhs) {return !(lhs == rhs);}

  // rational == K
  template<typename K>
  constexpr bool operator==(rational<K> const& lhs, K const& rhs) {
    return lhs.numerator() == rhs && lhs.denominator() == 1;
  }
  template<typename K> constexpr bool operator==(K const& lhs, rational<K> const& rhs) {return rhs == lhs;}
  template<typename K> constexpr bool operator!=(rational<K> const& lhs, K const& rhs) {return !(lhs == rhs);}
  template<typename K> constexpr bool operator!=(K const& lhs, rational<K> const& rhs) {return !(rhs == lhs);}

  // rational == int
  namespace rational_detail {
    template<typename K, typename Int>
    using enable_int_overloads_t = typename std::enable_if<
      std::is_same<Int, int>::value && !std::is_same<K, int>::value,
      std::nullptr_t>::type;
  }
  template<typename K, typename Int, rational_detail::enable_int_overloads_t<K, Int> = nullptr>
  constexpr bool operator==(rational<K> const& lhs, Int const& rhs) {
    return lhs.numerator() == rhs && lhs.denominator() == 1;
  }
  template<typename K, typename Int, rational_detail::enable_int_overloads_t<K, Int> = nullptr>
  constexpr bool operator==(Int const& lhs, rational<K> const& rhs) {return rhs == lhs;}
  template<typename K, typename Int, rational_detail::enable_int_overloads_t<K, Int> = nullptr>
  constexpr bool operator!=(rational<K> const& lhs, Int const& rhs) {return !(lhs == rhs);}
  template<typename K, typename Int, rational_detail::enable_int_overloads_t<K, Int> = nullptr>
  constexpr bool operator!=(Int const& lhs, rational<K> const& rhs) {return !(rhs == lhs);}

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

  // rational + rational
  template<typename K, typename Add, typename Neg>
  constexpr rational<K> _impl_add(rational<K> const& lhs, rational<K> const& rhs, Add add, Neg neg) {
    K const& a = lhs.numerator();
    K const& b = rhs.numerator();
    K const& c = lhs.denominator();
    K const& d = rhs.denominator();
    if (c == 0 || d == 0) {
      if (c == d)
        return {a == add(0, b)? a: 0, 0, _canonical_tag()};
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
  template<typename K>
  constexpr rational<K> operator+(rational<K> const& lhs, rational<K> const& rhs) {
    return _impl_add(lhs, rhs, std::plus<K>(), kashiwa::lambda::identity<rational<K>>());
  }
  template<typename K>
  constexpr rational<K> operator-(rational<K> const& lhs, rational<K> const& rhs) {
    return _impl_add(lhs, rhs, std::minus<K>(), std::negate<rational<K>>());
  }
  template<typename K>
  rational<K> operator+=(rational<K>& lhs, rational<K> const& rhs) {return lhs = lhs + rhs;}
  template<typename K>
  rational<K> operator-=(rational<K>& lhs, rational<K> const& rhs) {return lhs = lhs - rhs;}

  // rational * rational
  template<typename K>
  constexpr rational<K> operator*(rational<K> const& lhs, rational<K> const& rhs) {
    K a = lhs.numerator();
    K b = rhs.numerator();
    K c = lhs.denominator();
    K d = rhs.denominator();
    if (c == 0 || d == 0) {
      if (c == d)
        return {a * b, 0, _canonical_tag()}; // {nan, inf} * {nan, inf}
      else if (a == 0 || b == 0)
        return {0, 0, _canonical_tag()}; // 0 * {nan, inf} or nonzero * nan -> nan
      else
        return {(a > 0? 1: -1) * (b > 0? 1: -1), 0, _canonical_tag()}; // nonzero * inf -> inf
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

    return {a * b, c * d, _canonical_tag()};
  }
  template<typename K>
  constexpr rational<K> operator/(rational<K> const& lhs, rational<K> const& rhs) {
    return lhs * rational<K> {rhs.denominator(), rhs.numerator()};
  }
  template<typename K>
  rational<K> operator*=(rational<K>& lhs, rational<K> const& rhs) {return lhs = lhs * rhs;}
  template<typename K>
  rational<K> operator/=(rational<K>& lhs, rational<K> const& rhs) {return lhs = lhs / rhs;}


  // rational * scalar
  template<typename K>
  constexpr rational<K> operator*(rational<K> const& lhs, K const& rhs) {
    rational<K> const tmp {rhs, lhs.denominator()};
    return {lhs.numerator() * tmp.numerator(), tmp.denominator(), _canonical_tag()};
  }
  template<typename K>
  constexpr rational<K> operator/(rational<K> const& lhs, K const& rhs) {
    rational<K> const tmp {lhs.numerator(), rhs};
    return {tmp.numerator(), lhs.denominator() * tmp.denominator(), _canonical_tag()};
  }
  template<typename K>
  constexpr rational<K> operator*(K const& lhs, rational<K> const& rhs) {return rhs * lhs;}
  template<typename K>
  constexpr rational<K> operator/(K const& lhs, rational<K> const& rhs) {
    rational<K> const tmp {lhs, rhs.numerator()};
    return {tmp.numerator() * rhs.denominator(), tmp.denominator(), _canonical_tag()};
  }

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
