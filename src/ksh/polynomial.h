// -*- mode: c++ -*-
#ifndef KASHIWA_POLYNOMIAL_H
#define KASHIWA_POLYNOMIAL_H
#include <cstddef>
#include <vector>
#include <initializer_list>
#include <algorithm>
#include <utility>
#include <type_traits>
#include <ostream>
namespace kashiwa {

  // K に対する要求
  //
  // - K は int から初期化できる。
  // - 特に 0 は零元に対応する。
  // - K は int と直接比較できる。
  //
  template<typename K>
  struct polynomial {
    std::vector<K> m_data;

  public:
    polynomial() {}
    polynomial(std::initializer_list<K> list): m_data(list) {
      this->normalize();
    }

    std::vector<K> const& data() const {return this->m_data;}

    K assign(K const& value) const {
      K result = 0;
      for (std::size_t i = m_data.size(); i--;)
        result = result * value + m_data[i];
      return result;
    }

    polynomial& operator+=(polynomial<K> const& rhs) {
      if (m_data.size() < rhs.m_data.size())
        m_data.resize(rhs.m_data.size());
      for (std::size_t i = 0; i < rhs.m_data.size(); i++)
        m_data[i] += rhs.m_data[i];
      return *this;
    }
    polynomial& operator-=(polynomial<K> const& rhs) {
      if (m_data.size() < rhs.m_data.size())
        m_data.resize(rhs.m_data.size());
      for (std::size_t i = 0; i < rhs.m_data.size(); i++)
        m_data[i] -= rhs.m_data[i];
      return *this;
    }

    void normalize() {
      std::size_t i = this->m_data.size();
      if (this->m_data[i - 1] != 0) return;
      for (i--; this->m_data[i - 1] == 0; i--);
      this->m_data.resize(i);
    }

    polynomial const& operator+() const {return *this;}

    polynomial operator-() const {
      polynomial ret(*this);
      for (auto& e: ret.m_data) e *= -1;
      return ret;
    }
  };

  template<typename K>
  std::size_t deg(polynomial<K> const& poly) {return std::max(poly.data().size(), (std::size_t) 1) - 1;}

  template<typename K>
  bool operator==(polynomial<K> const& lhs, polynomial<K> const& rhs) {
    return lhs.data().size() == rhs.data().size()
      && std::equal(lhs.data().begin(), lhs.data().end(), rhs.data().begin());
  }
  template<typename K>
  bool operator!=(polynomial<K> const& lhs, polynomial<K> const& rhs) {return !(lhs == rhs);}

  // overloads compare to K values
  template<typename K>
  bool operator==(polynomial<K> const& lhs, K const& rhs) {
    std::vector<K> const& data = lhs.data();
    switch (data.size()) {
    case 0: return rhs == 0;
    case 1: return rhs == data[0];
    default: return false;
    }
  }
  template<typename K> bool operator==(K const& lhs, polynomial<K> const& rhs) {return rhs == lhs;}
  template<typename K> bool operator!=(polynomial<K> const& lhs, K const& rhs) {return !(lhs == rhs);}
  template<typename K> bool operator!=(K const& lhs, polynomial<K> const& rhs) {return !(rhs == lhs);}

  // overloads compare to int literals
  namespace polynomial_detail {
    template<typename K, typename Int>
    using enable_int_overloads_t = typename std::enable_if<
      std::is_same<Int, int>::value && !std::is_same<K, int>::value,
      std::nullptr_t>::type;
  }
  template<typename K, typename Int, polynomial_detail::enable_int_overloads_t<K, Int> = nullptr>
  bool operator==(polynomial<K> const& lhs, Int const& rhs) {
    std::vector<K> const& data = lhs.data();
    switch (data.size()) {
    case 0: return rhs == 0;
    case 1: return rhs == data[0];
    default: return false;
    }
  }
  template<typename K, typename Int, polynomial_detail::enable_int_overloads_t<K, Int> = nullptr>
  bool operator==(Int const& lhs, polynomial<K> const& rhs) {return rhs == lhs;}
  template<typename K, typename Int, polynomial_detail::enable_int_overloads_t<K, Int> = nullptr>
  bool operator!=(polynomial<K> const& lhs, Int const& rhs) {return !(lhs == rhs);}
  template<typename K, typename Int, polynomial_detail::enable_int_overloads_t<K, Int> = nullptr>
  bool operator!=(Int const& lhs, polynomial<K> const& rhs) {return !(rhs == lhs);}

  template<typename K, typename F>
  polynomial<K> merge(polynomial<K> const& lhs, polynomial<K> const& rhs, F f) {
    std::vector<K> const& ldata = lhs.data();
    std::vector<K> const& rdata = rhs.data();
    K const* pL = &ldata[0];
    K const* pR = &rdata[0];

    K const* p1;
    std::size_t iN1 = ldata.size();
    std::size_t iN2 = rdata.size();
    if (iN1 >= iN2) {
      p1 = pL;
    } else {
      p1 = pR;
      std::swap(iN1, iN2);
    }

    polynomial<K> result;
    result.m_data.reserve(iN1);
    std::size_t i = 0;
    for (; i < iN2; i++) result.m_data.emplace_back(f(pL[i], pR[i]));
    for (; i < iN1; i++) result.m_data.emplace_back(p1[i]);
    result.normalize();
    return result;
  }

  template<typename K>
  polynomial<K> operator+(polynomial<K> const& lhs, polynomial<K> const& rhs) {
    return merge(lhs, rhs, [](K const& lhs, K const& rhs) {return lhs + rhs;});
  }
  template<typename K>
  polynomial<K> operator-(polynomial<K> const& lhs, polynomial<K> const& rhs) {
    return merge(lhs, rhs, [](K const& lhs, K const& rhs) {return lhs - rhs;});
  }

  //
  // polynomial * polynomial
  //
  template<typename K>
  polynomial<K> operator*(polynomial<K> const& lhs, polynomial<K> const& rhs) {
    polynomial<K> result;
    if (lhs.m_data.size() == 0 || rhs.m_data.size() == 0) return result;

    std::vector<K>& data = result.m_data;
    std::vector<K> const& ldata = lhs.m_data;
    std::vector<K> const& rdata = rhs.m_data;

    std::size_t const iN = ldata.size();
    std::size_t const jN = rdata.size();
    data.resize(iN + jN - 1, (K) 0);

    for (std::size_t i = 0; i < iN; i++)
      for (std::size_t j = 0; j < jN; j++)
        data[i + j] += ldata[i] * rdata[j];

    return result;
  }
  template<typename K>
  polynomial<K>& operator*=(polynomial<K>& lhs, polynomial<K> const& rhs) {return lhs = lhs * rhs;}

  //
  // polynomial * scalar
  //
  template<typename K>
  polynomial<K>& operator*=(polynomial<K>& lhs, K const& rhs) {
    if (rhs == 0)
      lhs.m_data.clear();
    else
      for (auto& e: lhs.m_data) e *= rhs;
    return lhs;
  }
  template<typename K>
  polynomial<K> operator*(polynomial<K> const& lhs, K const& rhs) {
    polynomial<K> ret(lhs);
    ret *= rhs;
    return ret;
  }
  template<typename K>
  polynomial<K> operator*(K const& lhs, polynomial<K> const& rhs) {return rhs * lhs;}

  template<typename K>
  polynomial<K> pow(polynomial<K> const& lhs, unsigned exponent) {
    polynomial<K> result {1};
    if (!exponent) return result;

    polynomial<K> pow2 = lhs;
    for (;;) {
      unsigned const digit = exponent & 1;
      exponent >>= 1;
      if (digit) {
        result *= pow2;
        if (!exponent) break;
      }
      pow2 *= pow2;
    }

    return result;
  }

  // template<typename Ch, typename Tr, typename K>
  // std::basic_ostream<Ch, Tr>& print(std::basic_ostream& ostr, std::polynomial<K> const& poly) {
  //   std::ostream << (Ch) 'x'
  // }

  template<typename K>
  std::ostream& operator<<(std::ostream& ostr, polynomial<K> const& poly) {
    auto const& data = poly.data();
    if (data.size() == 0)
      return ostr << '0';
    else {
      std::size_t i = data.size() - 1;

      {
        K const& coefficient = data[i];
        if (i == 0)
          return ostr << coefficient;
        else {
          if (coefficient == -1)
            ostr << '-';
          else if (coefficient != 1)
            ostr << coefficient;
          ostr << 'x';
          if (i > 1) ostr << '^' << i;
        }
      }

      do {
        if (K coefficient = data[--i]) {
          ostr << ' ';
          if (coefficient < 0) {
            ostr << '-';
            coefficient = -coefficient;
          } else
            ostr << '+';
          ostr << ' ';

          if (i == 0)
            return ostr << coefficient;
          else {
            if (coefficient != 1)
              ostr << coefficient;
            ostr << 'x';
            if (i > 1) ostr << '^' << i;
          }
        }
      } while (i);
      return ostr;
    }
  }

}
#endif
