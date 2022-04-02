// -*- C++ -*-
#ifndef kashiwa_biginteger_hpp
#define kashiwa_biginteger_hpp
#include <cstddef>
#include <cstdint>
#include <limits>
#include <type_traits>
#include <vector>
#include <ostream>
#include <algorithm>
#include <mwg/except.h>
#include "def.hpp"
namespace kashiwa {

  template<
    typename ElemInt = std::uint32_t,
    typename CalcInt = std::uint64_t,
    CalcInt modulo = (CalcInt) std::numeric_limits<ElemInt>::max() + 1 >
  struct big_integer;

  template<typename ElemInt, typename CalcInt, CalcInt Modulo>
  struct big_integer {
    typedef ElemInt element_type;
    typedef CalcInt calculation_type;
    static constexpr element_type element_max = Modulo - 1;
    static constexpr calculation_type modulo = Modulo;

    static_assert(
      std::is_unsigned<element_type>::value &&
      std::is_unsigned<calculation_type>::value &&
      Modulo > 1 &&
      Modulo - 1 <= std::numeric_limits<element_type>::max() &&
      element_max <= std::numeric_limits<calculation_type>::max() / element_max);

    int sign;
    std::vector<element_type> data;

    big_integer(): sign(0) {}

    template<typename I, nullptr_if_t<std::is_integral<I>::value> = nullptr>
    big_integer(I const& value) {this->operator=(value);}

    big_integer const& operator+() const {return *this;}
    big_integer operator-() const {
      big_integer ret {*this};
      ret.sign = -ret.sign;
      return ret;
    }

  private:
    void abs_dec() {
      for (std::size_t i = 0; i < data.size(); i++) {
        if (data[i]) {
          data[i]--;
          if (i == data.size() - 1 && data[i] == 0) {
            data.pop_back();
            if (data.size() == 0) sign = 0;
          }
          return;
        } else {
          data[i] = element_max;
        }
      }
    }
    void abs_inc() {
      for (std::size_t i = 0; i < data.size(); i++) {
        if (data[i] < element_max) {
          data[i]++;
          return;
        } else {
          data[i] = 0;
        }
      }
      data.push_back(1);
    }

  public:
    big_integer& operator++() {
      if (sign == -1)
        abs_dec();
      else if (sign == 1)
        abs_inc();
      else {
        sign = 1;
        data.push_back(1);
      }
      return *this;
    }

    big_integer& operator--() {
      if (sign == -1)
        abs_inc();
      else if (sign == 1)
        abs_dec();
      else {
        sign = -1;
        data.push_back(1);
      }
      return *this;
    }

    big_integer operator++(int) {
      big_integer ret {*this};
      this->operator++();
      return ret;
    }

    big_integer operator--(int) {
      big_integer ret {*this};
      this->operator--();
      return ret;
    }

  public:
    template<typename I, nullptr_if_t<std::is_integral<I>::value> = nullptr>
    big_integer& operator=(I const& value) {
      data.clear();
      if (value == 0) {
        sign = 0;
        return *this;
      }

      typename std::make_unsigned<I>::type uvalue = value;
      if (std::is_signed<I>::value && value < 0) {
        sign = -1;
        uvalue = -uvalue;
      } else {
        sign = 1;
      }

      while (uvalue) {
        data.emplace_back((element_type) (uvalue % modulo));
        uvalue /= modulo;
      }
      return *this;
    }
  };

  namespace big_integer_detail {
    template<typename MpInteger, typename T>
    using enable_generic_operator_t = nullptr_if_t<(std::is_integral<T>::value || std::is_base_of<MpInteger, T>::value)>;
    template<typename T>
    using enable_scalar_operator_t = nullptr_if_t<std::is_integral<T>::value>;

    template<typename E, typename C, C M>
    int _sign(big_integer<E, C, M> const&  arg) {return arg.sign;}
    template<typename I, enable_scalar_operator_t<I> = nullptr>
    int _sign(I const&  arg) {
      if (arg > 0)
        return 1;
      else if (std::is_signed<I>::value && arg < 0)
        return -1;
      else
        return 0;
    }
  }

  //
  // compare(a, b), a == b, a != b, a < b, etc.
  //
  namespace big_integer_detail {
    template<typename I, nullptr_if_t<std::is_integral<I>::value && std::is_signed<I>::value> = nullptr>
    typename std::make_unsigned<I>::type _abs(I const& value) {
      typename std::make_unsigned<I>::type uvalue = value;
      if (value < 0) uvalue = -uvalue;
      return uvalue;
    }
    template<typename U, nullptr_if_t<std::is_integral<U>::value && std::is_unsigned<U>::value> = nullptr>
    U const& _abs(U const& value) {return value;}

    template<typename E, typename C, C M>
    int abs_compare(big_integer<E, C, M> const& lhs, big_integer<E, C, M> const& rhs) {
      if (lhs.data.size() != rhs.data.size())
        return lhs.data.size() < rhs.data.size()? -1: 1;

      for (std::size_t i = lhs.data.size(); i--; )
        if (lhs.data[i] != rhs.data[i])
          return lhs.data[i] < rhs.data[i]? -1: 1;

      return 0;
    }
    // lhs と rhs * M^exp を比較
    template<typename E, typename C, C M>
    int abs_compare(big_integer<E, C, M> const& lhs, big_integer<E, C, M> const& rhs, std::size_t exp) {
      if (lhs.data.size() != rhs.data.size() + exp)
        return lhs.data.size() < rhs.data.size() + exp? -1: 1;

      for (std::size_t i = rhs.data.size(); i--; )
        if (lhs.data[i + exp] != rhs.data[i])
          return lhs.data[i + exp] < rhs.data[i]? -1: 1;

      return 0;
    }

    template<typename E, typename C, C M, typename I, enable_scalar_operator_t<I> = nullptr>
    int abs_compare(big_integer<E, C, M> const& lhs, I const& rhs) {
      using integer_t = big_integer<E, C, M>;
      using elem_t = typename integer_t::element_type;

      typename std::make_unsigned<I>::type urhs = _abs(rhs);

      int ret = 0;
      for (std::size_t i = 0; i < lhs.data.size(); i++) {
        if (urhs == 0) return 1;
        elem_t const part = (elem_t) (urhs % integer_t::modulo);
        if (lhs.data[i] != part)
          ret = lhs.data[i] > part? 1: -1;
        urhs /= integer_t::modulo;
      }
      return urhs == 0? ret: -1;
    }

    template<typename E, typename C, C M, typename T>
    int impl_compare(big_integer<E, C, M> const& lhs, T const& rhs) {
      int const rsign = _sign(rhs);
      if (rsign > 0) {
        if (lhs.sign <= 0) return -1;
      } else if (rsign < 0) {
        if (lhs.sign >= 0) return 1;
      } else
        return lhs.sign;
      return rsign * abs_compare(lhs, rhs);
    }
    template<typename E, typename C, C M>
    int compare(big_integer<E, C, M> const& lhs, big_integer<E, C, M> const& rhs) {return impl_compare(lhs, rhs);}
    template<typename E, typename C, C M, typename I, enable_scalar_operator_t<I> = nullptr>
    int compare(big_integer<E, C, M> const& lhs, I const& rhs) {return impl_compare(lhs, rhs);}
    template<typename E, typename C, C M, typename I, enable_scalar_operator_t<I> = nullptr>
    int compare(I const& lhs, big_integer<E, C, M> const& rhs) {return -impl_compare(rhs, lhs);}

    template<typename LHS, typename RHS, nullptr_if_t<is_valid_v<decltype(compare(std::declval<LHS>(), std::declval<RHS>()))>> = nullptr>
    bool operator==(LHS const& lhs, RHS const& rhs) {return compare(lhs, rhs) == 0;}
    template<typename LHS, typename RHS, nullptr_if_t<is_valid_v<decltype(compare(std::declval<LHS>(), std::declval<RHS>()))>> = nullptr>
    bool operator!=(LHS const& lhs, RHS const& rhs) {return compare(lhs, rhs) != 0;}
    template<typename LHS, typename RHS, nullptr_if_t<is_valid_v<decltype(compare(std::declval<LHS>(), std::declval<RHS>()))>> = nullptr>
    bool operator< (LHS const& lhs, RHS const& rhs) {return compare(lhs, rhs) <  0;}
    template<typename LHS, typename RHS, nullptr_if_t<is_valid_v<decltype(compare(std::declval<LHS>(), std::declval<RHS>()))>> = nullptr>
    bool operator> (LHS const& lhs, RHS const& rhs) {return compare(lhs, rhs) >  0;}
    template<typename LHS, typename RHS, nullptr_if_t<is_valid_v<decltype(compare(std::declval<LHS>(), std::declval<RHS>()))>> = nullptr>
    bool operator<=(LHS const& lhs, RHS const& rhs) {return compare(lhs, rhs) <= 0;}
    template<typename LHS, typename RHS, nullptr_if_t<is_valid_v<decltype(compare(std::declval<LHS>(), std::declval<RHS>()))>> = nullptr>
    bool operator>=(LHS const& lhs, RHS const& rhs) {return compare(lhs, rhs) >= 0;}
  }
  using big_integer_detail::compare;
  using big_integer_detail::operator==;
  using big_integer_detail::operator!=;
  using big_integer_detail::operator> ;
  using big_integer_detail::operator< ;
  using big_integer_detail::operator>=;
  using big_integer_detail::operator<=;

  //
  // a += b, a -= b, a + b, a - b
  //
  namespace big_integer_detail {

    template<typename E, typename C, C M>
    void abs_add(big_integer<E, C, M>& lhs, big_integer<E, C, M> const& rhs) {
      typedef big_integer<E, C, M> integer_t;
      typedef typename integer_t::element_type  element_t;
      typedef typename integer_t::calculation_type  calc_t;

      if (lhs.data.size() < rhs.data.size())
        lhs.data.resize(rhs.data.size(), (element_t) 0u);

      calc_t carry = 0;
      std::size_t const iN = rhs.data.size();
      for (std::size_t i = 0; i < iN; i++) {
        calc_t const value = (calc_t) lhs.data[i] + (calc_t) rhs.data[i] + carry;
        carry       = value / integer_t::modulo;
        lhs.data[i] = (element_t) (value % integer_t::modulo);
      }

      if (carry) {
        mwg_assert(iN <= lhs.data.size());
        if (iN < lhs.data.size())
          lhs.data[iN] = carry;
        else
          lhs.data.push_back(carry);
      }
    }
    template<typename E, typename C, C M, typename I, enable_scalar_operator_t<I> = nullptr>
    void abs_add(big_integer<E, C, M>& lhs, I const& rhs) {
      using integer_t = big_integer<E, C, M>;
      using elem_t = typename integer_t::element_type;
      using calc_t = typename integer_t::calculation_type;

      typename std::make_unsigned<I>::type urhs = _abs(rhs);

      for (std::size_t i = 0; i < lhs.data.size(); i++) {
        if (urhs == 0) return;
        calc_t const sub = (calc_t) (urhs % integer_t::modulo) + (calc_t) lhs.data[i];
        urhs /= integer_t::modulo;
        if (sub >= integer_t::modulo) {
          urhs++;
          lhs.data[i] = (elem_t) (sub % integer_t::modulo);
        } else {
          lhs.data[i] = (elem_t) sub;
        }
      }

      while (urhs) {
        lhs.data.emplace_back((elem_t) (urhs % integer_t::modulo));
        urhs /= integer_t::modulo;
      }
    }

    template<typename E, typename C, C M>
    void abs_sub(big_integer<E, C, M>& lhs, big_integer<E, C, M> const& rhs) {
      typedef big_integer<E, C, M> integer_t;
      typedef typename integer_t::element_type elem_t;
      typedef typename integer_t::calculation_type calc_t;

      int const cmp = abs_compare(lhs, rhs);
      if (cmp == 0) {
        lhs.sign = 0;
        lhs.data.resize(1, (elem_t) 0);
      } else if (cmp != 0) {
        elem_t carry = 0;
        elem_t      * dst = &lhs.data[0];
        elem_t const* max;
        elem_t const* min;
        std::size_t maxN;
        std::size_t minN;

        if (cmp < 0) {
          lhs.sign = -lhs.sign;
          lhs.data.resize(rhs.data.size(), (elem_t) 0);
          max = &rhs.data[0];
          min = &lhs.data[0];
          maxN = rhs.data.size();
          minN = lhs.data.size();
        } else {
          max = &lhs.data[0];
          min = &rhs.data[0];
          maxN = lhs.data.size();
          minN = rhs.data.size();
        }

        for (std::size_t i = 0; i < minN; i++) {
          calc_t const sub = (calc_t) min[i] + carry;
          if (max[i] < sub) {
            dst[i] = (elem_t) (max[i] + (integer_t::modulo - sub));
            carry = 1;
          } else {
            dst[i] = (elem_t) (max[i] - sub);
            carry = 0;
          }
        }

        std::size_t highest = minN - 1;
        if (dst != max || carry) {
          highest = maxN - 1;
          for (std::size_t j = minN; j < maxN; j++) {
            if (max[j] < carry) {
              mwg_assert(max[j] == 0 && carry == 1);
              dst[j] = integer_t::element_max;
            } else {
              dst[j] = max[j] - carry;
              if (dst == max) {highest = j; break;}
              carry = 0;
            }
          }
        }

        if (highest + 1 == lhs.data.size()) {
          while (highest > 0 && lhs.data[highest] == 0) highest--;
          lhs.data.resize(highest + 1, (elem_t) 0);
        }
      }
    }
    template<typename E, typename C, C M, typename I, enable_scalar_operator_t<I> = nullptr>
    void abs_sub(big_integer<E, C, M>& lhs, I const& rhs) {
      typedef big_integer<E, C, M> integer_t;
      typedef typename integer_t::element_type elem_t;

      int const cmp = abs_compare(lhs, rhs);
      if (cmp == 0) {
        lhs.sign = 0;
        lhs.data.resize(1, (elem_t) 0);
        return;
      }

      typename std::make_unsigned<I>::type urhs = _abs(rhs);
      if (cmp > 0) {
        std::vector<E>& ldata = lhs.data;
        std::size_t const lsize = ldata.size();
        for (std::size_t i = 0; i < lsize; i++) {
          if (urhs == 0) return;
          elem_t const part = (elem_t) (urhs % integer_t::modulo);
          urhs /= integer_t::modulo;
          if (ldata[i] >= part)
            ldata[i] -= part;
          else {
            urhs++;
            ldata[i] += (elem_t) (integer_t::modulo - part);
          }
        }

        // Note: cmp > 0 なので非ゼロの要素が必ず見つかるはず。
        std::size_t i = lsize;
        while (ldata[i - 1] == 0) i--;
        if (i < lsize) ldata.resize(i);
      } else {
        lhs.sign = -lhs.sign;
        for (std::size_t i = 0; i < lhs.data.size(); i++) {
          elem_t const part = (elem_t) (urhs % integer_t::modulo);
          urhs /= integer_t::modulo;
          if (part >= lhs.data[i])
            lhs.data[i] = part - lhs.data[i];
          else {
            urhs--;
            lhs.data[i] = part + (elem_t) (integer_t::modulo - lhs.data[i]);
          }
        }
        while (urhs) {
          lhs.data.emplace_back((elem_t) (urhs % integer_t::modulo));
          urhs /= integer_t::modulo;
        }
      }
    }

    template<char Op, typename E, typename C, C M, typename T>
    big_integer<E, C, M>& impl_add_eq(big_integer<E, C, M>& lhs, T const& rhs) {
      int const rsign = _sign(rhs);
      if (rsign == 0) return lhs;
      if (lhs.sign == 0) {
        lhs = rhs;
        if (Op == '-')
          lhs.sign = -lhs.sign;
        return lhs;
      }

      if (lhs.sign == (Op == '-'? -rsign: rsign))
        abs_add(lhs, rhs);
      else
        abs_sub(lhs, rhs);
      return lhs;
    }

    template<typename E, typename C, C M>
    big_integer<E, C, M>& operator+=(big_integer<E, C, M>& lhs, big_integer<E, C, M> const& rhs) {return impl_add_eq<'+'>(lhs, rhs);}
    template<typename E, typename C, C M, typename I, enable_scalar_operator_t<I> = nullptr>
    big_integer<E, C, M>& operator+=(big_integer<E, C, M>& lhs, I const& rhs) {return impl_add_eq<'+'>(lhs, rhs);}
    template<typename E, typename C, C M>
    big_integer<E, C, M>& operator-=(big_integer<E, C, M>& lhs, big_integer<E, C, M> const& rhs) {return impl_add_eq<'-'>(lhs, rhs);}
    template<typename E, typename C, C M, typename I, enable_scalar_operator_t<I> = nullptr>
    big_integer<E, C, M>& operator-=(big_integer<E, C, M>& lhs, I const& rhs) {return impl_add_eq<'-'>(lhs, rhs);}


    template<typename E, typename C, C M, typename T, enable_generic_operator_t<big_integer<E, C, M>, T> = nullptr>
    big_integer<E, C, M> operator+(big_integer<E, C, M> const& lhs, T const& rhs) {
      big_integer<E, C, M> ret;
      if (lhs.sign != 0) ret = lhs;
      ret += rhs;
      return ret;
    }
    template<typename E, typename C, C M, typename I, enable_scalar_operator_t<I> = nullptr>
    big_integer<E, C, M> operator+(I const& lhs, big_integer<E, C, M> const& rhs) {
      big_integer<E, C, M> ret;
      if (rhs.sign != 0) ret = rhs;
      ret += lhs;
      return ret;
    }
    template<typename E, typename C, C M, typename T, enable_generic_operator_t<big_integer<E, C, M>, T> = nullptr>
    big_integer<E, C, M> operator-(big_integer<E, C, M> const& lhs, T const& rhs) {
      big_integer<E, C, M> ret;
      if (lhs.sign != 0) ret = lhs;
      ret -= rhs;
      return ret;
    }
    template<typename E, typename C, C M, typename I, enable_scalar_operator_t<I> = nullptr>
    big_integer<E, C, M> operator-(I const& lhs, big_integer<E, C, M> const& rhs) {
      big_integer<E, C, M> ret;
      if (rhs.sign != 0) {
        ret = rhs;
        ret.sign = -ret.sign;
      }
      ret -= lhs;
      return ret;
    }
  }

  using big_integer_detail::operator+=;
  using big_integer_detail::operator-=;
  using big_integer_detail::operator+;
  using big_integer_detail::operator-;

  //
  // a * b
  //
  // ToDo: Karatsuba, Toom-Cook, Schonhage-Strassen, Furer + 実測比較
  //
  namespace big_integer_detail {
    template<typename E, typename C, C M>
    void add_digit(big_integer<E, C, M>& num, std::size_t pos, typename big_integer<E, C, M>::calculation_type value) {
      using integer_t = big_integer<E, C, M>;
      using calc_t = typename integer_t::calculation_type;
      using elem_t = typename integer_t::element_type;

      if (num.data.size() <= pos) num.data.resize(pos, 0);
      for (std::size_t const posN = num.data.size(); pos < posN; pos++) {
        if (value == 0) return;
        calc_t const elem = num.data[pos] + value % integer_t::modulo;
        value /= integer_t::modulo;
        num.data[pos] = (elem_t) (elem % integer_t::modulo);
        value += elem / integer_t::modulo;
      }
      while (value) {
        num.data.emplace_back(value % integer_t::modulo);
        value /= integer_t::modulo;
      }
    }

    // 絶対値について num -= value * M^pos を計算する。
    // 前提: M^{pos-1} <= num
    template<typename E, typename C, C M>
    void sub_digit(big_integer<E, C, M>& num, std::size_t pos, typename big_integer<E, C, M>::calculation_type value) {
      using integer_t = big_integer<E, C, M>;
      using calc_t = typename integer_t::calculation_type;
      using elem_t = typename integer_t::element_type;

      // Note: 引く数が 0 の時は pos < num.data.size() は必ずしも成り立たない事に注意する。
      mwg_check(value == 0 || pos < num.data.size(), "pos=%zd |num.data|=%zd", pos, num.data.size());

      for (std::size_t const posN = num.data.size(); pos < posN; pos++) {
        if (value == 0) return;
        calc_t const elem = value % integer_t::modulo;
        value = value / integer_t::modulo;
        if (num.data[pos] >= elem) {
          num.data[pos] -= elem;
        } else {
          value++;
          num.data[pos] += (elem_t) (integer_t::modulo - elem);
        }
      }

      if (value != 0) {
        num.sign = -num.sign;
        if (--value) num.data.emplace_back((elem_t) value);
        for (std::size_t i = num.data.size(); --i > 0; )
          num.data[i] = (elem_t) (integer_t::modulo - 1 - num.data[i]);
        num.data[0] = (elem_t) (integer_t::modulo - num.data[0]);
      } else if (std::size_t sz = num.data.size()) {
        while (sz && num.data[sz - 1] == 0) sz--;
        num.data.resize(sz);
        if (sz == 0) num.sign = 0;
      }
    }

    template<typename E, typename C, C M>
    big_integer<E, C, M> operator*(big_integer<E, C, M> const& lhs, big_integer<E, C, M> const& rhs) {
      typedef big_integer<E, C, M> integer_t;
      using calc_t = typename integer_t::calculation_type;

      big_integer<E, C, M> ret;
      if (lhs.sign == 0 || rhs.sign == 0) return ret;
      ret.sign = lhs.sign * rhs.sign;

      std::size_t reservedSize = lhs.data.size() + rhs.data.size() - 1;
      if (((calc_t) lhs.data.back() + 1) * ((calc_t) rhs.data.back() + 1) > integer_t::modulo)
        reservedSize++;
      ret.data.reserve(reservedSize);

      std::size_t const rN = rhs.data.size();
      std::size_t const lN = lhs.data.size();
      for (std::size_t r = 0; r < rN; r++)
        for (std::size_t l = 0; l < lN; l++)
          add_digit(ret, l + r, (calc_t) lhs.data[l] * (calc_t) rhs.data[r]);

      return ret;
    }

    template<typename Integer, typename Modulo>
    constexpr int integral_ceil_log(Integer value, Modulo const& modulo) {
      int count = 0;
      while (value) {
        value /= modulo;
        count++;
      }
      return count;
    }

    template<typename Integer, typename Modulo>
    constexpr int integral_floor_log(Integer value, Modulo const& modulo) {
      int count = 0;
      while (value >= modulo) {
        value /= modulo;
        count++;
      }
      return count;
    }

    template<typename E, typename C, C M, typename I, nullptr_if_t<std::is_integral<I>::value> = nullptr>
    big_integer<E, C, M> operator*(big_integer<E, C, M> const& lhs, I const& rhs) {
      typedef big_integer<E, C, M> integer_t;
      using calc_t = typename integer_t::calculation_type;
      using urhs_t = typename std::make_unsigned<I>::type;

      big_integer<E, C, M> ret;
      if (lhs.sign == 0 || rhs == 0) return ret;
      ret.sign = lhs.sign;

      urhs_t urhs = rhs;
      if (std::is_signed<I>::value && rhs < 0) {
        ret.sign = -ret.sign;
        urhs = -urhs;
      }

      constexpr int offset = integral_ceil_log(std::numeric_limits<urhs_t>::max(), integer_t::modulo);
      ret.data.reserve(lhs.data.size() + offset);

      std::size_t const lN = lhs.data.size();
      for (std::size_t r = 0; urhs; r++) {
        calc_t const elem = urhs % integer_t::modulo;
        urhs /= integer_t::modulo;
        for (std::size_t l = 0; l < lN; l++)
          add_digit(ret, l + r, (calc_t) lhs.data[l] * elem);
      }

      return ret;
    }
    template<typename E, typename C, C M, typename I, nullptr_if_t<std::is_integral<I>::value> = nullptr>
    big_integer<E, C, M> operator*(I const& lhs, big_integer<E, C, M> const& rhs) {return rhs * lhs;}
    template<typename E, typename C, C M>
    big_integer<E, C, M> operator*=(big_integer<E, C, M>& lhs, big_integer<E, C, M> const& rhs) {return lhs = lhs * rhs;}
    template<typename E, typename C, C M, typename I, nullptr_if_t<std::is_integral<I>::value> = nullptr>
    big_integer<E, C, M> operator*=(big_integer<E, C, M>& lhs, I const& rhs) {return lhs = lhs * rhs;}
  }

  using big_integer_detail::operator*;
  using big_integer_detail::operator*=;

  //
  // pow(a, u)
  //
  template<typename E, typename C, C M>
  big_integer<E, C, M> pow(big_integer<E, C, M> const& lhs, unsigned exponent) {
    big_integer<E, C, M> result {1};
    if (!exponent) return result;

    big_integer<E, C, M> pow2 = lhs;
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

  namespace big_integer_detail {

    template<typename E, typename C, C M>
    E divide(big_integer<E, C, M> const& lhs, E const& rhs, big_integer<E, C, M>* pquo) {
      typedef big_integer<E, C, M> integer_t;
      using calc_t = typename integer_t::calculation_type;
      using elem_t = typename integer_t::element_type;

      if (rhs == 0) {
        throw mwg::except("big_integer: division by zero", mwg::ecode::EArgRange);
      } else if (rhs == 1) {
        if (pquo && pquo != &lhs) *pquo = lhs;
        return 0;
      }

      calc_t carry = 0;
      if (pquo) {
        pquo->sign = lhs.sign;
        std::vector<elem_t>& data = pquo->data;
        data.resize(lhs.data.size());
        for (std::size_t i = lhs.data.size(); i--; ) {
          calc_t const elem = carry * integer_t::modulo + lhs.data[i];
          data[i] = elem / rhs;
          carry = elem % rhs;
        }
        if (data.back() == 0) data.pop_back();
      } else {
        for (std::size_t i = lhs.data.size(); i--; )
          carry = (carry * integer_t::modulo + lhs.data[i]) % rhs;
      }
      return carry;
    }
    template<typename E, typename C, C M>
    E operator%(big_integer<E, C, M> const& lhs, E const& rhs) {
      return lhs.sign * divide(lhs, rhs, (big_integer<E, C, M>*) nullptr);
    }
    template<typename E, typename C, C M>
    big_integer<E, C, M> operator/(big_integer<E, C, M> const& lhs, E const& rhs) {
      big_integer<E, C, M> ret;
      divide(lhs, rhs, &ret);
      return ret;
    }
    template<typename E, typename C, C M>
    big_integer<E, C, M>& operator%=(big_integer<E, C, M>& lhs, E const& rhs) {
      E const rem = divide(lhs, rhs, (big_integer<E, C, M>*) nullptr);
      lhs = lhs.sign * rem;
      return lhs;
    }
    template<typename E, typename C, C M>
    big_integer<E, C, M>& operator/=(big_integer<E, C, M>& lhs, E const& rhs) {
      divide(lhs, rhs, &lhs);
      return lhs;
    }

    // num - div * M^exp * factor を計算する
    // 前提: 結果は正である。
    template<typename E, typename C, C M>
    void abs_submul(
      big_integer<E, C, M>& num, big_integer<E, C, M> const& div,
      std::size_t exp, typename big_integer<E, C, M>::calculation_type factor
    ) {
      for (std::size_t i = div.data.size(); i--; )
        sub_digit(num, i + exp, factor * div.data[i]);
    }

    // Knuth の方法
    template<typename E, typename C, C M>
    bool divide(
      big_integer<E, C, M> const& lhs, big_integer<E, C, M> const& rhs,
      big_integer<E, C, M>* pquo, big_integer<E, C, M>* prem
    ) {
      typedef big_integer<E, C, M> integer_t;
      using calc_t = typename integer_t::calculation_type;
      using elem_t = typename integer_t::element_type;

      if (rhs.sign == 0) {
        throw mwg::except("big_integer: division by zero", mwg::ecode::EArgRange);
      } else if (lhs.sign == 0 || abs_compare(lhs, rhs) < 0) {
        bool const ret = lhs.sign == 0;
        if (prem && prem != &lhs) *prem = lhs;
        if (pquo) *pquo = 0;
        return ret;
      } else if (rhs.data.size() == 1) {
        elem_t const rem = divide(lhs, rhs.data[0], pquo);
        if (prem) *prem = rem;
        return rem == 0;
      }

      std::vector<E> const& rhsdata = rhs.data;
      std::size_t const rhssize = rhs.data.size();
      mwg_assert(rhssize >= 2);
      calc_t const rhs2 = rhsdata[rhssize - 1] * integer_t::modulo + rhsdata[rhssize - 2];
      calc_t const factor =
        rhs2 == std::numeric_limits<calc_t>::max()? 1:
        (integer_t::modulo * integer_t::modulo - 1) / (rhs2 + 1);
      mwg_assert(1 <= factor && factor < integer_t::modulo);

      // Note: prem == &lhs/&rhs の場合、
      //   以下の ※ の操作で lhs/rhs が破壊される。
      //   従って ※ より後では lhs/rhs には触れない事に注意する。
      integer_t _rem, div;
      integer_t& rem = prem? *prem: _rem;
      if (factor == 1) {
        div = rhs;
        rem = lhs;
      } else {
        div = rhs * factor;
        rem = lhs * factor; // ※
      }

      std::vector<E> const& ddata = div.data;
      std::vector<E>      & ndata = rem.data;
      calc_t const arhs = (calc_t) ddata.back() + 1;
      std::size_t rpos = ndata.size() - 1;
      std::size_t qpos = ndata.size() - ddata.size();
      mwg_assert(arhs != 0);

      // 最上位の桁
      if (pquo) pquo->sign = div.sign * rem.sign;
      if (ndata[rpos] >= arhs) {
        if (pquo) {
          pquo->data.resize(qpos + 1, 0);
          pquo->data[qpos] = 1;
        }
        abs_submul(rem, div, qpos, 1);
      } else {
        if (pquo) pquo->data.resize(qpos, 0);
      }

      // 残りの桁
      while (rpos--, qpos--) {
        calc_t const aquot = (ndata[rpos + 1] * integer_t::modulo + ndata[rpos]) / arhs;
        if (pquo) add_digit(*pquo, qpos, aquot);
        abs_submul(rem, div, qpos, aquot);

        // 誤差修正
        if (abs_compare(rem, div, qpos) >= 0) {
          if (pquo) add_digit(*pquo, qpos, 1);
          abs_submul(rem, div, qpos, 1);
        }
      }

      // Note: roundup(div) を使って引き算していくので、
      //   div <= rem < roundup(div) なる rem が残ることがある。
      //   但し roundup() は有効数字第二位で切り上げ。
      //   残っていたら引き算する。
      if (abs_compare(rem, div) >= 0) {
        if (pquo) add_digit(*pquo, 0, 1);
        rem -= div;
      }

      if (prem && factor > 1) rem /= (elem_t) factor;
      return rem == 0;
    }

    template<typename E, typename C, C M>
    big_integer<E, C, M> operator%(big_integer<E, C, M> const& lhs, big_integer<E, C, M> const& rhs) {
      big_integer<E, C, M> ret;
      divide(lhs, rhs, (big_integer<E, C, M>*) nullptr, &ret);
      return ret;
    }
    template<typename E, typename C, C M>
    big_integer<E, C, M> operator/(big_integer<E, C, M> const& lhs, big_integer<E, C, M> const& rhs) {
      big_integer<E, C, M> ret;
      divide(lhs, rhs, &ret, (big_integer<E, C, M>*) nullptr);
      return ret;
    }
    template<typename E, typename C, C M>
    big_integer<E, C, M>& operator%=(big_integer<E, C, M>& lhs, big_integer<E, C, M> const& rhs) {
      divide(lhs, rhs, (big_integer<E, C, M>*) nullptr, &lhs);
      return lhs;
    }
    template<typename E, typename C, C M>
    big_integer<E, C, M>& operator/=(big_integer<E, C, M>& lhs, big_integer<E, C, M> const& rhs) {
      divide(lhs, rhs, &lhs, (big_integer<E, C, M>*) nullptr);
      return lhs;
    }
  }
  using big_integer_detail::operator%;
  using big_integer_detail::operator/;
  using big_integer_detail::operator%=;
  using big_integer_detail::operator/=;

  namespace big_integer_detail {

    template<typename E>
    constexpr E integral_pow(E base, unsigned index) {
      E result = 1;
      while (index--) result *= base;
      return result;
    }

    template<typename E, typename C, C M>
    std::ostream& operator<<(std::ostream& ostr, big_integer<E, C, M> const& value) {
      if (value.sign == 0) return ostr << '0';
      if (value.sign == -1) ostr << '-';

      using integer_t = big_integer<E, C, M>;
      using elem_t = typename integer_t::element_type;
      using calc_t = typename integer_t::calculation_type;

      constexpr int wgroup = integral_floor_log(integer_t::modulo, 10u);
      constexpr elem_t mod10 = integral_pow((elem_t) 10, wgroup);
      if (mod10 == integer_t::modulo) {
        // 元から 10^n の基数の場合には変換は不要
        std::vector<char> digits;
        digits.resize(wgroup);

        // 最初の桁
        std::size_t i = value.data.size() - 1;
        elem_t elem = value.data[i];
        int k = wgroup;
        while (elem && k--) {
          digits[k] = '0' + elem % 10;
          elem /= 10;
        }
        ostr.write(&digits[k], wgroup - k);

        // 残りの桁
        while (i--) {
          elem_t elem = value.data[i];
          for (int k = wgroup; k--; ) {
            digits[k] = '0' + elem % 10;
            elem /= 10;
          }
          ostr.write(&digits[0], digits.size());
        }
      } else {
        integer_t tmp = value;
        std::vector<char> digits;
        digits.reserve(tmp.data.size() * integral_ceil_log(integer_t::modulo, 10u));
        for (;;) {
          calc_t carry = big_integer_detail::divide(tmp, mod10, &tmp);
          if (tmp.data.size() != 0) {
            for (int i = 0; i < wgroup; i++) {
              digits.push_back((char)('0' + carry % 10));
              carry /= 10;
            }
          } else {
            while (carry) {
              digits.push_back((char)('0' + carry % 10));
              carry /= 10;
            }
            break;
          }
        }

        std::reverse(digits.begin(), digits.end());
        ostr.write(&digits[0], digits.size());
      }

      return ostr;
    }
  }

  using big_integer_detail::operator<<;

  typedef big_integer<> bigint;

}
#endif
