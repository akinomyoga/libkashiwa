// -*- C++ -*-
#ifndef KASHIWA_BIGINTEGER_H
#define KASHIWA_BIGINTEGER_H
#include <cstddef>
#include <cstdint>
#include <limits>
#include <type_traits>
#include <vector>
#include <ostream>
#include <algorithm>
#include <mwg/except.h>
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

    template<typename I, typename std::enable_if<std::is_integral<I>::value, std::nullptr_t>::type = nullptr>
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
    template<typename I, typename std::enable_if<std::is_integral<I>::value, std::nullptr_t>::type = nullptr>
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
    using enable_generic_operator_t = typename std::enable_if<
      (std::is_integral<T>::value || std::is_base_of<MpInteger, T>::value), std::nullptr_t>::type;
    template<typename T>
    using enable_scalar_operator_t = typename std::enable_if<
      std::is_integral<T>::value, std::nullptr_t>::type;

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
  // compare(a, b)
  //
  namespace big_integer_detail {
    template<typename I, typename std::enable_if<std::is_integral<I>::value && std::is_signed<I>::value, std::nullptr_t>::type = nullptr>
    typename std::make_unsigned<I>::type _abs(I const& value) {
      typename std::make_unsigned<I>::type uvalue = value;
      if (value < 0) uvalue = -uvalue;
      return uvalue;
    }
    template<typename U, typename std::enable_if<std::is_integral<U>::value && std::is_unsigned<U>::value, std::nullptr_t>::type = nullptr>
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
  }

  using big_integer_detail::compare;

  //
  // a == b, a < b, etc.
  //
  namespace big_integer_detail {
    template<typename E, typename C, C M>
    int operator==(big_integer<E, C, M> const& lhs, big_integer<E, C, M> const& rhs) {return compare(lhs, rhs) == 0;}
    template<typename E, typename C, C M, typename I, enable_scalar_operator_t<I> = nullptr>
    int operator==(big_integer<E, C, M> const& lhs, I const& rhs) {return compare(lhs, rhs) == 0;}
    template<typename E, typename C, C M, typename I, enable_scalar_operator_t<I> = nullptr>
    int operator==(I const& lhs, big_integer<E, C, M> const& rhs) {return compare(rhs, lhs) == 0;}
    template<typename E, typename C, C M>
    int operator!=(big_integer<E, C, M> const& lhs, big_integer<E, C, M> const& rhs) {return compare(lhs, rhs) != 0;}
    template<typename E, typename C, C M, typename I, enable_scalar_operator_t<I> = nullptr>
    int operator!=(big_integer<E, C, M> const& lhs, I const& rhs) {return compare(lhs, rhs) != 0;}
    template<typename E, typename C, C M, typename I, enable_scalar_operator_t<I> = nullptr>
    int operator!=(I const& lhs, big_integer<E, C, M> const& rhs) {return compare(rhs, lhs) != 0;}
    template<typename E, typename C, C M>
    int operator<(big_integer<E, C, M> const& lhs, big_integer<E, C, M> const& rhs) {return compare(lhs, rhs) < 0;}
    template<typename E, typename C, C M, typename I, enable_scalar_operator_t<I> = nullptr>
    int operator<(big_integer<E, C, M> const& lhs, I const& rhs) {return compare(lhs, rhs) < 0;}
    template<typename E, typename C, C M, typename I, enable_scalar_operator_t<I> = nullptr>
    int operator<(I const& lhs, big_integer<E, C, M> const& rhs) {return compare(rhs, lhs) < 0;}
    template<typename E, typename C, C M>
    int operator>(big_integer<E, C, M> const& lhs, big_integer<E, C, M> const& rhs) {return compare(lhs, rhs) > 0;}
    template<typename E, typename C, C M, typename I, enable_scalar_operator_t<I> = nullptr>
    int operator>(big_integer<E, C, M> const& lhs, I const& rhs) {return compare(lhs, rhs) > 0;}
    template<typename E, typename C, C M, typename I, enable_scalar_operator_t<I> = nullptr>
    int operator>(I const& lhs, big_integer<E, C, M> const& rhs) {return compare(rhs, lhs) > 0;}
    template<typename E, typename C, C M>
    int operator<=(big_integer<E, C, M> const& lhs, big_integer<E, C, M> const& rhs) {return compare(lhs, rhs) <= 0;}
    template<typename E, typename C, C M, typename I, enable_scalar_operator_t<I> = nullptr>
    int operator<=(big_integer<E, C, M> const& lhs, I const& rhs) {return compare(lhs, rhs) <= 0;}
    template<typename E, typename C, C M, typename I, enable_scalar_operator_t<I> = nullptr>
    int operator<=(I const& lhs, big_integer<E, C, M> const& rhs) {return compare(rhs, lhs) <= 0;}
    template<typename E, typename C, C M>
    int operator>=(big_integer<E, C, M> const& lhs, big_integer<E, C, M> const& rhs) {return compare(lhs, rhs) >= 0;}
    template<typename E, typename C, C M, typename I, enable_scalar_operator_t<I> = nullptr>
    int operator>=(big_integer<E, C, M> const& lhs, I const& rhs) {return compare(lhs, rhs) >= 0;}
    template<typename E, typename C, C M, typename I, enable_scalar_operator_t<I> = nullptr>
    int operator>=(I const& lhs, big_integer<E, C, M> const& rhs) {return compare(rhs, lhs) >= 0;}
  }
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
        for (std::size_t i = 0; i < lhs.data.size(); i++) {
          if (urhs == 0) return;
          elem_t const part = (elem_t) (urhs % integer_t::modulo);
          urhs /= integer_t::modulo;
          if (lhs.data[i] >= part)
            lhs.data[i] -= part;
          else {
            urhs++;
            lhs.data[i] += (elem_t) (integer_t::modulo - part);
          }
        }
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

      if (num.data.size() <= pos) num.data.resize(pos, 0);
      for (std::size_t const posN = num.data.size(); pos < posN; pos++) {
        if (value == 0) return;
        calc_t const elem = num.data[pos] + value % integer_t::modulo;
        value /= integer_t::modulo;
        num.data[pos] = elem % integer_t::modulo;
        value += elem / integer_t::modulo;
      }
      while (value) {
        num.data.emplace_back(value % integer_t::modulo);
        value /= integer_t::modulo;
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
    constexpr int integral_log(Integer value, Modulo const& modulo) {
      int count = 0;
      while (value) {
        value /= modulo;
        count++;
      }
      return count;
    }

    template<typename E, typename C, C M, typename I, typename std::enable_if<std::is_integral<I>::value, std::nullptr_t>::type = nullptr>
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

      constexpr int offset = integral_log(std::numeric_limits<urhs_t>::max(), integer_t::modulo);
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
    template<typename E, typename C, C M, typename I, typename std::enable_if<std::is_integral<I>::value, std::nullptr_t>::type = nullptr>
    big_integer<E, C, M> operator*(I const& lhs, big_integer<E, C, M> const& rhs) {return rhs * lhs;}
    template<typename E, typename C, C M>
    big_integer<E, C, M> operator*=(big_integer<E, C, M>& lhs, big_integer<E, C, M> const& rhs) {return lhs = lhs * rhs;}
    template<typename E, typename C, C M, typename I, typename std::enable_if<std::is_integral<I>::value, std::nullptr_t>::type = nullptr>
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
    E divide(big_integer<E, C, M> const& lhs, E const& rhs, big_integer<E, C, M>& quot) {
      typedef big_integer<E, C, M> integer_t;
      using calc_t = typename integer_t::calculation_type;
      using elem_t = typename integer_t::element_type;

      if (rhs == 0) {
        throw mwg::except("big_integer: division by zero", mwg::ecode::EArgRange);
      } else if (rhs == 1) {
        if (&lhs != &quot) quot.data = lhs.data;
        return 0;
      }

      std::vector<elem_t>& data = quot.data;
      if (&quot != &lhs) data.resize(lhs.data.size());
      calc_t carry = 0;
      for (std::size_t i = lhs.data.size(); i--; ) {
        calc_t const elem = carry * integer_t::modulo + lhs.data[i];
        data[i] = elem / 10;
        carry = elem % 10;
      }
      if (data.back() == 0) data.pop_back();
      return carry;
    }
    template<typename E, typename C, C M>
    E divide(big_integer<E, C, M> const& lhs, E const& rhs) {
      typedef big_integer<E, C, M> integer_t;
      using calc_t = typename integer_t::calculation_type;

      if (rhs == 0)
        throw mwg::except("big_integer: division by zero", mwg::ecode::EArgRange);
      else if (rhs == 1)
        return 0;

      calc_t carry = 0;
      for (std::size_t i = lhs.data.size(); i--; )
        carry = (carry * integer_t::modulo + lhs.data[i]) % 10;
      return carry;
    }
    template<typename E, typename C, C M>
    E operator%(big_integer<E, C, M> const& lhs, E const& rhs) {
      return divide(lhs, rhs);
    }
    template<typename E, typename C, C M>
    big_integer<E, C, M> operator/(big_integer<E, C, M> const& lhs, E const& rhs) {
      big_integer<E, C, M> ret;
      divide(lhs, rhs, ret);
      return ret;
    }
    template<typename E, typename C, C M>
    big_integer<E, C, M>& operator%=(big_integer<E, C, M>& lhs, E const& rhs) {
      lhs = divide(lhs, rhs);
      return lhs;
    }
    template<typename E, typename C, C M>
    big_integer<E, C, M>& operator/=(big_integer<E, C, M>& lhs, E const& rhs) {
      divide(lhs, rhs, lhs);
      return lhs;
    }
  }

  using big_integer_detail::operator%;
  using big_integer_detail::operator/;
  using big_integer_detail::operator%=;
  using big_integer_detail::operator/=;

  namespace big_integer_detail {
    template<typename E, typename C, C M>
    std::ostream& operator<<(std::ostream& ostr, big_integer<E, C, M> const& value) {
      if (value.sign == 0) return ostr << '0';
      if (value.sign == -1) ostr << '-';

      using integer_t = big_integer<E, C, M>;
      using elem_t = typename integer_t::element_type;
      using calc_t = typename integer_t::calculation_type;

      std::vector<elem_t> tmp = value.data;
      std::vector<char> digits;
      digits.reserve(tmp.size() * integral_log(integer_t::modulo, 10));
      while (tmp.size()) {
        calc_t carry = 0;
        for (std::size_t i = tmp.size(); i--; ) {
          calc_t const elem = carry * integer_t::modulo + tmp[i];
          tmp[i] = elem / 10;
          carry = elem % 10;
        }
        if (tmp.back() == 0) tmp.pop_back();
        digits.push_back((char)('0' + carry));
      }

      std::reverse(digits.begin(), digits.end());
      ostr.write(&digits[0], digits.size());
      return ostr;
    }
  }

  using big_integer_detail::operator<<;

  typedef big_integer<> bigint;

}
#endif
