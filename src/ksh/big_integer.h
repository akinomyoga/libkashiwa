// -*- C++ -*-
#ifndef KASHIWA_BIGINTEGER_H
#define KASHIWA_BIGINTEGER_H
#include <cstddef>
#include <cstdint>
#include <limits>
#include <type_traits>
#include <vector>
#include <mwg/except.h>
namespace kashiwa {

  template<
    typename StoreInt = std::uint32_t,
    typename CalcInt = std::uint64_t,
    CalcInt modulo = (CalcInt) std::numeric_limits<StoreInt>::max() + 1 >
  struct big_integer;

  template<typename StoreInt, typename CalcInt, CalcInt Modulo>
  struct big_integer {
    typedef StoreInt element_type;
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
    big_integer(int value) {this->_set_integral_value(value);}

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
    void _set_integral_value(I const& value) {
      data.clear();

      if (value == 0) {
        sign = 0;
        return;
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
    }
  };

  namespace big_integer_detail {
    template<typename MpInteger, typename T>
    using enable_generic_operator_t = typename std::enable_if<
      (std::is_integral<T>::value || std::is_base_of<MpInteger, T>::value), std::nullptr_t>::type;
    template<typename T>
    using enable_scalar_operator_t = typename std::enable_if<
      std::is_integral<T>::value, std::nullptr_t>::type;

    template<typename S, typename C, C M>
    int _sign(big_integer<S, C, M> const&  arg) {return arg.sign;}
    template<typename I, enable_scalar_operator_t<I> = nullptr>
    int _sign(I const&  arg) {
      if (arg > 0)
        return 1;
      else if (std::is_signed<I>::value && arg < 0)
        return -1;
      else
        return 0;
    }

    template<typename S, typename C, C M>
    void _set(big_integer<S, C, M>& lhs, big_integer<S, C, M> const& rhs) {
      lhs = rhs;
    }
    template<typename S, typename C, C M, typename I, enable_scalar_operator_t<I> = nullptr>
    void _set(big_integer<S, C, M>& lhs, I const& rhs) {
      lhs._set_integral_value(rhs);
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

    template<typename S, typename C, C M>
    int abs_compare(big_integer<S, C, M> const& lhs, big_integer<S, C, M> const& rhs) {
      if (lhs.data.size() != rhs.data.size())
        return lhs.data.size() < rhs.data.size()? -1: 1;

      for (std::size_t i = lhs.data.size(); i--; )
        if (lhs.data[i] != rhs.data[i])
          return lhs.data[i] < rhs.data[i]? -1: 1;

      return 0;
    }
    template<typename S, typename C, C M, typename I, enable_scalar_operator_t<I> = nullptr>
    int abs_compare(big_integer<S, C, M> const& lhs, I const& rhs) {
      using integer_t = big_integer<S, C, M>;
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

    template<typename S, typename C, C M, typename T>
    int impl_compare(big_integer<S, C, M> const& lhs, T const& rhs) {
      int const rsign = _sign(rhs);
      if (rsign > 0) {
        if (lhs.sign <= 0) return -1;
      } else if (rsign < 0) {
        if (lhs.sign >= 0) return 1;
      } else
        return lhs.sign;
      return rsign * abs_compare(lhs, rhs);
    }
    template<typename S, typename C, C M>
    int compare(big_integer<S, C, M> const& lhs, big_integer<S, C, M> const& rhs) {return impl_compare(lhs, rhs);}
    template<typename S, typename C, C M, typename I, enable_scalar_operator_t<I> = nullptr>
    int compare(big_integer<S, C, M> const& lhs, I const& rhs) {return impl_compare(lhs, rhs);}
    template<typename S, typename C, C M, typename I, enable_scalar_operator_t<I> = nullptr>
    int compare(I const& lhs, big_integer<S, C, M> const& rhs) {return -impl_compare(rhs, lhs);}
  }

  using big_integer_detail::compare;

  //
  // a == b, a < b, etc.
  //
  namespace big_integer_detail {
    template<typename S, typename C, C M>
    int operator==(big_integer<S, C, M> const& lhs, big_integer<S, C, M> const& rhs) {return compare(lhs, rhs) == 0;}
    template<typename S, typename C, C M, typename I, enable_scalar_operator_t<I> = nullptr>
    int operator==(big_integer<S, C, M> const& lhs, I const& rhs) {return compare(lhs, rhs) == 0;}
    template<typename S, typename C, C M, typename I, enable_scalar_operator_t<I> = nullptr>
    int operator==(I const& lhs, big_integer<S, C, M> const& rhs) {return compare(rhs, lhs) == 0;}
    template<typename S, typename C, C M>
    int operator!=(big_integer<S, C, M> const& lhs, big_integer<S, C, M> const& rhs) {return compare(lhs, rhs) != 0;}
    template<typename S, typename C, C M, typename I, enable_scalar_operator_t<I> = nullptr>
    int operator!=(big_integer<S, C, M> const& lhs, I const& rhs) {return compare(lhs, rhs) != 0;}
    template<typename S, typename C, C M, typename I, enable_scalar_operator_t<I> = nullptr>
    int operator!=(I const& lhs, big_integer<S, C, M> const& rhs) {return compare(rhs, lhs) != 0;}
    template<typename S, typename C, C M>
    int operator<(big_integer<S, C, M> const& lhs, big_integer<S, C, M> const& rhs) {return compare(lhs, rhs) < 0;}
    template<typename S, typename C, C M, typename I, enable_scalar_operator_t<I> = nullptr>
    int operator<(big_integer<S, C, M> const& lhs, I const& rhs) {return compare(lhs, rhs) < 0;}
    template<typename S, typename C, C M, typename I, enable_scalar_operator_t<I> = nullptr>
    int operator<(I const& lhs, big_integer<S, C, M> const& rhs) {return compare(rhs, lhs) < 0;}
    template<typename S, typename C, C M>
    int operator>(big_integer<S, C, M> const& lhs, big_integer<S, C, M> const& rhs) {return compare(lhs, rhs) > 0;}
    template<typename S, typename C, C M, typename I, enable_scalar_operator_t<I> = nullptr>
    int operator>(big_integer<S, C, M> const& lhs, I const& rhs) {return compare(lhs, rhs) > 0;}
    template<typename S, typename C, C M, typename I, enable_scalar_operator_t<I> = nullptr>
    int operator>(I const& lhs, big_integer<S, C, M> const& rhs) {return compare(rhs, lhs) > 0;}
    template<typename S, typename C, C M>
    int operator<=(big_integer<S, C, M> const& lhs, big_integer<S, C, M> const& rhs) {return compare(lhs, rhs) <= 0;}
    template<typename S, typename C, C M, typename I, enable_scalar_operator_t<I> = nullptr>
    int operator<=(big_integer<S, C, M> const& lhs, I const& rhs) {return compare(lhs, rhs) <= 0;}
    template<typename S, typename C, C M, typename I, enable_scalar_operator_t<I> = nullptr>
    int operator<=(I const& lhs, big_integer<S, C, M> const& rhs) {return compare(rhs, lhs) <= 0;}
    template<typename S, typename C, C M>
    int operator>=(big_integer<S, C, M> const& lhs, big_integer<S, C, M> const& rhs) {return compare(lhs, rhs) >= 0;}
    template<typename S, typename C, C M, typename I, enable_scalar_operator_t<I> = nullptr>
    int operator>=(big_integer<S, C, M> const& lhs, I const& rhs) {return compare(lhs, rhs) >= 0;}
    template<typename S, typename C, C M, typename I, enable_scalar_operator_t<I> = nullptr>
    int operator>=(I const& lhs, big_integer<S, C, M> const& rhs) {return compare(rhs, lhs) >= 0;}
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

    template<typename S, typename C, C M>
    void abs_add(big_integer<S, C, M>& lhs, big_integer<S, C, M> const& rhs) {
      typedef big_integer<S, C, M> integer_t;
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
    template<typename S, typename C, C M, typename I, enable_scalar_operator_t<I> = nullptr>
    void abs_add(big_integer<S, C, M>& lhs, I const& rhs) {
      using integer_t = big_integer<S, C, M>;
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

    template<typename S, typename C, C M>
    void abs_sub(big_integer<S, C, M>& lhs, big_integer<S, C, M> const& rhs) {
      typedef big_integer<S, C, M> integer_t;
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
    template<typename S, typename C, C M, typename I, enable_scalar_operator_t<I> = nullptr>
    void abs_sub(big_integer<S, C, M>& lhs, I const& rhs) {
      typedef big_integer<S, C, M> integer_t;
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

    template<char Op, typename S, typename C, C M, typename T>
    big_integer<S, C, M>& impl_add_eq(big_integer<S, C, M>& lhs, T const& rhs) {
      int const rsign = _sign(rhs);
      if (rsign == 0) return lhs;
      if (lhs.sign == 0) {
        _set(lhs, rhs);
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

    template<typename S, typename C, C M>
    big_integer<S, C, M>& operator+=(big_integer<S, C, M>& lhs, big_integer<S, C, M> const& rhs) {return impl_add_eq<'+'>(lhs, rhs);}
    template<typename S, typename C, C M, typename I, enable_scalar_operator_t<I> = nullptr>
    big_integer<S, C, M>& operator+=(big_integer<S, C, M>& lhs, I const& rhs) {return impl_add_eq<'+'>(lhs, rhs);}
    template<typename S, typename C, C M>
    big_integer<S, C, M>& operator-=(big_integer<S, C, M>& lhs, big_integer<S, C, M> const& rhs) {return impl_add_eq<'-'>(lhs, rhs);}
    template<typename S, typename C, C M, typename I, enable_scalar_operator_t<I> = nullptr>
    big_integer<S, C, M>& operator-=(big_integer<S, C, M>& lhs, I const& rhs) {return impl_add_eq<'-'>(lhs, rhs);}


    template<typename S, typename C, C M, typename T, enable_generic_operator_t<big_integer<S, C, M>, T> = nullptr>
    big_integer<S, C, M> operator+(big_integer<S, C, M> const& lhs, T const& rhs) {
      big_integer<S, C, M> ret;
      if (lhs.sign != 0) _set(ret, lhs);
      ret += rhs;
      return ret;
    }
    template<typename S, typename C, C M, typename I, enable_scalar_operator_t<I> = nullptr>
    big_integer<S, C, M> operator+(I const& lhs, big_integer<S, C, M> const& rhs) {
      big_integer<S, C, M> ret;
      if (rhs.sign != 0) _set(ret, rhs);
      ret += lhs;
      return ret;
    }
    template<typename S, typename C, C M, typename T, enable_generic_operator_t<big_integer<S, C, M>, T> = nullptr>
    big_integer<S, C, M> operator-(big_integer<S, C, M> const& lhs, T const& rhs) {
      big_integer<S, C, M> ret;
      if (lhs.sign != 0) _set(ret, lhs);
      ret -= rhs;
      return ret;
    }
    template<typename S, typename C, C M, typename I, enable_scalar_operator_t<I> = nullptr>
    big_integer<S, C, M> operator-(I const& lhs, big_integer<S, C, M> const& rhs) {
      big_integer<S, C, M> ret;
      if (rhs.sign != 0) {
        _set(ret, rhs);
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

}
#endif
