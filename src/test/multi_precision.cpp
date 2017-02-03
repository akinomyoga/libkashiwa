#include <cstddef>
#include <cstdio>
#include <cstdint>
#include <limits>
#include <vector>
#include <algorithm>
#include <type_traits>
#include <mwg/except.h>

namespace multi_precision_i1 {

  template<
    typename StoreInt = std::uint32_t,
    typename CalcInt = std::uint64_t,
    CalcInt modulo = (CalcInt) std::numeric_limits<StoreInt>::max() + 1 >
  struct mp_integer;

  template<typename StoreInt, typename CalcInt, CalcInt Modulo>
  struct mp_integer {
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

    mp_integer(): sign(0) {}

    mp_integer(int value) {this->_set_integral_value(value);}

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

  namespace mp_integer_detail {
    template<typename MpInteger, typename T>
    using enable_generic_operator_t = typename std::enable_if<
      (std::is_integral<T>::value || std::is_base_of<MpInteger, T>::value), std::nullptr_t>::type;
    template<typename T>
    using enable_scalar_operator_t = typename std::enable_if<
      std::is_integral<T>::value, std::nullptr_t>::type;

    template<typename S, typename C, C M>
    int _sign(mp_integer<S, C, M> const&  arg) {return arg.sign;}
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
    void _set(mp_integer<S, C, M>& lhs, mp_integer<S, C, M> const& rhs) {
      lhs = rhs;
    }
    template<typename S, typename C, C M, typename I, enable_scalar_operator_t<I> = nullptr>
    void _set(mp_integer<S, C, M>& lhs, I const& rhs) {
      lhs._set_integral_value(rhs);
    }
  }

  //
  // compare(a, b)
  //
  namespace mp_integer_detail {
    template<typename I, typename std::enable_if<std::is_integral<I>::value && std::is_signed<I>::value, std::nullptr_t>::type = nullptr>
    typename std::make_unsigned<I>::type _abs(I const& value) {
      typename std::make_unsigned<I>::type uvalue = value;
      if (value < 0) uvalue = -uvalue;
      return uvalue;
    }
    template<typename U, typename std::enable_if<std::is_integral<U>::value && std::is_unsigned<U>::value, std::nullptr_t>::type = nullptr>
    U const& _abs(U const& value) {return value;}

    template<typename S, typename C, C M>
    int abs_compare(mp_integer<S, C, M> const& lhs, mp_integer<S, C, M> const& rhs) {
      if (lhs.data.size() != rhs.data.size())
        return lhs.data.size() < rhs.data.size()? -1: 1;

      for (std::size_t i = lhs.data.size(); i--; )
        if (lhs.data[i] != rhs.data[i])
          return lhs.data[i] < rhs.data[i]? -1: 1;

      return 0;
    }
    template<typename S, typename C, C M, typename I, enable_scalar_operator_t<I> = nullptr>
    int abs_compare(mp_integer<S, C, M> const& lhs, I const& rhs) {
      using integer_t = mp_integer<S, C, M>;
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
    int impl_compare(mp_integer<S, C, M> const& lhs, T const& rhs) {
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
    int compare(mp_integer<S, C, M> const& lhs, mp_integer<S, C, M> const& rhs) {return impl_compare(lhs, rhs);}
    template<typename S, typename C, C M, typename I, enable_scalar_operator_t<I> = nullptr>
    int compare(mp_integer<S, C, M> const& lhs, I const& rhs) {return impl_compare(lhs, rhs);}
    template<typename S, typename C, C M, typename I, enable_scalar_operator_t<I> = nullptr>
    int compare(I const& lhs, mp_integer<S, C, M> const& rhs) {return -impl_compare(rhs, lhs);}
  }

  using mp_integer_detail::compare;

  //
  // a == b, a < b, etc.
  //
  namespace mp_integer_detail {
    template<typename S, typename C, C M>
    int operator==(mp_integer<S, C, M> const& lhs, mp_integer<S, C, M> const& rhs) {return compare(lhs, rhs) == 0;}
    template<typename S, typename C, C M, typename I, enable_scalar_operator_t<I> = nullptr>
    int operator==(mp_integer<S, C, M> const& lhs, I const& rhs) {return compare(lhs, rhs) == 0;}
    template<typename S, typename C, C M, typename I, enable_scalar_operator_t<I> = nullptr>
    int operator==(I const& lhs, mp_integer<S, C, M> const& rhs) {return compare(rhs, lhs) == 0;}
    template<typename S, typename C, C M>
    int operator!=(mp_integer<S, C, M> const& lhs, mp_integer<S, C, M> const& rhs) {return compare(lhs, rhs) != 0;}
    template<typename S, typename C, C M, typename I, enable_scalar_operator_t<I> = nullptr>
    int operator!=(mp_integer<S, C, M> const& lhs, I const& rhs) {return compare(lhs, rhs) != 0;}
    template<typename S, typename C, C M, typename I, enable_scalar_operator_t<I> = nullptr>
    int operator!=(I const& lhs, mp_integer<S, C, M> const& rhs) {return compare(rhs, lhs) != 0;}
    template<typename S, typename C, C M>
    int operator<(mp_integer<S, C, M> const& lhs, mp_integer<S, C, M> const& rhs) {return compare(lhs, rhs) < 0;}
    template<typename S, typename C, C M, typename I, enable_scalar_operator_t<I> = nullptr>
    int operator<(mp_integer<S, C, M> const& lhs, I const& rhs) {return compare(lhs, rhs) < 0;}
    template<typename S, typename C, C M, typename I, enable_scalar_operator_t<I> = nullptr>
    int operator<(I const& lhs, mp_integer<S, C, M> const& rhs) {return compare(rhs, lhs) < 0;}
    template<typename S, typename C, C M>
    int operator>(mp_integer<S, C, M> const& lhs, mp_integer<S, C, M> const& rhs) {return compare(lhs, rhs) > 0;}
    template<typename S, typename C, C M, typename I, enable_scalar_operator_t<I> = nullptr>
    int operator>(mp_integer<S, C, M> const& lhs, I const& rhs) {return compare(lhs, rhs) > 0;}
    template<typename S, typename C, C M, typename I, enable_scalar_operator_t<I> = nullptr>
    int operator>(I const& lhs, mp_integer<S, C, M> const& rhs) {return compare(rhs, lhs) > 0;}
    template<typename S, typename C, C M>
    int operator<=(mp_integer<S, C, M> const& lhs, mp_integer<S, C, M> const& rhs) {return compare(lhs, rhs) <= 0;}
    template<typename S, typename C, C M, typename I, enable_scalar_operator_t<I> = nullptr>
    int operator<=(mp_integer<S, C, M> const& lhs, I const& rhs) {return compare(lhs, rhs) <= 0;}
    template<typename S, typename C, C M, typename I, enable_scalar_operator_t<I> = nullptr>
    int operator<=(I const& lhs, mp_integer<S, C, M> const& rhs) {return compare(rhs, lhs) <= 0;}
    template<typename S, typename C, C M>
    int operator>=(mp_integer<S, C, M> const& lhs, mp_integer<S, C, M> const& rhs) {return compare(lhs, rhs) >= 0;}
    template<typename S, typename C, C M, typename I, enable_scalar_operator_t<I> = nullptr>
    int operator>=(mp_integer<S, C, M> const& lhs, I const& rhs) {return compare(lhs, rhs) >= 0;}
    template<typename S, typename C, C M, typename I, enable_scalar_operator_t<I> = nullptr>
    int operator>=(I const& lhs, mp_integer<S, C, M> const& rhs) {return compare(rhs, lhs) >= 0;}
  }
  using mp_integer_detail::operator==;
  using mp_integer_detail::operator!=;
  using mp_integer_detail::operator> ;
  using mp_integer_detail::operator< ;
  using mp_integer_detail::operator>=;
  using mp_integer_detail::operator<=;

  //
  // a += b, a -= b, a + b, a - b
  //
  namespace mp_integer_detail {

    template<typename S, typename C, C M>
    void abs_add(mp_integer<S, C, M>& lhs, mp_integer<S, C, M> const& rhs) {
      typedef mp_integer<S, C, M> integer_t;
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
    void abs_add(mp_integer<S, C, M>& lhs, I const& rhs) {
      using integer_t = mp_integer<S, C, M>;
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
    void abs_sub(mp_integer<S, C, M>& lhs, mp_integer<S, C, M> const& rhs) {
      typedef mp_integer<S, C, M> integer_t;
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
    void abs_sub(mp_integer<S, C, M>& lhs, I const& rhs) {
      typedef mp_integer<S, C, M> integer_t;
      typedef typename integer_t::element_type elem_t;
      typedef typename integer_t::calculation_type calc_t;

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
    mp_integer<S, C, M>& impl_add_eq(mp_integer<S, C, M>& lhs, T const& rhs) {
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
    mp_integer<S, C, M>& operator+=(mp_integer<S, C, M>& lhs, mp_integer<S, C, M> const& rhs) {return impl_add_eq<'+'>(lhs, rhs);}
    template<typename S, typename C, C M, typename I, enable_scalar_operator_t<I> = nullptr>
    mp_integer<S, C, M>& operator+=(mp_integer<S, C, M>& lhs, I const& rhs) {return impl_add_eq<'+'>(lhs, rhs);}
    template<typename S, typename C, C M>
    mp_integer<S, C, M>& operator-=(mp_integer<S, C, M>& lhs, mp_integer<S, C, M> const& rhs) {return impl_add_eq<'-'>(lhs, rhs);}
    template<typename S, typename C, C M, typename I, enable_scalar_operator_t<I> = nullptr>
    mp_integer<S, C, M>& operator-=(mp_integer<S, C, M>& lhs, I const& rhs) {return impl_add_eq<'-'>(lhs, rhs);}


    template<typename S, typename C, C M, typename T, enable_generic_operator_t<mp_integer<S, C, M>, T> = nullptr>
    mp_integer<S, C, M> operator+(mp_integer<S, C, M> const& lhs, T const& rhs) {
      mp_integer<S, C, M> ret;
      if (lhs.sign != 0) _set(ret, lhs);
      ret += rhs;
      return ret;
    }
    template<typename S, typename C, C M, typename I, enable_scalar_operator_t<I> = nullptr>
    mp_integer<S, C, M> operator+(I const& lhs, mp_integer<S, C, M> const& rhs) {
      mp_integer<S, C, M> ret;
      if (rhs.sign != 0) _set(ret, rhs);
      ret += lhs;
      return ret;
    }
    template<typename S, typename C, C M, typename T, enable_generic_operator_t<mp_integer<S, C, M>, T> = nullptr>
    mp_integer<S, C, M> operator-(mp_integer<S, C, M> const& lhs, T const& rhs) {
      mp_integer<S, C, M> ret;
      if (lhs.sign != 0) _set(ret, lhs);
      ret -= rhs;
      return ret;
    }
    template<typename S, typename C, C M, typename I, enable_scalar_operator_t<I> = nullptr>
    mp_integer<S, C, M> operator-(I const& lhs, mp_integer<S, C, M> const& rhs) {
      mp_integer<S, C, M> ret;
      if (rhs.sign != 0) {
        _set(ret, rhs);
        ret.sign = -ret.sign;
      }
      ret -= lhs;
      return ret;
    }
  }

  using mp_integer_detail::operator+=;
  using mp_integer_detail::operator-=;
  using mp_integer_detail::operator+;
  using mp_integer_detail::operator-;

  template<typename S, typename C, C M>
  void dump(mp_integer<S, C, M> const& value) {
    std::printf("sign = %d, data = [ ", value.sign);
    for (std::size_t i = value.data.size(); i--; ) {
      std::printf("%llx", (unsigned long long) value.data[i]);
      if (i != 0) std::printf(", ");
    }
    std::printf(" ]\n");
  }
}

namespace multi_precision_i1 {
  void test1() {
    mp_integer<std::uint32_t, std::uint64_t, 0x100> a = 1234;
    dump(a);
    for (int i = 0; i <= 100; i++) {
      a += a;
      a += 1;
      dump(a);
    }
  }
}

int main() {
  std::printf("hello\n");

  multi_precision_i1::test1();

  return 0;
}
