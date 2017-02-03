#include <cstdio>
#include <cstdint>
#include <limits>
#include <vector>
#include <algorithm>
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
      Modulo > 1 &&
      Modulo - 1 <= std::numeric_limits<element_type>::max() &&
      element_max <= std::numeric_limits<calculation_type>::max() / element_max);

    int sign;
    std::vector<element_type> data;

    mp_integer(): sign(0) {}

    mp_integer(int value) {
      if (value == 0) {
        sign = 0;
        return;
      }

      unsigned int uvalue = value;
      if (value < 0) {
        sign = -1;
        uvalue = -uvalue;
      } else {
        sign = 1;
      }

      while (uvalue) {
        data.push_back(uvalue % modulo);
        uvalue = uvalue / modulo;
      }

      std::reverse(data.begin(), data.end());
    }
  };

  template<typename S, typename C, C M>
  void abs_set(mp_integer<S, C, M>& lhs, mp_integer<S, C, M> const& rhs) {
    std::vector<typename mp_integer<S, C, M>::element_type>& ddata = lhs.data;
    ddata.resize(0);
    ddata.reserve(rhs.data.size());
    ddata.insert(ddata.begin(), rhs.data.begin(), rhs.data.end());
  }

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

  template<typename S, typename C, C M>
  int abs_compare(mp_integer<S, C, M> const& lhs, mp_integer<S, C, M> const& rhs) {
    if (lhs.data.size() != rhs.data.size())
      return lhs.data.size() < rhs.data.size()? -1: 1;

    for (std::size_t i = lhs.data.size(); i--; )
      if (lhs.data[i] != rhs.data[i])
        return lhs.data[i] < rhs.data[i]? -1: 1;

    return 0;
  }

  template<typename S, typename C, C M>
  void abs_sub(mp_integer<S, C, M>& lhs, mp_integer<S, C, M> const& rhs) {
    typedef mp_integer<S, C, M> integer_t;
    typedef typename integer_t::element_type element_t;
    typedef typename integer_t::calculation_type calc_t;

    int const cmp = abs_compare(lhs, rhs);
    if (cmp == 0) {
      lhs.sign = 0;
      lhs.data.resize(1, (element_t) 0);
    } else if (cmp != 0) {
      element_t carry = 0;
      element_t      * dst = &lhs.data[0];
      element_t const* max;
      element_t const* min;
      std::size_t maxN;
      std::size_t minN;

      if (cmp < 0) {
        lhs.sign = -lhs.sign;
        lhs.data.resize(rhs.data.size(), (element_t) 0);
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
          dst[i] = (element_t) (max[i] + (integer_t::modulo - sub));
          carry = 1;
        } else {
          dst[i] = (element_t) (max[i] - sub);
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
        lhs.data.resize(highest + 1, (element_t) 0);
      }
    }
  }

  template<typename S, typename C, C M>
  mp_integer<S, C, M>& operator+=(mp_integer<S, C, M>& lhs, mp_integer<S, C, M> const& rhs) {
    if (rhs.sign == 0) return lhs;
    if (lhs.sign == 0) {
      lhs.sign = -rhs.sign;
      abs_set(lhs, rhs);
      return lhs;
    }

    if (lhs.sign == rhs.sign)
      abs_add(lhs, rhs);
    else
      abs_sub(lhs, rhs);
    return lhs;
  }

  template<typename S, typename C, C M>
  mp_integer<S, C, M> operator+(mp_integer<S, C, M> const& lhs, mp_integer<S, C, M> const& rhs) {
    mp_integer<S, C, M> ret;
    if (lhs.sign != 0)
      abs_set(ret, lhs);
    ret += rhs;
    return ret;
  }

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
      a += (decltype(a)) 1;
      dump(a);
    }
  }
}

int main() {
  std::printf("hello\n");

  multi_precision_i1::test1();

  return 0;
}
