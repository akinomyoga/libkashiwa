#include <cstdio>
#include <cstdint>
#include <vector>
#include <algorithm>
#include <mwg/except.h>

namespace multi_precision_i1 {

  struct mp_integer {
    typedef std::uint32_t element_type;
    static constexpr element_type modulo = 0x10000;

    int sign;
    std::vector<element_type> data;

    mp_integer(): sign(0) {}
    mp_integer(int value): sign(0) {
      if (value == 0) return;

      unsigned int uvalue = value;
      if (value < 0) {
        sign = -1;
        uvalue = - uvalue;
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

  void abs_set(mp_integer& lhs, mp_integer const& rhs) {
    std::vector<mp_integer::element_type>& ddata = lhs.data;
    ddata.resize(0);
    ddata.reserve(rhs.data.size());
    ddata.insert(ddata.begin(), rhs.data.begin(), rhs.data.end());
  }

  void abs_add(mp_integer& lhs, mp_integer const& rhs) {
    if (lhs.data.size() < rhs.data.size())
      lhs.data.resize(rhs.data.size(), (mp_integer::element_type) 0u);

    mp_integer::element_type carry = 0;
    std::size_t const iN = rhs.data.size();
    for (std::size_t i = 0; i < iN; i++) {
      mp_integer::element_type const value = lhs.data[i] + rhs.data[i] + carry;
      carry       = value / mp_integer::modulo;
      lhs.data[i] = value % mp_integer::modulo;
    }

    if (carry) {
      mwg_assert(iN <= lhs.data.size());
      if (iN < lhs.data.size())
        lhs.data[iN] = carry;
      else
        lhs.data.push_back(carry);
    }
  }

  int abs_compare(mp_integer const& lhs, mp_integer const& rhs) {
    if (lhs.data.size() != rhs.data.size())
      return lhs.data.size() < rhs.data.size()? -1: 1;

    for (std::size_t i = lhs.data.size(); i--; )
      if (lhs.data[i] != rhs.data[i])
        return lhs.data[i] < rhs.data[i]? -1: 1;

    return 0;
  }

  void abs_sub(mp_integer& lhs, mp_integer const& rhs) {
    int const cmp = abs_compare(lhs, rhs);
    if (cmp == 0) {
      lhs.sign = 0;
      lhs.data.resize(1, (mp_integer::element_type) 0);
    } else if (cmp != 0) {
      mp_integer::element_type carry = 0;
      mp_integer::element_type      * dst = &lhs.data[0];
      mp_integer::element_type const* max;
      mp_integer::element_type const* min;
      std::size_t maxN;
      std::size_t minN;

      if (cmp < 0) {
        lhs.sign = -lhs.sign;
        lhs.data.resize(rhs.data.size(), (mp_integer::element_type) 0);
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
        mp_integer::element_type const sub = min[i] + carry;
        if (max[i] < sub) {
          dst[i] = max[i] + (mp_integer::modulo - sub);
          carry = 1;
        } else {
          dst[i] = max[i] - sub;
          carry = 0;
        }
      }

      std::size_t highest = minN - 1;
      if (dst != max || carry) {
        highest = maxN - 1;
        for (std::size_t j = minN; j < maxN; j++) {
          if (max[j] < carry) {
            dst[j] = max[j] + (mp_integer::modulo - carry);
            carry = 1;
          } else {
            dst[j] = max[j] - carry;
            if (dst == max) {highest = j; break;}
            carry = 0;
          }
        }
      }

      if (highest + 1 == lhs.data.size()) {
        while (highest > 0 && lhs.data[highest] == 0) highest--;
        lhs.data.resize(highest + 1, (mp_integer::element_type) 0);
      }
    }
  }

  mp_integer& operator+=(mp_integer& lhs, mp_integer const& rhs) {
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

  mp_integer operator+(mp_integer const& lhs, mp_integer const& rhs) {
    mp_integer ret;
    if (lhs.sign != 0)
      abs_set(ret, lhs);
    ret += rhs;
    return ret;
  }

  void dump(mp_integer const& value) {
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
    mp_integer a = 1234;
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
