#include <cstdio>
#include <cstdint>
#include <iostream>
#include <sstream>
#include <mwg/except.h>
#include <ksh/polynomial.h>
#include <ksh/rational.h>

namespace ksh = kashiwa;

using polynomial_t = ksh::polynomial<int>;

void check_print(polynomial_t const& poly, const char* expected) {
  std::ostringstream str;
  str << poly;
  mwg_check(str.str() == expected, "result = \"%s\" expected = \"%s\"", str.str().c_str(), expected);
}

int main() {

  // comparison operator

  mwg_check((polynomial_t {} == 0));
  mwg_check((polynomial_t {1} == 1));
  mwg_check((polynomial_t {} != 1));
  mwg_check((polynomial_t {1} != 0));
  mwg_check((ksh::polynomial<ksh::rational<int>> {} == 0));
  mwg_check((ksh::polynomial<ksh::rational<int>> {1} == 1));
  mwg_check((ksh::polynomial<ksh::rational<int>> {} != 1));
  mwg_check((ksh::polynomial<ksh::rational<int>> {1} != 0));
  mwg_check((ksh::polynomial<ksh::rational<std::uint16_t>> {} == 0));
  mwg_check((ksh::polynomial<ksh::rational<std::uint16_t>> {1} == 1));
  mwg_check((ksh::polynomial<ksh::rational<std::uint16_t>> {} != 1));
  mwg_check((ksh::polynomial<ksh::rational<std::uint16_t>> {1} != 0));

  // initialization

  mwg_check((polynomial_t {0, 0} == polynomial_t {}));
  mwg_check((polynomial_t {1, 0} == polynomial_t {1}));

  mwg_check((deg(polynomial_t {}) == 0));
  mwg_check((deg(polynomial_t {0}) == 0));
  mwg_check((deg(polynomial_t {1}) == 0));
  mwg_check((deg(polynomial_t {0, 1}) == 1));
  mwg_check((deg(polynomial_t {0, 1, 0}) == 1));
  mwg_check((deg(polynomial_t {2, 1, 1}) == 2));

  // operator: -polynomial
  mwg_check((-polynomial_t {1, 2, 1} == polynomial_t {-1, -2, -1}));
  mwg_check((+polynomial_t {1, 2, 1} == polynomial_t {1, 2, 1}));

  // operator: polynomial + polynomial
  polynomial_t p1 {1, 1};
  mwg_check((p1 * p1 == polynomial_t {1, 2, 1}));
  mwg_check((p1 * p1 + p1 == polynomial_t {2, 3, 1}));
  mwg_check((p1 * p1 - p1 == polynomial_t {0, 1, 1}));
  mwg_check((p1 - p1 == polynomial_t()));

  // operator: polynomial + scalar
  mwg_check((p1 + 1 == polynomial_t {2, 1}));
  mwg_check((polynomial_t {2} - 2 == 0));


  polynomial_t p2 = p1 * p1;
  mwg_check(((p2 -= p1) == polynomial_t {0, 1, 1}));
  mwg_check(((p2 += p1) == polynomial_t {1, 2, 1}));
  mwg_check(((p2 *= p1) == polynomial_t {1, 3, 3, 1}));

  polynomial_t p3 {1, 1};
  mwg_check((pow(p3, 0) == polynomial_t {1}));
  mwg_check((pow(p3, 1) == polynomial_t {1, 1}));
  mwg_check((pow(p3, 2) == polynomial_t {1, 2, 1}));
  mwg_check((pow(p3, 3) == polynomial_t {1, 3, 3, 1}));
  mwg_check((pow(p3, 4) == polynomial_t {1, 4, 6, 4, 1}));
  mwg_check((pow(p3, 5) == polynomial_t {1, 5, 10, 10, 5, 1}));
  mwg_check((pow(p3, 6) == polynomial_t {1, 6, 15, 20, 15, 6, 1}));
  mwg_check((pow(p3, 7) == polynomial_t {1, 7, 21, 35, 35, 21, 7, 1}));

  check_print(p2, "x^3 + 3x^2 + 3x + 1");
  check_print(polynomial_t {1, 1, 1} , "x^2 + x + 1");
  check_print(polynomial_t {0, -1, 2}, "2x^2 - x");
  check_print(polynomial_t {1, -1}   , "-x + 1");
  check_print(polynomial_t {3}       , "3");

  return 0;
}
