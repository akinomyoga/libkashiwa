#include <iostream>
#include <cstdint>
#include <ksh/rational.hpp>
#include <ksh/polynomial.hpp>
#include <ksh/big_integer.hpp>
#include <mwg/except.h>

int main() {
  using int_t = kashiwa::bigint;
  //using int_t = std::int64_t;
  using rational_t = kashiwa::rational<int_t>;
  using poly_t = kashiwa::polynomial<rational_t>;
  poly_t p1 = {1};
  poly_t p2 = {0, 1};
  std::cout << "\x1b[1mP0\x1b[m = " << p1 << std::endl;
  std::cout << "\x1b[1mP1\x1b[m = " << p2 << std::endl;

  auto const x = poly_t {0, 1};
  for (int n = 1; n < 30; n++) {
    auto p3 = ((2 * n + 1) * x * p2 - n * p1) / (n + 1);

    int_t lcm = 1;
    for (rational_t const& coeff: p3.data())
      lcm = kashiwa::lcm(lcm, coeff.denominator());
    std::cout << "\x1b[1mP" << n + 1 << "\x1b[m = (1/" << lcm << ")(" << p3 * lcm << ")" << std::endl;

    p1 = p2;
    p2 = p3;
  }
  return 0;
}
