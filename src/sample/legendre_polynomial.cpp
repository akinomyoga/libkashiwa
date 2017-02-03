#include <iostream>
#include <cstdint>
#include <ksh/rational.h>
#include <ksh/polynomial.h>

int main() {
  using int_t = std::int64_t;
  using rational_t = kashiwa::rational<int_t>;
  using poly_t = kashiwa::polynomial<rational_t>;
  poly_t p1 = {1};
  poly_t p2 = {0, 1};
  std::cout << "P0 = " << p1 << std::endl;
  std::cout << "P1 = " << p2 << std::endl;

  auto const x = poly_t {0, 1};
  for (int n = 1; n < 30; n++) {
    auto p3 = (rational_t(2 * n + 1) * x * p2 - rational_t(n) * p1) / rational_t(n + 1);

    int_t lcm = 1;
    for (rational_t const& coeff: p3.data())
      lcm = kashiwa::lcm(lcm, coeff.denominator());
    std::cout << "P" << n + 1 << " = (1/" << lcm << ")(" << p3 * rational_t(lcm) << ")" << std::endl;

    p1 = p2;
    p2 = p3;
  }
  return 0;
}
