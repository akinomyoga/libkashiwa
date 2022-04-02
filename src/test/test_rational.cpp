#include <cstdio>
#include <sstream>
#include <mwg/except.h>
#include <ksh/rational.hpp>

namespace ksh = kashiwa;

using rational_t = ksh::rational<int>;

void check_print(rational_t const& value, const char* expected) {
  std::ostringstream str;
  str << value;
  mwg_check(str.str() == expected, "result = \"%s\" expected = \"%s\"", str.str().c_str(), expected);
}


int main() {

  // test initialization
  {
    constexpr rational_t r1 {242, 484};
    mwg_check((r1.numerator() == 1 && r1.denominator() == 2));
    constexpr rational_t r2 {0, 4};
    mwg_check((r2.numerator() == 0 && r2.denominator() == 1));
    constexpr rational_t r3 {3, 0};
    mwg_check((r3.numerator() == 1 && r3.denominator() == 0));
    constexpr rational_t r4 {1, -2};
    mwg_check((r4.numerator() == -1 && r4.denominator() == 2));
  }

  constexpr rational_t zero;
  constexpr rational_t nan {0, 0};
  constexpr rational_t inf {1, 0};
  mwg_check((zero == 0));

  // ostr << rational
  {
    check_print(zero, "0");
    check_print(nan, "nan");
    check_print(inf, "inf");
    check_print(-inf, "-inf");
    check_print(rational_t {1, 2}, "(1/2)");
    check_print(rational_t {2, 1}, "2");
    check_print(rational_t {2, -1}, "-2");
    check_print(rational_t {2, -3}, "(-2/3)");
  }

  // +rational
  {
    constexpr rational_t a = -rational_t {1, 2};
    constexpr rational_t b = +rational_t {1, 2};
    mwg_check((a == rational_t {-1, 2}));
    mwg_check((b == rational_t {1, 2}));

    constexpr rational_t pinf {1, 0};
    constexpr rational_t ninf {-1, 0};
    mwg_check((pinf == inf));
    mwg_check((ninf == -inf));

    mwg_check((+zero == 0));
    mwg_check((-zero == 0));
    mwg_check((+pinf == inf));
    mwg_check((-pinf == -inf));
    mwg_check((+ninf == -inf));
    mwg_check((-ninf == inf));
    mwg_check((+nan == nan));
    mwg_check((-nan == nan));
  }

  // rational + rational
  {
    constexpr rational_t c = rational_t {1, 2} + rational_t {1, 3};
    constexpr rational_t d = rational_t {1, 2} - rational_t {1, 3};
    mwg_check((c == rational_t {5, 6}));
    mwg_check((d == rational_t {1, 6}));

    mwg_check((zero + zero == zero));
    mwg_check((zero - zero == zero));
    mwg_check((zero +  nan == nan));
    mwg_check((zero -  nan == nan));
    mwg_check((zero +  inf == inf));
    mwg_check((zero -  inf == -inf));
    mwg_check(( nan +  inf == nan));
    mwg_check(( nan -  inf == nan));
    mwg_check(( inf +  inf == inf));
    mwg_check(( inf - -inf == inf));
    mwg_check((-inf + -inf == -inf));
    mwg_check((-inf -  inf == -inf));
    mwg_check(( inf + -inf == nan));
    mwg_check(( inf -  inf == nan));
    mwg_check((-inf +  inf == nan));
    mwg_check((-inf - -inf == nan));
  }

  // rational * rational
  {
    constexpr rational_t e = rational_t {1, 2} * rational_t {2, 3};
    constexpr rational_t f = rational_t {1, 2} / rational_t {2, 3};
    mwg_check((e == rational_t {1, 3}));
    mwg_check((f == rational_t {3, 4}));

    mwg_check((zero * zero == 0));
    mwg_check((zero *  inf == nan));
    mwg_check((zero * -inf == nan));
    mwg_check((zero / zero == nan));
    mwg_check((zero /  inf == 0));
    mwg_check((zero / -inf == 0));
    mwg_check((zero *  nan == nan));
    mwg_check((zero /  nan == nan));
    mwg_check(( inf / zero ==  inf));
    mwg_check((-inf / zero == -inf));
    mwg_check(( inf *  inf ==  inf));
    mwg_check(( inf * -inf == -inf));
    mwg_check((-inf * -inf ==  inf));
    mwg_check(( inf /  inf == nan));
    mwg_check(( inf / -inf == nan));
    mwg_check((-inf / -inf == nan));
    mwg_check(( inf *  nan == nan));
    mwg_check((-inf *  nan == nan));
    mwg_check(( inf /  nan == nan));
    mwg_check((-inf /  nan == nan));
    mwg_check(( nan / zero == nan));
    mwg_check(( nan /  inf == nan));
    mwg_check(( nan / -inf == nan));
    mwg_check(( nan *  nan == nan));
    mwg_check(( nan /  nan == nan));
  }

  // rational * scalar
  {
    constexpr rational_t r1 = 1 / rational_t{3, 2};
    mwg_check((r1 == rational_t {2, 3}));

    mwg_check((zero * 0 == 0));
    mwg_check((zero * 2 == 0));
    mwg_check((zero * -5 == 0));
    mwg_check((nan * 2 == nan));
    mwg_check((nan * -5 == nan));
    mwg_check((nan * 0 == nan));
    mwg_check((inf * 2 == inf));
    mwg_check((inf * -5 == -inf));
    mwg_check((inf * 0 == nan));

    mwg_check((zero / 0 == nan));
    mwg_check((zero / 2 == 0));
    mwg_check((zero / -5 == 0));
    mwg_check((nan / 2 == nan));
    mwg_check((nan / -5 == nan));
    mwg_check((nan / 0 == nan));
    mwg_check((inf / 2 == inf));
    mwg_check((inf / -5 == -inf));
    mwg_check((inf / 0 == inf));

    mwg_check((0 / inf == 0));
    mwg_check((5 / inf == 0));
    mwg_check((-5 / inf == 0));
    mwg_check((0 / nan == nan));
    mwg_check((5 / nan == nan));
    mwg_check((-5 / nan == nan));
    mwg_check((0 / zero == nan));
    mwg_check((5 / zero == inf));
    mwg_check((-5 / zero == -inf));
  }

  // rational += rational, etc.
  {
    rational_t r1 {1, 3};
    r1 += rational_t {1, 2};
    mwg_check((r1 == rational_t {5, 6}));
    rational_t r2 {1, 2};
    r2 -= rational_t {1, 3};
    mwg_check((r2 == rational_t {1, 6}));
    rational_t r3 {1, 2};
    r3 *= rational_t {1, 3};
    mwg_check((r3 == rational_t {1, 6}));
    rational_t r4 {1, 2};
    r4 /= rational_t {3, 2};
    mwg_check((r4 == rational_t {1, 3}));
  }

  return 0;
}
