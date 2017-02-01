#include <ctime>
#include <cstdlib>
#include <algorithm>
#include <limits>
#include <mwg/except.h>
#include <ksh/integrator.h>

namespace kashiwa {
  template<typename T, std::size_t N>
  constexpr std::size_t size(T const (&)[N]) {return N;}
}

namespace ksh = kashiwa;

int main(){
  constexpr int KMAX = 10;
  constexpr int order = 256;

  double value[KMAX];
  //ksh::gauss_chebyshev_quadrature(KMAX, value, order, 0.0, 1.0, [=](double* integ, double t) {
  ksh::gauss_legendre_quadrature<order>(KMAX, value, 0.0, 1.0, [=](double* integ, double t) {
    double const x = 0.1 * (1.0 / t - 1.0);
    double const J = 0.1 / (t * t);
    double const n = 1.0 / std::expm1(x);
    for (int k = 0; k < KMAX; k++)
      integ[k] = J * n * std::pow(x, k + 1);
  });

  double const correct[] = {
    0.0, 6.0, 90.0, 945.0, 9450.0,
    93555.0, 638512875.0 / 691.0,
  };

  double fact = 1;
  for (int k = 0; k < KMAX; k++) {
    int const s = k + 2;
    fact *= k + 1;
    value[k] /= fact;

    double const rationalPart = std::pow(M_PI, s) / value[k];
    std::printf("zeta(%d) = %.15g = pi^%d / %.15g\n", s, value[k], s, rationalPart);
    if (s % 2 == 0) {
      double const expected = correct[s / 2];
      double const error = std::abs(rationalPart - expected) / std::max({std::abs(rationalPart), std::abs(expected), 1.0});
      double const thresh = order * 1e2 * std::numeric_limits<double>::epsilon();
      mwg_check_nothrow(
        error < thresh,
        "%g (expected: %g, error: %g > %g)", rationalPart, expected, error, thresh);
    }
  }


  return 0;
}
