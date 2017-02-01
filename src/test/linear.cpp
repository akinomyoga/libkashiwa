#include <ctime>
#include <cstdlib>
#include <algorithm>
#include <limits>
#include <mwg/except.h>
#include <ksh/linear_lu.h>

namespace kashiwa {
  template<typename T, std::size_t N>
  constexpr std::size_t size(T const (&)[N]) {return N;}
}

namespace ksh = kashiwa;

int main(){
  std::srand(std::time(NULL));

  constexpr std::size_t N = 100;

  double mat1[N*N];
  std::generate(std::begin(mat1), std::end(mat1), [](){return (std::rand() + 0.5) / (RAND_MAX + 1.0);});

  double vec1[N];
  std::generate(std::begin(vec1), std::end(vec1), [](){return (std::rand() + 0.5) / (RAND_MAX + 1.0);});

  double vec2[N];
  for (std::size_t i = 0; i < N; i++) {
    double sum = 0.0;
    for (std::size_t j = 0; j < N; j++)
      sum += mat1[i * N + j] * vec1[j];
    vec2[i] = sum;
  }

  ksh::working_buffer buffer;
  ksh::solve_linear_equation_lu(N, vec2, mat1, vec2, buffer);

  constexpr double tol = N * 100.0 * std::numeric_limits<double>::epsilon();
  for (std::size_t i = 0; i < N; i++) {
    mwg_check(
      std::abs(vec1[i] - vec2[i]) < tol,
      "error = %g (tol = %g)", std::abs(vec1[i] - vec2[i]), tol);
  }

  return 0;
}
