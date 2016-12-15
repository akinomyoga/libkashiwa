#include <cstdio>
#include <cmath>
#include <algorithm>
#include <utility>
#include <vector>
#include <iterator>
//#include <mwg/except.h>
#include "linear_lu.h"

namespace {

  class lu_decomposer {
    std::size_t N;
    double* arr;
    int* imap;

    double& A(int i, int j) {
      return arr[imap[i] * N + j];
    }

    double determine_pivot(int i) {
      int imax = i;
      double vmax = std::abs(A(i, i));
      for (std::size_t icand = i + 1; icand < N; icand++) {
        double const vcand = std::abs(A(icand, i));
        if (vcand > vmax) {
          imax = icand;
          vmax = vcand;
        }
      }
      if (imax > i) {
        using namespace std;
        swap(imap[i], imap[imax]);
      }
      //mwg_assert_nothrow(vmax!=0);
      return A(i, i); // pivot 交換後の対角成分
    }

  public:
    void decompose(double* arr, int* imap) {
      this->arr = arr;
      this->imap = imap;

      for (std::size_t i = 0; i < N; i++) imap[i] = i;

      for (std::size_t i = 0; i < N; i++) {
        double const scal = 1.0 / this->determine_pivot(i);

        for (std::size_t ii = i + 1; ii < N; ii++) {
          double const l = A(ii, i) *= scal;
          for (std::size_t jj = i + 1; jj < N; jj++) {
            double const u = A(i, jj);
            A(ii, jj) -= l * u;
          }
        }
      }
    }

  public:
    lu_decomposer(std::size_t N): N(N) {}
  };

  class lu_equation_solver {
    std::size_t N;
    double const* arr;
    double* vec;
    int const* imap;
    std::vector<double> vtmp;

    double const& A(int i, int j) {
      return arr[imap[i] * N + j];
    }
    double& b(int i) {
      return vec[imap[i]];
    }

    void forward_substitute() {
      for (std::size_t i = 0; i < N; i++) {
        double& vv(vtmp[i] = b(i));
        for (std::size_t j = 0; j < i; j++)
          vv -= A(i, j) * vtmp[j];
      }
    }
    void backward_substitute() {
      for (std::size_t i = N; i--; ) {
        double& vv(vec[i] = vtmp[i]);
        for (std::size_t j = i + 1; j < N; j++)
          vv -= A(i, j) * vec[j];
        vv /= A(i, i);
      }
    }

  public:
    void solve(double const* arr, int const* imap, double* vec) {
      this->arr = const_cast<double*>(arr);
      this->imap = const_cast<int*>(imap);
      this->vec = vec;
      this->forward_substitute();
      this->backward_substitute();
    }

  public:
    lu_equation_solver(std::size_t N): N(N) {
      this->vtmp.resize(N, 0.0);
    }
  };


  // template<int N>
  // void solve_linear_equation_lu(double* arr, double* vec) {
  //   // 取り敢えずの実装:
  //   //   arr を破壊的に使用する
  //   //   vec ベクトルを指定し結果を格納する。
  //   LinearEquationSolver_LU(N, arr, vec);
  // }

}

namespace kashiwa {
  void lu_decompose(std::size_t N, double* arr, int* imap) {
    lu_decomposer calculator(N);
    calculator.decompose(arr, imap);
  }
  void solve_lu_equation(std::size_t N, double const* luMatrix, int const* imap, double* vec) {
    lu_equation_solver calculator(N);
    calculator.solve(luMatrix, imap, vec);
  }
  void solve_linear_equation_lu(std::size_t N, double* arr, double* vec) {
    std::vector<int> imap((std::size_t) N);
    lu_decompose(N, arr, &imap[0]);
    solve_lu_equation(N, arr, &imap[0], vec);
  }
  void solve_linear_equation(std::size_t N, double const* mat, double* vec, double* tmpMat) {
    if (tmpMat != mat) std::copy(mat, mat + N * N, tmpMat);
    solve_linear_equation_lu(N, tmpMat, vec);
  }
}
