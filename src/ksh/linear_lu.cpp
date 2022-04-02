#include <cstdio>
#include <cmath>
#include <algorithm>
#include <utility>
#include <vector>
#include <iterator>
#include <mwg/except.h>
#include "buffer.hpp"
#include "linear_lu.hpp"

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

  class lu_solver {
    std::size_t N;
    double* result;
    double const* lumat;
    int    const* lupiv;
    double const* vec;
    double* vtmp;

    double A(int i, int j) const {return lumat[lupiv[i] * N + j];}
    double b(int i) const {return vec[lupiv[i]];}

    void forward_substitute() const {
      for (std::size_t i = 0; i < N; i++) {
        double value = b(i);
        for (std::size_t j = 0; j < i; j++)
          value -= A(i, j) * vtmp[j];
        vtmp[i] = value;
      }
    }
    void backward_substitute() const {
      for (std::size_t i = N; i--; ) {
        double value = vtmp[i];
        for (std::size_t j = i + 1; j < N; j++)
          value -= A(i, j) * result[j];
        result[i] = value / A(i, i);
      }
    }

  public:
    void solve(double* result, double const* lumat, int const* lupiv, double const* vec, double* vtmp) {
      mwg_assert(vec != vtmp);
      this->result = result;
      this->lumat  = lumat;
      this->lupiv  = lupiv;
      this->vec    = vec;
      this->vtmp   = vtmp;
      this->forward_substitute();
      this->backward_substitute();
    }

  public:
    lu_solver(std::size_t N): N(N) {}
  };

}

namespace kashiwa {

  void lu_decompose(std::size_t N, double* lumat, int* imap) {
    lu_decomposer(N).decompose(lumat, imap);
  }

  void solve_lu_equation(std::size_t N, double* result, double const* lumat, int const* lupiv, double const* vec, working_buffer& buffer) {
    double* vtmp = result;
    if (vtmp == vec) {
      buffer.ensure<double>(N);
      vtmp = buffer.ptr<double>();
    }

    lu_solver(N).solve(result, lumat, lupiv, vec, vtmp);
  }

  void solve_linear_equation_lu(std::size_t N, double* result, double const* mat, double const* vec, working_buffer& buffer, double* tmpMat) {
    {
      std::size_t sz = N * sizeof(int);
      if (!tmpMat) sz += N * N * sizeof(double);
      if (result == vec) sz += N * sizeof(double);
      buffer.ensure(sz);
    }

    double* pbuf = buffer.ptr<double>();
    double* lumat = tmpMat;
    if (!tmpMat) {
      lumat = pbuf;
      pbuf += N * N;
    }
    double* vtmp = result;
    if (result == vec) {
      vtmp = pbuf;
      pbuf += N;
    }
    int* const lupiv = reinterpret_cast<int*>(pbuf);

    if (lumat != mat) std::copy(mat, mat + N * N, lumat);
    lu_decompose(N, lumat, lupiv);

    lu_solver(N).solve(result, lumat, lupiv, vec, vtmp);
  }

}
