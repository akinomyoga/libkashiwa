// -*- mode: c++; coding: utf-8 -*-
#ifndef KASHIWA_LINEAR_FIT_H
#define KASHIWA_LINEAR_FIT_H
#include <cstdlib>
#include <numeric>
#include "buffer.h"
#include "linear_lu.h"
namespace kashiwa {

  inline void linear_fit(int NParam, double* param, int NData, double const* data, double const* basis, working_buffer& buffer) {
    buffer.ensure<double>(NParam + NParam * NParam);
    double* const contra     = buffer.ptr<double>();
    double* const covariance = contra + NParam;

    for (int i = 0; i < NParam; i++) {
      double const* const fi = basis + i * NData;
      contra[i] = std::inner_product(fi, fi + NData, data, 0.0);
      for (int j = 0; j <= i; j++) {
        double const* const fj = basis + j * NData;
        double const value = std::inner_product(fi, fi + NData, fj, 0.0);
        covariance[i * NParam + j] = value;
        covariance[j * NParam + i] = value;
      }
    }

    working_buffer buffer2;
    kashiwa::solve_linear_equation_lu(NParam, param, covariance, contra, buffer2, covariance);
  }

}

#endif
