// -*- mode: C++; coding: utf-8 -*-
#pragma once
#ifndef IDT_RFH_LINEAR_I1_H
#define IDT_RFH_LINEAR_I1_H
#include <cstdlib>
namespace idt{
namespace rfh{

  /// @fn void LUDecompose(int N, double* arr, int* imap);
  ///   正方行列 A を LU 分解して下三角行列 L と 上三角行列 U に分解します。
  ///   @param[in    ] N    行列の次元を指定します。
  ///   @param[in,out] arr  入力の行列 \f$A\f$ を指定します。
  ///     結果として得られた下三角行列と上三角行列の情報を格納します。
  ///     下三角行列 L_{ij} = 
  ///   @param[   out] imap pivot 選択による行の順序を格納します。
  void LUDecompose(int N, double* arr, int* imap);
  /// @fn void LUEquationSolve(int N, double const* arr, int const* imap, double* vec);
  void LUEquationSolve(int N, double const* arr, int const* imap, double* vec);
  /// @fn void SolveLinearEquationLU(std::size_t N, double* arr, double* vec);
  void SolveLinearEquationLU(std::size_t N, double* arr, double* vec);

  /// @fn void SolveLinearEquation(int N, double* arr, double* vec);
  ///   連立線形方程式 \f$A x = b\f$ を LU 分解によって解きます。
  ///   @param[in    ] N   解く方程式の次元を指定します。
  ///   @param[in,out] arr 方程式の係数行列 \f$A\f$ を指定します。
  ///     <del>LU 分解した上三角行列と下三角行列の情報を格納します。</del>内部的に使用します。
  ///     \f$A_{ij}\f$ は arr[i*N+j] に対応します。
  ///   @param[in,out] vec 方程式の定数項のベクトル \f$b\f$ を指定します。
  ///     線形方程式の根 \f$x\f$ を格納します。
  void SolveLinearEquation(int N, double* arr, double* vec);

}
}

#endif
