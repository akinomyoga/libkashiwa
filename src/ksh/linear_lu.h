// -*- mode: C++; coding: utf-8 -*-
#ifndef KASHIWA_LINEAR_LU_H
#define KASHIWA_LINEAR_LU_H
#include <cstdlib>
namespace kashiwa {

  /// @fn void lu_decompose(int N, double* arr, int* imap);
  ///   正方行列 A を LU 分解して下三角行列 L と 上三角行列 U に分解します。
  ///   @param[in    ] N    行列の次元を指定します。
  ///   @param[in,out] arr  入力の行列 \f$A\f$ を指定します。
  ///     結果として得られた下三角行列と上三角行列の情報を格納します。
  ///     下三角行列 L_{ij} =
  ///   @param[   out] imap pivot 選択による行の順序を格納します。
  void lu_decompose(std::size_t N, double* arr, int* imap);
  /// @fn void solve_lu_equation(int N, double const* luMatrix, int const* imap, double* vec);
  void solve_lu_equation(std::size_t N, double const* luMatrix, int const* imap, double* vec);
  /// @fn void solve_linear_equation_lu(std::size_t N, double* arr, double* vec);
  void solve_linear_equation_lu(std::size_t N, double* arr, double* vec);

  /*?lwiki
   * @fn void solve_linear_equation(std::size_t N, double* mat, double* vec, double* tmpMat);
   *   連立線形方程式 $A x = b$ を LU 分解によって解きます。
   *   @param[in    ] N
   *     解く方程式の次元を指定します。
   *   @param[in    ] arr
   *     方程式の係数行列 $A$ を指定します。
   *     <del>LU 分解した上三角行列と下三角行列の情報を格納します。</del>内部的に使用します。
   *     $A_{ij}$ は `arr[i*N+j]` に対応します。
   *   @param[in,out] vec
   *     方程式の定数項のベクトル \f$b\f$ を指定します。
   *     線形方程式の根 \f$x\f$ を格納します。
   *   @param[   out] tmpMat
   *     内部的に使用するバッファを指定します。少なくとも `N*N` のサイズが必要です。
   */
  void solve_linear_equation(std::size_t N, double const* mat, double* vec, double* tmpMat);

}

#endif
