// -*- mode: C++; coding: utf-8 -*-
#ifndef kashiwa_linear_lu_hpp
#define kashiwa_linear_lu_hpp
#include <cstdlib>
#include "buffer.hpp"
namespace kashiwa {

  /*?lwiki @fn void lu_decompose(int N, double* lumat, int* lupiv);
   *   正方行列 $A$ を LU 分解して下三角行列 $L$ と 上三角行列 $U$ に分解します。
   *   $L$ の対角成分は 1 に規格化されます。
   *
   *   @param[in] std::size_t N;
   *     行列の次元を指定します。
   *
   *   @param[in,out] double* lumat;
   *     行列 $A$ を指定します。
   *     結果として得られた下三角行列 $L$ と上三角行列 $U$ の情報を格納します。
   *     下三角行列は $L_{ij} = \texttt{lumat[lupiv[i] * N + j]}, \quad(j < i)$ で与えられます。
   *     上三角行列は $U_{ij} = \texttt{lumat[lupiv[i] * N + j]}, \quad(j \ge i)$ で与えられます。
   *
   *   @param[out] int* lupiv;
   *     pivot 選択による行の順序を格納します。
   *
   */
  void lu_decompose(std::size_t N, double* lumat, int* lupiv);

  /*?lwiki
   * @fn void solve_lu_equation(std::size_t N, double* result, double const* lumat, int const* lupiv, double const* vec, working_buffer& buffer);
   *   LU 分解された行列を用いて線形方程式 $LUx = b$ を解きます。
   *
   *   @param[in] std::size_t N;
   *     行列の次元を指定します。
   *
   *   @param[out] double* result;
   *     解 $x$ の格納先を指定します。
   *
   *   @param[in] double const* lumat;
   *   @param[in] int const* lupiv;
   *     `lu_decompose` によって LU 分解された行列データを指定します。
   *
   *   @param[in] double const* vec;
   *     定数項 $b$ を指定します。
   *
   *   @param[in] working_buffer& buffer;
   *     作業用領域として用いるバッファを指定します。
   *
   */
  void solve_lu_equation(std::size_t N, double* result, double const* lumat, int const* lupiv, double const* vec, working_buffer& buffer);

  /*?lwiki
   * @fn void solve_linear_equation_lu(std::size_t N, double* result, double const* mat, double const* vec, working_buffer& buffer, double* tmpMat = nullptr);
   *   連立線形方程式 $A x = b$ を LU 分解によって解きます。
   *
   *   @param[in] std::size_t N;
   *     方程式の次元を指定します。
   *
   *   @param[out] double* result;
   *     解 $x$ の格納先を指定します。
   *
   *   @param[in] double const* mat;
   *     方程式の係数行列 $A$ を指定します。$A_{ij}$ は `mat[i * N + j]` に対応します。
   *
   *   @param[in] double const* vec;
   *     定数項 $b$ を指定します。
   *
   *   @param[in] working_buffer& buffer;
   *     作業用領域として用いるバッファを指定します。
   *
   *   @param[out,opt] double* tmpMat;
   *     作業用行列を指定します。`N * N` のサイズが必要です。
   *     `nullptr` を指定すると作業用配列は `buffer` 内に確保します。
   *     `mat` を指定すると `mat` を破壊的に使用します。
   *     それ以外を指定した場合は `mat` との重複はないものとして扱います。
   *
   */
  void solve_linear_equation_lu(std::size_t N, double* result, double const* mat, double const* vec, working_buffer& buffer, double* tmpMat = nullptr);
}

#endif
