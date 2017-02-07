// -*- mode:c++ -*-
#ifndef KASHIWA_EXPLICIT_RUNGE_KUTTA_H
#define KASHIWA_EXPLICIT_RUNGE_KUTTA_H
#ifdef _MSC_VER
# define _USE_MATH_DEFINES
#endif
#include <cstddef>
#include <cmath>
#include "def.h"
#include "buffer.h"

namespace kashiwa {
namespace runge_kutta {

  // オイラー(1768) Institutiones Calculi Integralis
  struct euler_integrator {
    static const int stage = 1;
    static const int order = 1;
    mutable working_buffer buffer;

    template<typename F>
    void operator()(double& time, double* value, std::size_t size, F const& f, double h) const {
      buffer.ensure<double>(size);
      double* ksh_restrict knode = buffer.ptr<double>();
      f(knode, time, value);
      for (std::size_t i = 0; i < size; i++)
        value[i] += h * knode[i];

      time += h;
    }
  };

  // 中点法、改良オイラー法、ルンゲ2次公式(1895)
  // - wikipedia には modified/improved のどちらも書かれていない。
  // - http://www.mymathlib.com/diffeq/runge-kutta/improved_euler_method.html
  //   によると中点法は improved Euler's method である。
  // - スペクトル法の本によるとこれは改良オイラー法である。
  struct midpoint_integrator {
    static const int stage = 2;
    static const int order = 2;
    mutable working_buffer buffer;

    template<typename F>
    void operator()(double& time, double* value, std::size_t size, F const& f, double h) const {
      buffer.ensure<double>(2 * size);
      double* ksh_restrict knode = buffer.ptr<double>();
      double* ksh_restrict xnode = buffer.ptr<double>() + size;

      f(knode, time, value);
      for (std::size_t i = 0; i < size; i++)
        xnode[i] = value[i] + 0.5 * h * knode[i];

      f(knode, time + 0.5 * h, xnode);
      for (std::size_t i = 0; i < size; i++)
        value[i] += h * knode[i];

      time += h;
    }
  };

  // ホイン法、modified Euler's method
  // - en.wikipedia には improved/modified の両方の名称が載っている。
  // - スペクトル法の本によるとこれは修正オイラー法という名前である。
  // - http://pc-physics.com/syuseieuler1.html これも修正オイラー法。
  // - http://detail.chiebukuro.yahoo.co.jp/qa/question_detail/q1091336470 このページも
  struct heun_integrator {
    static const int stage = 2;
    static const int order = 2;
    mutable working_buffer buffer;

    template<typename F>
    void operator()(double& time, double* value, std::size_t size, F const& f, double h) const {
      buffer.ensure<double>(2 * size);
      double* ksh_restrict k = buffer.ptr<double>();
      double* ksh_restrict x = buffer.ptr<double>() + size;

      f(k, time, value);
      for (std::size_t i = 0; i < size; i++) {
        x[i] = value[i] + h * k[i];
        value[i] += (1.0 / 2.0) * h * k[i];
      }

      f(k, time + h, x);
      for (std::size_t i = 0; i < size; i++)
        value[i] += (1.0 / 2.0) * h * k[i];

      time += h;
    }
  };

  struct ralston_integrator {
    static const int stage = 2;
    static const int order = 2;
    mutable working_buffer buffer;

    template<typename F>
    void operator()(double& time, double* value, std::size_t size, F const& f, double h) const {
      buffer.ensure<double>(2 * size);
      double* ksh_restrict k = buffer.ptr<double>();
      double* ksh_restrict x = buffer.ptr<double>() + size;

      f(k, time, value);
      for (std::size_t i = 0; i < size; i++) {
        x[i] = value[i] + (2.0 / 3.0) * h * k[i];
        value[i] += (1.0 / 4.0) * h * k[i];
      }

      f(k, time + (2.0 / 3.0) * h, x);
      for (std::size_t i = 0; i < size; i++)
        value[i] += (3.0 / 4.0) * h * k[i];

      time += h;
    }
  };

  //---------------------------------------------------------------------------
  // 3次公式

  /* Runge の3次公式 (4段3次の公式なので非効率)
   * - ハイラーの本
   */
  struct runge3_integrator {
    static const int stage = 4;
    static const int order = 3;
    mutable working_buffer buffer;

    template<typename F>
    void operator()(double& time, double* value, std::size_t size, F const& f, double h) const {
      buffer.ensure<double>(3 * size);
      double* ksh_restrict k = buffer.ptr<double>();
      double* ksh_restrict x = buffer.ptr<double>() + size;
      double* ksh_restrict y = buffer.ptr<double>() + size * 2;

      f(k, time, value);
      for (std::size_t i = 0; i < size; i++) {
        x[i] = value[i] + (1.0 / 2.0) * h * k[i];
        y[i] = value[i] + (1.0 / 6.0) * h * k[i];
      }

      f(k, time + (1.0 / 2.0) * h, x);
      for (std::size_t i = 0; i < size; i++) {
        x[i] = value[i] + h * k[i];
        y[i] += (2.0 / 3.0) * h * k[i];
      }

      f(k, time + h, x);
      for (std::size_t i = 0; i < size; i++)
        x[i] = value[i] + h * k[i];

      f(k, time + h, x);
      for (std::size_t i = 0; i < size; i++)
        value[i] = y[i] + (1.0 / 6.0) * h * k[i];

      time += h;
    }
  };

  /* Heun 3次法
   * - ハイラーの本
   * - http://www.330k.info/essay/Explicit-Runge-Kutta-Butcher-Tableau
   * - http://www.mymathlib.com/diffeq/runge-kutta/runge_kutta_v2_3.html
   *   には "Runge-Kutta 3次法 v2" として紹介されている。
   */
  struct heun3_integrator {
    static const int stage = 3;
    static const int order = 3;
    mutable working_buffer buffer;

    template<typename F>
    void operator()(double& time, double* value, std::size_t size, F const& f, double h) const {
      buffer.ensure<double>(2 * size);
      double* ksh_restrict k = buffer.ptr<double>();
      double* ksh_restrict x = buffer.ptr<double>() + size;

      f(k, time, value);
      for (std::size_t i = 0; i < size; i++) {
        x[i] = value[i] + (1.0 / 3.0) * h * k[i];
        value[i] += (1.0 / 4.0) * h * k[i];
      }

      f(k, time + (1.0 / 3.0) * h, x);
      for (std::size_t i = 0; i < size; i++)
        x[i] = 4.0 * value[i] - 3.0 * x[i] + (2.0 / 3.0) * h * k[i];

      f(k, time + (2.0 / 3.0) * h, x);
      for (std::size_t i = 0; i < size; i++)
        value[i] += (3.0 / 4.0) * h * k[i];

      time += h;
    }
  };

  // Ralston 3次法
  // * https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods によると、
  //   "Kutta's third-order method" である。
  struct ralston3_integrator {
    static const int stage = 3;
    static const int order = 3;
    mutable working_buffer buffer;

    template<typename F>
    void operator()(double& time, double* value, std::size_t size, F const& f, double h) const {
      buffer.ensure<double>(2 * size);
      double* ksh_restrict k = buffer.ptr<double>();
      double* ksh_restrict x = buffer.ptr<double>() + size;

      f(k, time, value);
      for (std::size_t i = 0; i < size; i++) {
        x[i] = value[i] + (1.0 / 2.0) * h * k[i];
        value[i] += (2.0 / 9.0) * h * k[i];
      }

      f(k, time + (1.0 / 2.0) * h, x);
      for (std::size_t i = 0; i < size; i++) {
        x[i] = (-4.0 / 5.0) * x[i] + (9.0 / 5.0) * value[i] + (3.0 / 4.0) * h * k[i];
        value[i] += (1.0 / 3.0) * h * k[i];
      }

      f(k, time + (3.0 / 4.0) * h, x);
      for (std::size_t i = 0; i < size; i++)
        value[i] += (4.0 / 9.0) * h * k[i];

      time += h;
    }
  };


  // Kutta 3次法、Classical RK3
  // * http://www.mymathlib.com/diffeq/runge-kutta/runge_kutta_v1_3.html
  //   によると単に "Third-order method v1" ということになっている。
  // * https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods によると、
  //   "Kutta's third-order method" である。
  // * http://www.330k.info/essay/Explicit-Runge-Kutta-Butcher-Tableau によると、
  //   Kutta 3次法、または、Classical RK3 と書かれている。
  struct kutta3_integrator {
    static const int stage = 3;
    static const int order = 3;
    mutable working_buffer buffer;

    template<typename F>
    void operator()(double& time, double* value, std::size_t size, F const& f, double h) const {
      buffer.ensure<double>(2 * size);
      double* ksh_restrict k = buffer.ptr<double>();
      double* ksh_restrict x = buffer.ptr<double>() + size;

      f(k, time, value);
      for (std::size_t i = 0; i < size; i++) {
        x[i] = value[i] + (1.0 / 2.0) * h * k[i];
        value[i] += (1.0 / 6.0) * h * k[i];
      }

      f(k, time + (1.0 / 2.0) * h, x);
      for (std::size_t i = 0; i < size; i++) {
        x[i] = (-7.0 / 2.0) * x[i] + (9.0 / 2.0) * value[i] + 2.0 * h * k[i];
        value[i] += (2.0 / 3.0) * h * k[i];
      }

      f(k, time + h, x);
      for (std::size_t i = 0; i < size; i++)
        value[i] += (1.0 / 6.0) * h * k[i];

      time += h;
    }
  };

  //---------------------------------------------------------------------------
  // 4段4次公式

  // RK4 (classical Runge-Kutta method)
  struct rk4_integrator {
    static const int stage = 4;
    static const int order = 4;
    mutable working_buffer buffer;

    template<typename F>
    void operator()(double& time, double* value, std::size_t size, F const& f, double h) const {
      buffer.ensure<double>(3 * size);
      double* ksh_restrict knode = buffer.ptr<double>();
      double* ksh_restrict xnode = buffer.ptr<double>() + size;
      double* ksh_restrict delta = buffer.ptr<double>() + size * 2;

      // k1
      f(knode, time, value);
      for (std::size_t i = 0; i < size; i++) {
        delta[i] = (1.0 / 6.0) * h * knode[i];
        xnode[i] = value[i] + 0.5 * h * knode[i];
      }

      // k2
      f(knode, time + 0.5 * h, xnode);
      for (std::size_t i = 0; i < size; i++) {
        delta[i] += (2.0 / 6.0) * h * knode[i];
        xnode[i] = value[i] + 0.5 * h * knode[i];
      }

      // k3
      f(knode, time + 0.5 * h, xnode);
      for (std::size_t i = 0; i < size; i++) {
        delta[i] += (2.0 / 6.0) * h * knode[i];
        xnode[i] = value[i] + h * knode[i];
      }

      // k4
      f(knode, time + h, xnode);
      for (std::size_t i = 0; i < size; i++) {
        double const a = delta[i] + (1.0 / 6.0) * h * knode[i];
        value[i] += a;
      }

      time += h;
    }
  };

  // Kutta 3/8-rule
  struct kutta_3_8_integrator {
    static const int stage = 4;
    static const int order = 4;
    mutable working_buffer buffer;

    template<typename F>
    void operator()(double& time, double* value, std::size_t size, F const& f, double h) const {
      buffer.ensure<double>(3 * size);
      double* ksh_restrict k  = buffer.ptr<double>();
      double* ksh_restrict xi = buffer.ptr<double>() + size;
      double* ksh_restrict x4 = buffer.ptr<double>() + size * 2;

      // k1
      f(k, time, value);
      for (std::size_t i = 0; i < size; i++) {
        xi[i] = value[i] + (1.0 / 3.0) * h * k[i];
        x4[i] = value[i] + h * k[i];
        value[i] += (1.0 / 8.0) * h * k[i];
      }

      // k2
      f(k, time + (1.0 / 3.0) * h, xi);
      for (std::size_t i = 0; i < size; i++) {
        xi[i] = 2.0 * xi[i] - x4[i] + h * k[i];
        x4[i] -= h * k[i];
        value[i] += (3.0 / 8.0) * h * k[i];
      }

      // k3
      f(k, time + (2.0 / 3.0) * h, xi);
      for (std::size_t i = 0; i < size; i++) {
        x4[i] += h * k[i];
        value[i] += (3.0 / 8.0) * h * k[i];
      }

      // k4
      f(k, time + h, x4);
      for (std::size_t i = 0; i < size; i++)
        value[i] += (1.0 / 8.0) * h * k[i];

      time += h;
    }
  };

  // Runge-Kutta Gill method
  //   Ref. gill.1 によると、丸め誤差の補正を但しく実行する為には計算手順があるとの事。
  //   [gill.1] [[Runge-Kutta-Gill法について - あらきけいすけの雑記帳>http://d.hatena.ne.jp/arakik10/20091004/p1]]
  struct gill_integrator {
    static const int stage = 4;
    static const int order = 4;
    mutable working_buffer buffer;

    template<typename F>
    void operator()(double& time, double* value, std::size_t size, F const& f, double h) const {
      buffer.ensure<double>(2 * size);
      double* ksh_restrict k  = buffer.ptr<double>();
      double* ksh_restrict q  = buffer.ptr<double>() + size;

      static constexpr double alpha1 = 1.0 / 2.0;
      static constexpr double alpha2 = 1.0 - std::sqrt(1.0 / 2.0);
      static constexpr double alpha3 = 1.0 + std::sqrt(1.0 / 2.0);
      static constexpr double alpha4 = 1.0 / 2.0;

      // k1
      f(k, time, value);
      for (std::size_t i = 0; i < size; i++) {
        double const y = value[i] + alpha1 * (h * k[i] - 2.0 * q[i]);

        // ※丸め誤差をキャンセルする為に r = ynew - yold とする必要がある。
        //   先に r を計算してから y に足すのでは駄目らしい。
        //   http://d.hatena.ne.jp/arakik10/20091004/p1
        //   http://ci.nii.ac.jp/naid/110002718589/
        double const r = y - value[i];

        q[i] = q[i] + 3.0 * r - alpha1 * h * k[i];
        value[i] = y;
      }

      // k2
      f(k, time + 0.5 * h, value);
      for (std::size_t i = 0; i < size; i++) {
        double const y = value[i] + alpha2 * (h * k[i] - q[i]);
        double const r = y - value[i];
        q[i] = q[i] + 3.0 * r - alpha2 * h * k[i];
        value[i] = y;
      }

      // k3
      f(k, time + 0.5 * h, value);
      for (std::size_t i = 0; i < size; i++) {
        double const y = value[i] + alpha3 * (h * k[i] - q[i]);
        double const r = y - value[i];
        q[i] = q[i] + 3.0 * r - alpha3 * h * k[i];
        value[i] = y;
      }

      // k4
      f(k, time + h, value);
      for (std::size_t i = 0; i < size; i++) {
        double const y = value[i] + (alpha4 / 3.0) * (h * k[i] - 2.0 * q[i]);
        double const r = y - value[i];
        q[i] = q[i] + 3.0 * r - alpha4 * h * k[i]; // ※次のステップで使う
        value[i] = y;
      }

      time += h;
    }
  };


  //---------------------------------------------------------------------------
  // 6段5次公式

  // http://www.330k.info/essay/Explicit-Runge-Kutta-Butcher-Tableau
  struct butcher5v1_integrator {
    static const int stage = 6;
    static const int order = 5;
    mutable working_buffer buffer;

    template<typename F>
    void operator()(double& time, double* ksh_restrict value, std::size_t size, F const& f, double h) const {
      buffer.ensure<double>(6 * size);
      double* ksh_restrict  x  = buffer.ptr<double>();
      double* ksh_restrict  k1 = buffer.ptr<double>() + size * 1;
      double* ksh_restrict  k2 = buffer.ptr<double>() + size * 2;
      double* ksh_restrict  k3 = buffer.ptr<double>() + size * 3;
      double* ksh_restrict  k4 = buffer.ptr<double>() + size * 4;
      double* ksh_restrict  k5 = buffer.ptr<double>() + size * 5;
      double* ksh_restrict& k6 = k2;

      // k1
      f(k1, time, value);

      // k2
      static constexpr double a21 = 1.0 / 8.0;
      static constexpr double c2  = 1.0 / 8.0;
      for (std::size_t i = 0; i < size; i++)
        x[i] = value[i] + a21 * h * k1[i];
      f(k2, time + c2 * h, x);

      // k3
      static constexpr double a32 = 1.0 / 4.0;
      static constexpr double c3  = 1.0 / 4.0;
      for (std::size_t i = 0; i < size; i++)
        x[i] = value[i] + a32 * h * k2[i];
      f(k3, time + c3 * h, x);

      // k4
      static constexpr double a41 = -1.0 / 2.0;
      static constexpr double a42 =  1.0;
      static constexpr double c4  =  1.0 / 2.0;
      for (std::size_t i = 0; i < size; i++)
        x[i] = value[i] + h * (a41 * k1[i] + a42 * k2[i]);
      f(k4, time + c4 * h, x);

      // k5
      static constexpr double a51 = 15.0 / 16.0;
      static constexpr double a52 = -3.0 /  2.0;
      static constexpr double a53 =  3.0 /  4.0;
      static constexpr double a54 =  9.0 / 16.0;
      static constexpr double c5  =  3.0 /  4.0;
      for (std::size_t i = 0; i < size; i++)
        x[i] = value[i] + h * (a51 * k1[i] + a52 * k2[i] + a53 * k3[i] + a54 * k4[i]);
      f(k5, time + c5 * h, x);

      // k6
      static constexpr double a61 = -17.0 / 7.0;
      static constexpr double a62 =   4.0;
      static constexpr double a64 = -12.0 / 7.0;
      static constexpr double a65 =   8.0 / 7.0;
      static constexpr double c6  =   1.0;
      for (std::size_t i = 0; i < size; i++)
        x[i] = value[i] + h * (a61 * k1[i] + a62 * k2[i] + a64 * k4[i] + a65 * k5[i]);
      f(k6, time + c6 * h, x);

      // increment
      static constexpr double b1  =  7.0 / 90.0;
      static constexpr double b3  = 16.0 / 45.0;
      static constexpr double b4  =  2.0 / 15.0;
      static constexpr double b5  = 16.0 / 45.0;
      static constexpr double b6  =  7.0 / 90.0;
      for (std::size_t i = 0; i < size; i++)
        value[i] += h * (b1 * k1[i] + b3 * k3[i] + b4 * k4[i] + b5 * k5[i] + b6 * k6[i]);

      time += h;
    }
  };

  // http://www.330k.info/essay/Explicit-Runge-Kutta-Butcher-Tableau
  struct butcher5v2_integrator {
    static const int stage = 6;
    static const int order = 5;
    mutable working_buffer buffer;

    template<typename F>
    void operator()(double& time, double* ksh_restrict value, std::size_t size, F const& f, double h) const {
      buffer.ensure<double>(5 * size);
      double* ksh_restrict  x  = buffer.ptr<double>();
      double* ksh_restrict  y  = buffer.ptr<double>() + size * 1;
      double* ksh_restrict  k1 = buffer.ptr<double>() + size * 2;
      double* ksh_restrict  k2 = buffer.ptr<double>() + size * 3;
      double* ksh_restrict& k3 = k2;
      double* ksh_restrict& k4 = k2;
      double* ksh_restrict& k5 = k2;
      double* ksh_restrict  k6 = buffer.ptr<double>() + size * 4;

      static constexpr double a21 = 1.0 / 4.0;
      static constexpr double c2  = 1.0 / 4.0;

      static constexpr double a31 = 1.0 / 8.0;
      static constexpr double a32 = 1.0 / 8.0;
      static constexpr double c3  = 1.0 / 4.0;

      static constexpr double a42 = -1.0 / 2.0;
      static constexpr double a43 =  1.0;
      static constexpr double c4  =  1.0 / 2.0;

      static constexpr double a51 = 3.0 / 16.0;
      static constexpr double a54 = 9.0 / 16.0;
      static constexpr double c5  = 3.0 /  4.0;

      static constexpr double a61 =  -3.0 / 7.0;
      static constexpr double a62 =   2.0 / 7.0;
      static constexpr double a63 =  12.0 / 7.0;
      static constexpr double a64 = -12.0 / 7.0;
      static constexpr double a65 =   8.0 / 7.0;
      static constexpr double c6  =   1.0;

      static constexpr double b1  =  7.0 / 90.0;
      static constexpr double b3  = 16.0 / 45.0;
      static constexpr double b4  =  2.0 / 15.0;
      static constexpr double b5  = 16.0 / 45.0;
      static constexpr double b6  =  7.0 / 90.0;

      // k1
      f(k1, time, value);

      // k2
      for (std::size_t i = 0; i < size; i++) {
        x[i]  = value[i] + a21 * h * k1[i];
        k6[i] = value[i] + a61 * h * k1[i];
        y[i]  = value[i] + b1 * h * k1[i];
      }
      f(k2, time + c2 * h, x);

      // k3
      for (std::size_t i = 0; i < size; i++) {
        x[i] = value[i] + h * (a31 * k1[i] + a32 * k2[i]);
        k6[i] += a62 * h * k2[i];
      }
      f(k3, time + c3 * h, x);

      // k4
      for (std::size_t i = 0; i < size; i++) {
        x[i] = a42 / a32 * x[i] + (1.0 - a42 / a32) * value[i] - ((a42 / a32) * a31) * h * k1[i] + a43 * h * k3[i];
        k6[i] += a63 * h * k3[i];
        y[i] += b3 * h * k3[i];
      }
      f(k4, time + c4 * h, x);

      // k5
      for (std::size_t i = 0; i < size; i++) {
        x[i] = value[i] + h * (a51 * k1[i] + a54 * k4[i]);
        k6[i] += a64 * h * k4[i];
        y[i] += b4 * h * k4[i];
      }
      f(k5, time + c5 * h, x);

      // k6
      for (std::size_t i = 0; i < size; i++) {
        x[i] = k6[i] + a65 * h * k5[i];
        y[i] += b5 * h * k5[i];
      }
      f(k6, time + c6 * h, x);

      // increment
      for (std::size_t i = 0; i < size; i++)
        value[i] = y[i] + b6 * h * k6[i];

      time += h;

    }
  };

  // http://www.330k.info/essay/Explicit-Runge-Kutta-Butcher-Tableau
  struct butcher5v3_integrator {
    static const int stage = 6;
    static const int order = 5;
    mutable working_buffer buffer;

    template<typename F>
    void operator()(double& time, double* ksh_restrict value, std::size_t size, F const& f, double h) const {
      buffer.ensure<double>(6 * size);
      double* ksh_restrict  x  = buffer.ptr<double>();
      double* ksh_restrict  k1 = buffer.ptr<double>() + size * 1;
      double* ksh_restrict  k2 = buffer.ptr<double>() + size * 2;
      double* ksh_restrict  k3 = buffer.ptr<double>() + size * 3;
      double* ksh_restrict  k4 = buffer.ptr<double>() + size * 4;
      double* ksh_restrict  k5 = buffer.ptr<double>() + size * 5;
      double* ksh_restrict& k6 = k2;

      // k1
      f(k1, time, value);

      // k2
      static constexpr double a21 = -1.0 / 2.0;
      static constexpr double c2  = -1.0 / 2.0;
      for (std::size_t i = 0; i < size; i++)
        x[i] = value[i] + a21 * h * k1[i];
      f(k2, time + c2 * h, x);

      // k3
      static constexpr double a31 =  5.0 / 16.0;
      static constexpr double a32 = -1.0 / 16.0;
      static constexpr double c3  =  1.0 /  4.0;
      for (std::size_t i = 0; i < size; i++)
        x[i] = value[i] + h * (a31 * k1[i] + a32 * k2[i]);
      f(k3, time + c3 * h, x);

      // k4
      static constexpr double a41 = -3.0 / 4.0;
      static constexpr double a42 =  1.0 / 4.0;
      static constexpr double a43 =  1.0;
      static constexpr double c4  =  1.0 / 2.0;
      for (std::size_t i = 0; i < size; i++)
        x[i] = value[i] + h * (a41 * k1[i] + a42 * k2[i] + a43 * k3[i]);
      f(k4, time + c4 * h, x);

      // k5
      static constexpr double a51 = 3.0 / 16.0;
      static constexpr double a54 = 9.0 / 16.0;
      static constexpr double c5  = 3.0 /  4.0;
      for (std::size_t i = 0; i < size; i++)
        x[i] = value[i] + h * (a51 * k1[i] + a54 * k4[i]);
      f(k5, time + c5 * h, x);

      // k6
      static constexpr double a62 =  -1.0 / 7.0;
      static constexpr double a63 =  12.0 / 7.0;
      static constexpr double a64 = -12.0 / 7.0;
      static constexpr double a65 =   8.0 / 7.0;
      static constexpr double c6  = 1.0;
      for (std::size_t i = 0; i < size; i++)
        x[i] = value[i] + h * (a62 * k2[i] + a63 * k3[i] + a64 * k4[i] + a65 * k5[i]);
      f(k6, time + c6 * h, x);

      // increment
      static constexpr double b1  =  7.0 / 90.0;
      static constexpr double b3  = 16.0 / 45.0;
      static constexpr double b4  =  2.0 / 15.0;
      static constexpr double b5  = 16.0 / 45.0;
      static constexpr double b6  =  7.0 / 90.0;
      for (std::size_t i = 0; i < size; i++)
        value[i] += h * (b1 * k1[i] + b3 * k3[i] + b4 * k4[i] + b5 * k5[i] + b6 * k6[i]);

      time += h;
    }
  };

  //---------------------------------------------------------------------------
  // 高次公式

  // http://www.330k.info/essay/Explicit-Runge-Kutta-Butcher-Tableau
  struct hammud6_integrator {
    static const int stage = 7;
    static const int order = 6;
    mutable working_buffer buffer;

    static constexpr double sqrt21 = std::sqrt(21.0); // Ref [cv8.2] では -sqrt(21.0). どちらでも OK.

    template<typename F>
    void operator()(double& time, double* ksh_restrict value, std::size_t size, F const& f, double h) const {
      buffer.ensure<double>(7 * size);
      double* ksh_restrict  x  = buffer.ptr<double>();
      double* ksh_restrict  k1 = buffer.ptr<double>() + size * 1;
      double* ksh_restrict  k2 = buffer.ptr<double>() + size * 2;
      double* ksh_restrict  k3 = buffer.ptr<double>() + size * 3;
      double* ksh_restrict  k4 = buffer.ptr<double>() + size * 4;
      double* ksh_restrict  k5 = buffer.ptr<double>() + size * 5;
      double* ksh_restrict  k6 = buffer.ptr<double>() + size * 6;
      double* ksh_restrict& k7 = k2;

      static constexpr double sqrt5 = std::sqrt(5);
      static constexpr double b1 = 1.0 / 12.0;
      static constexpr double b5 = 5.0 / 12.0;
      static constexpr double b6 = 5.0 / 12.0;
      static constexpr double b7 = 1.0 / 12.0;


      // k1
      f(k1, time, value);

      // k2
      static constexpr double a21 = 4.0 / 7.0;
      static constexpr double c2  = 4.0 / 7.0;
      for (std::size_t i = 0; i < size; i++)
        x[i] = value[i] + a21 * h * k1[i];
      f(k2, time + c2 * h, x);

      // k3
      static constexpr double a31 = 115.0 / 112.0;
      static constexpr double a32 =  -5.0 /  16.0;
      static constexpr double c3  =   5.0 /   7.0;
      for (std::size_t i = 0; i < size; i++)
        x[i] = value[i] + h * (a31 * k1[i] + a32 * k2[i]);
      f(k3, time + c3 * h, x);

      // k4
      static constexpr double a41 = 589.0 / 630.0;
      static constexpr double a42 =   5.0 /  18.0;
      static constexpr double a43 = -16.0 /  45.0;
      static constexpr double c4  =   6.0 /   7.0;
      for (std::size_t i = 0; i < size; i++)
        x[i] = value[i] + h * (a41 * k1[i] + a42 * k2[i] + a43 * k3[i]);
      f(k4, time + c4 * h, x);

      // k5
      static constexpr double a51 = 229.0 / 1200.0 -  29.0 * sqrt5 / 6000.0;
      static constexpr double a52 = 119.0 /  240.0 - 187.0 * sqrt5 / 1200.0;
      static constexpr double a53 = -14.0 /   75.0 +  34.0 * sqrt5 /  375.0;
      static constexpr double a54 =                   -3.0 * sqrt5 /  100.0;
      static constexpr double c5  = (5.0 - sqrt5) / 10.0;
      for (std::size_t i = 0; i < size; i++)
        x[i] = value[i] + h * (a51 * k1[i] + a52 * k2[i] + a53 * k3[i] + a54 * k4[i]);
      f(k5, time + c5 * h, x);

      // k6
      static constexpr double a61 =  71.0 / 2400.0 - 587.0 * sqrt5 / 12000.0;
      static constexpr double a62 = 187.0 / 480.0  - 391.0 * sqrt5 /  2400.0;
      static constexpr double a63 = -38.0 / 75.0   +  26.0 * sqrt5 /   375.0;
      static constexpr double a64 =  27.0 / 80.0   -   3.0 * sqrt5 /   400.0;
      static constexpr double a65 =   1.0 / 4.0    +         sqrt5 /     4.0;
      static constexpr double c6  = (5.0 + sqrt5) / 10.0;
      for (std::size_t i = 0; i < size; i++)
        x[i] = value[i] + h * (a61 * k1[i] + a62 * k2[i] + a63 * k3[i] + a64 * k4[i] + a65 * k5[i]);
      f(k6, time + c6 * h, x);

      // k7 <= k2
      static constexpr double a71 = -49.0 / 480.0 + 43.0 * sqrt5 / 160.0;
      static constexpr double a72 = -425.0 / 96.0 + 51.0 * sqrt5 /  32.0;
      static constexpr double a73 =   52.0 / 15.0 -  4.0 * sqrt5 /   5.0;
      static constexpr double a74 =  -27.0 / 16.0 +  3.0 * sqrt5 /  16.0;
      static constexpr double a75 =    5.0 / 4.0  -  3.0 * sqrt5 /   4.0;
      static constexpr double a76 =    5.0 / 2.0  -        sqrt5 /   2.0;
      static constexpr double c7 = 1.0;
      for (std::size_t i = 0; i < size; i++)
        x[i] = value[i] + h * (a71 * k1[i] + a72 * k2[i] + a73 * k3[i] + a74 * k4[i] + a75 * k5[i] + a76 * k6[i]);
      f(k7, time + c7 * h, x);

      // increment
      for (std::size_t i = 0; i < size; i++)
        value[i] += h * (b1 * k1[i] + b5 * k5[i] + b6 * k6[i] + b7 * k7[i]);

      time += h;
    }
  };

  // http://www.330k.info/essay/Explicit-Runge-Kutta-Butcher-Tableau
  struct shanks7_integrator {
    static const int stage = 9;
    static const int order = 7;
    mutable working_buffer buffer;

    static constexpr double sqrt21 = std::sqrt(21.0); // Ref [cv8.2] では -sqrt(21.0). どちらでも OK.

    template<typename F>
    void operator()(double& time, double* ksh_restrict value, std::size_t size, F const& f, double h) const {
      buffer.ensure<double>(9 * size);
      double* ksh_restrict  x  = buffer.ptr<double>();
      double* ksh_restrict  k1 = buffer.ptr<double>() + size * 1;
      double* ksh_restrict  k2 = buffer.ptr<double>() + size * 2;
      double* ksh_restrict& k3 = k2;
      double* ksh_restrict  k4 = buffer.ptr<double>() + size * 4;
      double* ksh_restrict  k5 = buffer.ptr<double>() + size * 5;
      double* ksh_restrict  k6 = buffer.ptr<double>() + size * 6;
      double* ksh_restrict  k7 = buffer.ptr<double>() + size * 7;
      double* ksh_restrict  k8 = buffer.ptr<double>() + size * 8;
      double* ksh_restrict& k9 = k2;

      // k1
      f(k1, time, value);

      // k2
      static constexpr double a21 = 2.0 / 9.0;
      static constexpr double c2  = 2.0 / 9.0;
      for (std::size_t i = 0; i < size; i++)
        x[i] = value[i] + a21 * h * k1[i];
      f(k2, time + c2 * h, x);

      // k3
      static constexpr double a31 = 1.0 / 12.0;
      static constexpr double a32 = 1.0 / 4.0;
      static constexpr double c3  = 1.0 / 3.0;
      for (std::size_t i = 0; i < size; i++)
        x[i] = value[i] + h * (a31 * k1[i] + a32 * k2[i]);
      f(k3, time + c3 * h, x);

      // k4
      static constexpr double a41 = 1.0 / 8.0;
      static constexpr double a43 = 3.0 / 8.0;
      static constexpr double c4  = 1.0 / 2.0;
      for (std::size_t i = 0; i < size; i++)
        x[i] = value[i] + h * (a41 * k1[i] + a43 * k3[i]);
      f(k4, time + c4 * h, x);

      // k5
      static constexpr double a51 = 23.0 / 216.0;
      static constexpr double a53 =  7.0 / 72.0;
      static constexpr double a54 = -1.0 / 27.0;
      static constexpr double c5  =  1.0 / 6.0;
      for (std::size_t i = 0; i < size; i++)
        x[i] = value[i] + h * (a51 * k1[i] + a53 * k3[i] + a54 * k4[i]);
      f(k5, time + c5 * h, x);

      // k6
      static constexpr double a61 = -4136.0 / 729.0;
      static constexpr double a63 = -4528.0 / 243.0;
      static constexpr double a64 =  5264.0 / 729.0;
      static constexpr double a65 =  1456.0 / 81.0;
      static constexpr double c6  =     8.0 / 9.0;
      for (std::size_t i = 0; i < size; i++)
        x[i] = value[i] + h * (a61 * k1[i] + a63 * k3[i] + a64 * k4[i] + a65 * k5[i]);
      f(k6, time + c6 * h, x);

      // k7
      static constexpr double a71 = 8087.0 / 11664.0;
      static constexpr double a73 =  484.0 / 243.0;
      static constexpr double a74 = -518.0 / 729.0;
      static constexpr double a75 = -658.0 / 351.0;
      static constexpr double a76 =    7.0 / 624.0;
      static constexpr double c7  =    1.0 / 9.0;
      for (std::size_t i = 0; i < size; i++)
        x[i] = value[i] + h * (a71 * k1[i] + a73 * k3[i] + a74 * k4[i] + a75 * k5[i] + a76 * k6[i]);
      f(k7, time + c7 * h, x);

      // k8
      static constexpr double a81 = -1217.0 / 2160.0;
      static constexpr double a83 =  -145.0 / 72.0;
      static constexpr double a84 =  8342.0 / 6615.0;
      static constexpr double a85 =   361.0 / 195.0;
      static constexpr double a86 =  3033.0 / 50960.0;
      static constexpr double a87 =   117.0 / 490.0;
      static constexpr double c8  =     5.0 / 6.0;
      for (std::size_t i = 0; i < size; i++)
        x[i] = value[i] + h * (a81 * k1[i] + a83 * k3[i] + a84 * k4[i] + a85 * k5[i] + a86 * k6[i] + a87 * k7[i]);
      f(k8, time + c8 * h, x);

      // k9
      static constexpr double a91 =    259.0 / 2768.0;
      static constexpr double a93 =    -84.0 / 173.0;
      static constexpr double a94 =    -14.0 / 173.0;
      static constexpr double a95 =   6210.0 / 2249.0;
      static constexpr double a96 = -99873.0 / 251888.0;
      static constexpr double a97 = -29160.0 / 15743.0;
      static constexpr double a98 =   2160.0 / 2249.0;
      static constexpr double c9  = 1.0;
      for (std::size_t i = 0; i < size; i++)
        x[i] = value[i] + h * (a91 * k1[i] + a93 * k3[i] + a94 * k4[i] + a95 * k5[i] + a96 * k6[i] + a97 * k7[i] + a98 * k8[i]);
      f(k9, time + c9 * h, x);

      // increment
      static constexpr double b1 =    173.0 / 3360.0;
      static constexpr double b4 =   1846.0 / 5145.0;
      static constexpr double b5 =     27.0 / 91.0;
      static constexpr double b6 = -19683.0 / 713440.0;
      static constexpr double b7 = -19683.0 / 713440.0;
      static constexpr double b8 =     27.0 / 91.0;
      static constexpr double b9 =    173.0 / 3360.0;
      for (std::size_t i = 0; i < size; i++)
        value[i] += h * (b1 * k1[i] + b4 * k4[i] + b5 * k5[i] + b6 * k6[i] + b7 * k7[i] + b8 * k8[i] + b9 * k9[i]);

      time += h;
    }
  };

  // http://www.330k.info/essay/Explicit-Runge-Kutta-Butcher-Tableau
  struct cooper_verner7_integrator {
    static const int stage = 9;
    static const int order = 7;
    mutable working_buffer buffer;

    static constexpr double sqrt21 = std::sqrt(21.0);

    template<typename F>
    void operator()(double& time, double* ksh_restrict value, std::size_t size, F const& f, double h) const {
      buffer.ensure<double>(9 * size);
      double* ksh_restrict  x  = buffer.ptr<double>();
      double* ksh_restrict  k1 = buffer.ptr<double>() + size * 1;
      double* ksh_restrict  k2 = buffer.ptr<double>() + size * 2;
      double* ksh_restrict& k3 = k2;
      double* ksh_restrict  k4 = buffer.ptr<double>() + size * 4;
      double* ksh_restrict  k5 = buffer.ptr<double>() + size * 5;
      double* ksh_restrict  k6 = buffer.ptr<double>() + size * 6;
      double* ksh_restrict  k7 = buffer.ptr<double>() + size * 7;
      double* ksh_restrict  k8 = buffer.ptr<double>() + size * 8;
      double* ksh_restrict& k9 = k2;

      // k1
      f(k1, time, value);

      // k2
      static constexpr double a21 = (7.0 + 1.0 * sqrt21) / 42.0;
      static constexpr double c2  = (7.0 +       sqrt21) / 42.0;
      for (std::size_t i = 0; i < size; i++)
        x[i] = value[i] + a21 * h * k1[i];
      f(k2, time + c2 * h, x);

      // k3
      static constexpr double a32 = (7.0 + 1.0 * sqrt21) / 21.0;
      static constexpr double c3  = (7.0 +       sqrt21) / 21.0;
      for (std::size_t i = 0; i < size; i++)
        x[i] = value[i] + a32 * h * k2[i];
      f(k3, time + c3 * h, x);

      // k4
      static constexpr double a41 = ( 7.0 + 1.0 * sqrt21) / 56.0;
      static constexpr double a43 = (21.0 + 3.0 * sqrt21) / 56.0;
      static constexpr double c4  = ( 7.0 +       sqrt21) / 14.0;
      for (std::size_t i = 0; i < size; i++)
        x[i] = value[i] + h * (a41 * k1[i] + a43 * k3[i]);
      f(k4, time + c4 * h, x);

      // k5
      static constexpr double a51 = (  8.0 - 1.0 * sqrt21) / 16.0;
      static constexpr double a53 = (-21.0 + 6.0 * sqrt21) / 16.0;
      static constexpr double a54 = ( 21.0 - 5.0 * sqrt21) / 16.0;
      static constexpr double c5 = 1.0 / 2.0;
      for (std::size_t i = 0; i < size; i++)
        x[i] = value[i] + h * (a51 * k1[i] + a53 * k3[i] + a54 * k4[i]);
      f(k5, time + c5 * h, x);

      // k6
      static constexpr double a61 = (-1687.0 + 374.0 * sqrt21) / 196.0;
      static constexpr double a63 = (  969.0 - 210.0 * sqrt21) /  28.0;
      static constexpr double a64 = ( -381.0 +  83.0 * sqrt21) /  14.0;
      static constexpr double a65 = (   84.0 -  20.0 * sqrt21) /  49.0;
      static constexpr double c6  = (    7.0 -         sqrt21) /  14.0;
      for (std::size_t i = 0; i < size; i++)
        x[i] = value[i] + h * (a61 * k1[i] + a63 * k3[i] + a64 * k4[i] + a65 * k5[i]);
      f(k6, time + c6 * h, x);

      // k7
      static constexpr double a71 = (  583.0 - 131.0 * sqrt21) / 128.0;
      static constexpr double a73 = (-2373.0 + 501.0 * sqrt21) / 128.0;
      static constexpr double a74 = ( 4221.0 - 914.0 * sqrt21) / 288.0;
      static constexpr double a75 = (   -9.0 +   4.0 * sqrt21) /  18.0;
      static constexpr double a76 = (  189.0 +  35.0 * sqrt21) / 576.0;
      static constexpr double c7 = 1.0 / 2.0;
      for (std::size_t i = 0; i < size; i++)
        x[i] = value[i] + h * (a71 * k1[i] + a73 * k3[i] + a74 * k4[i] + a75 * k5[i] + a76 * k6[i]);
      f(k7, time + c7 * h, x);

      // k8
      static constexpr double a81 = ( -623.0 +  169.0 * sqrt21) /  392.0;
      static constexpr double a83 = (  435.0 -   81.0 * sqrt21) /   56.0;
      static constexpr double a84 = (-1437.0 +  307.0 * sqrt21) /  252.0;
      static constexpr double a85 = (-2028.0 - 1468.0 * sqrt21) / 7497.0;
      static constexpr double a86 = (  -21.0 -    4.0 * sqrt21) /  126.0;
      static constexpr double a87 = (  384.0 +   80.0 * sqrt21) /  833.0;
      static constexpr double c8  = (    7.0 +          sqrt21) /   14.0;
      for (std::size_t i = 0; i < size; i++)
        x[i] = value[i] + h * (a81 * k1[i] + a83 * k3[i] + a84 * k4[i] + a85 * k5[i] + a86 * k6[i] + a87 * k7[i]);
      f(k8, time + c8 * h, x);

      // k9
      static constexpr double a91 = (  579.0 -  131.0 * sqrt21) /  24.0;
      static constexpr double a93 = ( -791.0 +  167.0 * sqrt21) /   8.0;
      static constexpr double a94 = ( 8099.0 - 1765.0 * sqrt21) / 108.0;
      static constexpr double a95 = (-1976.0 +  784.0 * sqrt21) / 459.0;
      static constexpr double a96 = (   70.0 +    7.0 * sqrt21) /  54.0;
      static constexpr double a97 = (  160.0 -   80.0 * sqrt21) / 153.0;
      static constexpr double a98 = (   49.0 -    7.0 * sqrt21) /  18.0;
      static constexpr double c9 = 1.0;
      for (std::size_t i = 0; i < size; i++)
        x[i] = value[i] + h * (a91 * k1[i] + a93 * k3[i] + a94 * k4[i] + a95 * k5[i] + a96 * k6[i] + a97 * k7[i] + a98 * k8[i]);
      f(k9, time + c9 * h, x);

      // increment
      static constexpr double b1 =  1.0 /  20.0;
      static constexpr double b6 = 49.0 / 180.0;
      static constexpr double b7 = 16.0 /  45.0;
      static constexpr double b8 = 49.0 / 180.0;
      static constexpr double b9 =  1.0 /  20.0;
      for (std::size_t i = 0; i < size; i++)
        value[i] += h * (b1 * k1[i] + b6 * k6[i] + b7 * k7[i] + b8 * k8[i] + b9 * k9[i]);

      time += h;
    }
  };

  // Cooper-Verner 8th order
  //   Ref. cv8.3 によると "Verner's 8th order method" という名前である。
  //   Ref. cv8.2 に載っているのは sqrt(21) の前の符号が反転している。どちらでも良いようだ。
  //   [cv8.1] E.ハイラー, 三井 斌友, 『常微分方程式の数値解法 I 』, 丸善出版 (2012/7/17).
  //   [cv8.2] http://www.330k.info/essay/Explicit-Runge-Kutta-Butcher-Tableau
  //   [cv8.3] http://www.mymathlib.com/diffeq/runge-kutta/runge_kutta_verner.html
  struct cooper_verner8_integrator {
    static const int stage = 11;
    static const int order = 8;
    mutable working_buffer buffer;

    static constexpr double sqrt21 = std::sqrt(21.0); // Ref [cv8.2] では -sqrt(21.0). どちらでも OK.

    template<typename F>
    void operator()(double& time, double* ksh_restrict value, std::size_t size, F const& f, double h) const {
      buffer.ensure<double>(8 * size);
      double* ksh_restrict  delta = buffer.ptr<double>();
      double* ksh_restrict  xnode = buffer.ptr<double>() + size;
      double* ksh_restrict  k1    = buffer.ptr<double>() + size * 2;
      double* ksh_restrict  k2    = buffer.ptr<double>() + size * 3;
      double* ksh_restrict  k3    = buffer.ptr<double>() + size * 4;
      double* ksh_restrict& k4    = k2;
      double* ksh_restrict  k5    = buffer.ptr<double>() + size * 5;
      double* ksh_restrict  k6    = buffer.ptr<double>() + size * 6;
      double* ksh_restrict& k7    = k3;
      double* ksh_restrict& k8    = k4;
      double* ksh_restrict  k9    = buffer.ptr<double>() + size * 7;
      double* ksh_restrict& kA    = k1;
      double* ksh_restrict& kB    = k5;

      static constexpr double b1 =  1.0 /  20.0;
      static constexpr double b8 = 49.0 /  180.0;
      static constexpr double b9 = 16.0 /  45.0;
      static constexpr double bA = 49.0 / 180.0;
      static constexpr double bB =  1.0 /  20.0;

      // k1
      f(k1, time, value);

      // k2
      static constexpr double a21 = 0.5;
      static constexpr double c20 = a21;
      for (std::size_t i = 0; i < size; i++) {
        delta[i] = b1 * h * k1[i];
        xnode[i] = value[i] + a21 * h * k1[i];
      }
      f(k2, time + c20 * h, xnode);

      // k3
      static constexpr double a31 = 0.25;
      static constexpr double a32 = 0.25;
      static constexpr double c30 = a31 + a32;
      for (std::size_t i = 0; i < size; i++)
        xnode[i] = value[i] + a31 * h * k1[i] + a32 * h * k2[i];
      f(k3, time + c30 * h, xnode);

      // k4 <= k2
      static constexpr double a41 = (1.0 /  7.0);
      static constexpr double a42 = (1.0 / 98.0) * (-7 - 3 * sqrt21);
      static constexpr double a43 = (1.0 / 49.0) * (21 + 5 * sqrt21);
      static constexpr double c40 = a41 + a42 + a43;
      for (std::size_t i = 0; i < size; i++)
        xnode[i] = value[i] + a41 * h * k1[i] + a42 * h * k2[i] + a43 * h * k3[i];
      f(k4, time + c40 * h, xnode);

      // k5
      static constexpr double a51 = (1.0 /  84.0) * (11 + 1 * sqrt21);
      static constexpr double a53 = (1.0 /  63.0) * (18 + 4 * sqrt21);
      static constexpr double a54 = (1.0 / 252.0) * (21 - 1 * sqrt21);
      static constexpr double c50 = a51 + a53 + a54;
      for (std::size_t i = 0; i < size; i++)
        xnode[i] = value[i] + a51 * h * k1[i] + a53 * h * k3[i] + a54 * h * k4[i];
      f(k5, time + c50 * h, xnode);

      // k6
      static constexpr double a61 = (1.0 /  48.0) * (   5 +  1 * sqrt21);
      static constexpr double a63 = (1.0 /  36.0) * (   9 +  1 * sqrt21);
      static constexpr double a64 = (1.0 / 360.0) * (-231 + 14 * sqrt21);
      static constexpr double a65 = (1.0 /  80.0) * (  63 -  7 * sqrt21);
      static constexpr double c60 = a61 + a63 + a64 + a65;
      for (std::size_t i = 0; i < size; i++)
        xnode[i] = value[i] + a61 * h * k1[i] + a63 * h * k3[i] + a64 * h * k4[i] + a65 * h * k5[i];
      f(k6, time + c60 * h, xnode);

      // k7 <= k3
      static constexpr double a71 = (1.0 /  42.0) * (  10 -   1 * sqrt21);
      static constexpr double a73 = (1.0 / 315.0) * (-432 +  92 * sqrt21);
      static constexpr double a74 = (1.0 /  90.0) * ( 633 - 145 * sqrt21);
      static constexpr double a75 = (1.0 /  70.0) * (-504 + 115 * sqrt21);
      static constexpr double a76 = (1.0 /  35.0) * (  63 -  13 * sqrt21);
      static constexpr double c70 = a71 + a73 + a74 + a75 + a76;
      for (std::size_t i = 0; i < size; i++)
        xnode[i] = value[i] + a71 * h * k1[i] + a73 * h * k3[i] + a74 * h * k4[i] + a75 * h * k5[i] + a76 * h * k6[i];
      f(k7, time + c70 * h, xnode);

      // k8 <= k4
      static constexpr double a81 =  1.0 /  14.0;
      static constexpr double a85 = (1.0 / 126.0) * (14 - 3 * sqrt21);
      static constexpr double a86 = (1.0 /  63.0) * (13 - 3 * sqrt21);
      static constexpr double a87 =  1.0 /   9.0;
      static constexpr double c80 = a81 + a85 + a86 + a87;
      for (std::size_t i = 0; i < size; i++)
        xnode[i] = value[i] + a81 * h * k1[i] + a85 * h * k5[i] + a86 * h * k6[i] + a87 * h * k7[i];
      f(k8, time + c80 * h, xnode);

      // k9
      static constexpr double a91 =   1.0 /   32.0;
      static constexpr double a95 = ( 1.0 /  576.0) * (  91.0 - 21.0 * sqrt21);
      static constexpr double a96 =  11.0 /   72.0;
      static constexpr double a97 = ( 1.0 / 1152.0) * (-385.0 - 75.0 * sqrt21);
      static constexpr double a98 = ( 1.0 /  128.0) * (  63.0 + 13.0 * sqrt21);
      static constexpr double c90 = a91 + a95 + a96 + a97 + a98;
      for (std::size_t i = 0; i < size; i++) {
        delta[i] += b8 * h * k8[i];
        xnode[i] = value[i] + a91 * h * k1[i] + a95 * h * k5[i] + a96 * h * k6[i] + a97 * h * k7[i] + a98 * h * k8[i];
      }
      f(k9, time + c90 * h, xnode);

      // kA <= k1
      static constexpr double aA1 =  1.0 /   14.0;
      static constexpr double aA5 =  1.0 /    9.0;
      static constexpr double aA6 = (1.0 / 2205.0) * (-733.0 - 147.0 * sqrt21);
      static constexpr double aA7 = (1.0 /  504.0) * ( 515.0 + 111.0 * sqrt21);
      static constexpr double aA8 = (1.0 /   56.0) * (- 51.0 -  11.0 * sqrt21);
      static constexpr double aA9 = (1.0 /  245.0) * ( 132.0 +  28.0 * sqrt21);
      static constexpr double cA0 = aA1 + aA5 + aA6 + aA7 + aA8 + aA9;
      for (std::size_t i = 0; i < size; i++) {
        delta[i] += b9 * h * k9[i];
        xnode[i] = value[i] + aA1 * h * k1[i] + aA5 * h * k5[i] + aA6 * h * k6[i] + aA7 * h * k7[i] + aA8 * h * k8[i] + aA9 * h * k9[i];
      }
      f(kA, time + cA0 * h, xnode);

      // kB <= k5
      static constexpr double aB5 = (1.0 / 18.0) * (- 42.0 +  7.0 * sqrt21);
      static constexpr double aB6 = (1.0 / 45.0) * (- 18.0 + 28.0 * sqrt21);
      static constexpr double aB7 = (1.0 / 72.0) * (-273.0 - 53.0 * sqrt21);
      static constexpr double aB8 = (1.0 / 72.0) * ( 301.0 + 53.0 * sqrt21);
      static constexpr double aB9 = (1.0 / 45.0) * (  28.0 - 28.0 * sqrt21);
      static constexpr double aBA = (1.0 / 18.0) * (  49.0 -  7.0 * sqrt21);
      static constexpr double cB0 = aB5 + aB6 + aB7 + aB8 + aB9 + aBA;
      for (std::size_t i = 0; i < size; i++) {
        delta[i] += bA * h * kA[i];
        xnode[i] = value[i] + aB5 * h * k5[i] + aB6 * h * k6[i] + aB7 * h * k7[i] + aB8 * h * k8[i] + aB9 * h * k9[i] + aBA * h * kA[i];
      }
      f(kB, time + cB0 * h, xnode);

      // increment
      for (std::size_t i = 0; i < size; i++) {
        double const a = delta[i] + bB * h * kB[i];
        value[i] += a;
      }

      time += h;
    }
  };

  // 他
  // * http://www.mymathlib.com/diffeq/runge-kutta/runge_kutta_ralston_4.html Ralston's 4th
  // * ハイラーの Ralson(1962), Hull(1967) (u, v) = (0.4, 0.45) と同じ物か?
  // * http://www.mymathlib.com/diffeq/runge-kutta/runge_kutta_nystrom.html Nystrom's 5th
  // * http://www.mymathlib.com/diffeq/runge-kutta/runge_kutta_butcher.html Butcher's 6th

}
}

#endif
