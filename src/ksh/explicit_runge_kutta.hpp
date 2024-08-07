// -*- mode:c++ -*-
#ifndef kashiwa_explicit_runge_kutta_hpp
#define kashiwa_explicit_runge_kutta_hpp
#ifdef _MSC_VER
# define _USE_MATH_DEFINES
#endif
#include <cstddef>
#include <cmath>
#include "def.hpp"
#include "buffer.hpp"
#include "ode.hpp"

// References
//
// [1] [陽的Runge-Kutta法のButcher tableau - 330k info](http://www.330k.info/essay/Explicit-Runge-Kutta-Butcher-Tableau)
//   以下に移動した様だ: https://www.330k.info/essay/explicit-runge-kutta-butcher-tableau/
//   Web Archive にも残っている https://web.archive.org/web/20220121115842/https://www.330k.info/essay/explicit-runge-kutta-butcher-tableau/
//

namespace kashiwa {
namespace runge_kutta {

  // オイラー(1768) Institutiones Calculi Integralis
  struct euler_integrator: explicit_integrator_base<euler_integrator> {
    static const int stage = 1;
    static const int order = 1;
    mutable working_buffer buffer;

    template<typename F>
    void operator()(double& time, double* value, std::size_t size, F const& f, double h) const {
      buffer.ensure<double>(size);
      double* ksh_restrict knode = buffer.ptr<double>();
      f(knode, size, time, value);
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
  struct midpoint_integrator: explicit_integrator_base<midpoint_integrator> {
    static const int stage = 2;
    static const int order = 2;
    mutable working_buffer buffer;

    template<typename F>
    void operator()(double& time, double* value, std::size_t size, F const& f, double h) const {
      buffer.ensure<double>(2 * size);
      double* ksh_restrict knode = buffer.ptr<double>();
      double* ksh_restrict xnode = buffer.ptr<double>() + size;

      f(knode, size, time, value);
      for (std::size_t i = 0; i < size; i++)
        xnode[i] = value[i] + 0.5 * h * knode[i];

      f(knode, size, time + 0.5 * h, xnode);
      for (std::size_t i = 0; i < size; i++)
        value[i] += h * knode[i];

      time += h;
    }
  };

  // ホイン法、または modified Euler's method。文献 [1] によると
  // improved/modified の両方の名前が挙げられている。文献[2,3,4]ではこれを修正
  // オイラー法と呼んでいる。TVD-RK3 (SSP-RK2) はこれに等価である[5,6]。
  // - [1] en.wikipedia
  // - [2] スペクトル法の本
  // - [3] http://pc-physics.com/syuseieuler1.html
  // - [4] http://detail.chiebukuro.yahoo.co.jp/qa/question_detail/q1091336470
  // - [5] https://www.slis.tsukuba.ac.jp/~fujisawa.makoto.fu/cgi-bin/wiki/index.php?TVD+RK
  // - [6] https://gkeyll.readthedocs.io/en/latest/dev/ssp-rk.html#ssp-rk2
  struct heun_integrator: explicit_integrator_base<heun_integrator> {
    static const int stage = 2;
    static const int order = 2;
    mutable working_buffer buffer;

    template<typename F>
    void operator()(double& time, double* value, std::size_t size, F const& f, double h) const {
      buffer.ensure<double>(2 * size);
      double* ksh_restrict k = buffer.ptr<double>();
      double* ksh_restrict x = buffer.ptr<double>() + size;

      f(k, size, time, value);
      for (std::size_t i = 0; i < size; i++) {
        x[i] = value[i] + h * k[i];
        value[i] += (1.0 / 2.0) * h * k[i];
      }

      f(k, size, time + h, x);
      for (std::size_t i = 0; i < size; i++)
        value[i] += (1.0 / 2.0) * h * k[i];

      time += h;
    }
  };

  // https://en.wikipedia.org/wiki/Heun%27s_method によるとこれも
  // Heun's method と呼ばれることがあるらしい。
  struct ralston_integrator: explicit_integrator_base<ralston_integrator> {
    static const int stage = 2;
    static const int order = 2;
    mutable working_buffer buffer;

    template<typename F>
    void operator()(double& time, double* value, std::size_t size, F const& f, double h) const {
      buffer.ensure<double>(2 * size);
      double* ksh_restrict k = buffer.ptr<double>();
      double* ksh_restrict x = buffer.ptr<double>() + size;

      f(k, size, time, value);
      for (std::size_t i = 0; i < size; i++) {
        x[i] = value[i] + (2.0 / 3.0) * h * k[i];
        value[i] += (1.0 / 4.0) * h * k[i];
      }

      f(k, size, time + (2.0 / 3.0) * h, x);
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
  struct runge3_integrator: explicit_integrator_base<runge3_integrator> {
    static const int stage = 4;
    static const int order = 3;
    mutable working_buffer buffer;

    template<typename F>
    void operator()(double& time, double* value, std::size_t size, F const& f, double h) const {
      buffer.ensure<double>(3 * size);
      double* ksh_restrict k = buffer.ptr<double>();
      double* ksh_restrict x = buffer.ptr<double>() + size;
      double* ksh_restrict y = buffer.ptr<double>() + size * 2;

      f(k, size, time, value);
      for (std::size_t i = 0; i < size; i++) {
        x[i] = value[i] + (1.0 / 2.0) * h * k[i];
        y[i] = value[i] + (1.0 / 6.0) * h * k[i];
      }

      f(k, size, time + (1.0 / 2.0) * h, x);
      for (std::size_t i = 0; i < size; i++) {
        x[i] = value[i] + h * k[i];
        y[i] += (2.0 / 3.0) * h * k[i];
      }

      f(k, size, time + h, x);
      for (std::size_t i = 0; i < size; i++)
        x[i] = value[i] + h * k[i];

      f(k, size, time + h, x);
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
  struct heun3_integrator: explicit_integrator_base<heun3_integrator> {
    static const int stage = 3;
    static const int order = 3;
    mutable working_buffer buffer;

    template<typename F>
    void operator()(double& time, double* value, std::size_t size, F const& f, double h) const {
      buffer.ensure<double>(2 * size);
      double* ksh_restrict k = buffer.ptr<double>();
      double* ksh_restrict x = buffer.ptr<double>() + size;

      f(k, size, time, value);
      for (std::size_t i = 0; i < size; i++) {
        x[i] = value[i] + (1.0 / 3.0) * h * k[i];
        value[i] += (1.0 / 4.0) * h * k[i];
      }

      f(k, size, time + (1.0 / 3.0) * h, x);
      for (std::size_t i = 0; i < size; i++)
        x[i] = 4.0 * value[i] - 3.0 * x[i] + (2.0 / 3.0) * h * k[i];

      f(k, size, time + (2.0 / 3.0) * h, x);
      for (std::size_t i = 0; i < size; i++)
        value[i] += (3.0 / 4.0) * h * k[i];

      time += h;
    }
  };

  // Ralston 3次法
  // * https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods によると、
  //   "Kutta's third-order method" である。
  struct ralston3_integrator: explicit_integrator_base<ralston3_integrator> {
    static const int stage = 3;
    static const int order = 3;
    mutable working_buffer buffer;

    template<typename F>
    void operator()(double& time, double* value, std::size_t size, F const& f, double h) const {
      buffer.ensure<double>(2 * size);
      double* ksh_restrict k = buffer.ptr<double>();
      double* ksh_restrict x = buffer.ptr<double>() + size;

      f(k, size, time, value);
      for (std::size_t i = 0; i < size; i++) {
        x[i] = value[i] + (1.0 / 2.0) * h * k[i];
        value[i] += (2.0 / 9.0) * h * k[i];
      }

      f(k, size, time + (1.0 / 2.0) * h, x);
      for (std::size_t i = 0; i < size; i++) {
        x[i] = (-4.0 / 5.0) * x[i] + (9.0 / 5.0) * value[i] + (3.0 / 4.0) * h * k[i];
        value[i] += (1.0 / 3.0) * h * k[i];
      }

      f(k, size, time + (3.0 / 4.0) * h, x);
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
  struct kutta3_integrator: explicit_integrator_base<kutta3_integrator> {
    static const int stage = 3;
    static const int order = 3;
    mutable working_buffer buffer;

    template<typename F>
    void operator()(double& time, double* value, std::size_t size, F const& f, double h) const {
      buffer.ensure<double>(2 * size);
      double* ksh_restrict k = buffer.ptr<double>();
      double* ksh_restrict x = buffer.ptr<double>() + size;

      f(k, size, time, value);
      for (std::size_t i = 0; i < size; i++) {
        x[i] = value[i] + (1.0 / 2.0) * h * k[i];
        value[i] += (1.0 / 6.0) * h * k[i];
      }

      f(k, size, time + (1.0 / 2.0) * h, x);
      for (std::size_t i = 0; i < size; i++) {
        x[i] = (-7.0 / 2.0) * x[i] + (9.0 / 2.0) * value[i] + 2.0 * h * k[i];
        value[i] += (2.0 / 3.0) * h * k[i];
      }

      f(k, size, time + h, x);
      for (std::size_t i = 0; i < size; i++)
        value[i] += (1.0 / 6.0) * h * k[i];

      time += h;
    }
  };

  // TVD-RK3, 3rd-order Strong Stability Preserving RK
  // https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods
  // https://www.slis.tsukuba.ac.jp/~fujisawa.makoto.fu/cgi-bin/wiki/index.php?TVD+RK
  struct tvdrk3_integrator: explicit_integrator_base<tvdrk3_integrator> {
    static const int stage = 3;
    static const int order = 3;
    mutable working_buffer buffer;

    template<typename F>
    void operator()(double& time, double* value, std::size_t size, F const& f, double h) const {
      buffer.ensure<double>(3 * size);
      double* ksh_restrict x = buffer.ptr<double>();
      double* ksh_restrict k1 = buffer.ptr<double>() + size;
      double* ksh_restrict k2 = buffer.ptr<double>() + size * 2;
      double* ksh_restrict& k1_plus_k2 = k1;
      double* ksh_restrict& k3 = k2;

      f(k1, size, time, value);

      for (std::size_t i = 0; i < size; i++)
        x[i] = value[i] + h * k1[i];
      f(k2, size, time + h, x);

      for (std::size_t i = 0; i < size; i++) {
        k1_plus_k2[i] = k1[i] + k2[i];
        x[i] = value[i] + 0.25 * h * k1_plus_k2[i];
      }
      f(k3, size, time + 0.5 * h, x);

      for (std::size_t i = 0; i < size; i++)
        value[i] += h * ((1.0 / 6.0) * k1_plus_k2[i] + (2.0 / 3.0) * k3[i]);

      time += h;
    }
  };

  // 4-stage 3rd-order SSP RK
  // https://gist.github.com/ketch/764e8dd4a91399eef5be
  struct tvdrk43_integrator: explicit_integrator_base<tvdrk43_integrator> {
    static const int stage = 4;
    static const int order = 3;
    mutable working_buffer buffer;

    template<typename F>
    void operator()(double& time, double* value, std::size_t size, F const& f, double h) const {
      buffer.ensure<double>(3 * size);
      double* ksh_restrict x = buffer.ptr<double>();
      double* ksh_restrict ksum = buffer.ptr<double>() + size;
      double* ksh_restrict knew = buffer.ptr<double>() + size * 2;

      f(ksum, size, time, value);

      for (std::size_t i = 0; i < size; i++)
        x[i] = value[i] + 0.5 * h * ksum[i];
      f(knew, size, time + 0.5 * h, x);

      for (std::size_t i = 0; i < size; i++) {
        ksum[i] += knew[i];
        x[i] = value[i] + 0.5 * h * ksum[i];
      }
      f(knew, size, time + h, x);

      for (std::size_t i = 0; i < size; i++) {
        ksum[i] += knew[i];
        x[i] = value[i] + (1.0 / 6.0) * h * ksum[i];
      }
      f(knew, size, time + (1.0 / 2.0) * h, x);

      for (std::size_t i = 0; i < size; i++)
        value[i] += h * ((1.0 / 6.0) * ksum[i] + (1.0 / 2.0) * knew[i]);

      time += h;
    }
  };

  //---------------------------------------------------------------------------
  // 4段4次公式

  // RK4 (classical Runge-Kutta method)。TVD RK4 はこれに等価である。
  // [1] https://www.slis.tsukuba.ac.jp/~fujisawa.makoto.fu/cgi-bin/wiki/index.php?TVD+RK#f225c799
  struct rk4_integrator: explicit_integrator_base<rk4_integrator> {
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
      f(knode, size, time, value);
      for (std::size_t i = 0; i < size; i++) {
        delta[i] = (1.0 / 6.0) * h * knode[i];
        xnode[i] = value[i] + 0.5 * h * knode[i];
      }

      // k2
      f(knode, size, time + 0.5 * h, xnode);
      for (std::size_t i = 0; i < size; i++) {
        delta[i] += (2.0 / 6.0) * h * knode[i];
        xnode[i] = value[i] + 0.5 * h * knode[i];
      }

      // k3
      f(knode, size, time + 0.5 * h, xnode);
      for (std::size_t i = 0; i < size; i++) {
        delta[i] += (2.0 / 6.0) * h * knode[i];
        xnode[i] = value[i] + h * knode[i];
      }

      // k4
      f(knode, size, time + h, xnode);
      for (std::size_t i = 0; i < size; i++) {
        double const a = delta[i] + (1.0 / 6.0) * h * knode[i];
        value[i] += a;
      }

      time += h;
    }
  };

  // Kutta 3/8-rule
  struct kutta_3_8_integrator: explicit_integrator_base<kutta_3_8_integrator> {
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
      f(k, size, time, value);
      for (std::size_t i = 0; i < size; i++) {
        xi[i] = value[i] + (1.0 / 3.0) * h * k[i];
        x4[i] = value[i] + h * k[i];
        value[i] += (1.0 / 8.0) * h * k[i];
      }

      // k2
      f(k, size, time + (1.0 / 3.0) * h, xi);
      for (std::size_t i = 0; i < size; i++) {
        xi[i] = 2.0 * xi[i] - x4[i] + h * k[i];
        x4[i] -= h * k[i];
        value[i] += (3.0 / 8.0) * h * k[i];
      }

      // k3
      f(k, size, time + (2.0 / 3.0) * h, xi);
      for (std::size_t i = 0; i < size; i++) {
        x4[i] += h * k[i];
        value[i] += (3.0 / 8.0) * h * k[i];
      }

      // k4
      f(k, size, time + h, x4);
      for (std::size_t i = 0; i < size; i++)
        value[i] += (1.0 / 8.0) * h * k[i];

      time += h;
    }
  };

  // Runge-Kutta Gill method
  //   Ref. gill.1 によると、丸め誤差の補正を但しく実行する為には計算手順があるとの事。
  //   [gill.1] [[Runge-Kutta-Gill法について - あらきけいすけの雑記帳>http://d.hatena.ne.jp/arakik10/20091004/p1]]
  struct gill_integrator: explicit_integrator_base<gill_integrator> {
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
      f(k, size, time, value);
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
      f(k, size, time + 0.5 * h, value);
      for (std::size_t i = 0; i < size; i++) {
        double const y = value[i] + alpha2 * (h * k[i] - q[i]);
        double const r = y - value[i];
        q[i] = q[i] + 3.0 * r - alpha2 * h * k[i];
        value[i] = y;
      }

      // k3
      f(k, size, time + 0.5 * h, value);
      for (std::size_t i = 0; i < size; i++) {
        double const y = value[i] + alpha3 * (h * k[i] - q[i]);
        double const r = y - value[i];
        q[i] = q[i] + 3.0 * r - alpha3 * h * k[i];
        value[i] = y;
      }

      // k4
      f(k, size, time + h, value);
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
  struct butcher5v1_integrator: explicit_integrator_base<butcher5v1_integrator> {
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
      f(k1, size, time, value);

      // k2
      static constexpr double a21 = 1.0 / 8.0;
      static constexpr double c2  = 1.0 / 8.0;
      for (std::size_t i = 0; i < size; i++)
        x[i] = value[i] + a21 * h * k1[i];
      f(k2, size, time + c2 * h, x);

      // k3
      static constexpr double a32 = 1.0 / 4.0;
      static constexpr double c3  = 1.0 / 4.0;
      for (std::size_t i = 0; i < size; i++)
        x[i] = value[i] + a32 * h * k2[i];
      f(k3, size, time + c3 * h, x);

      // k4
      static constexpr double a41 = -1.0 / 2.0;
      static constexpr double a42 =  1.0;
      static constexpr double c4  =  1.0 / 2.0;
      for (std::size_t i = 0; i < size; i++)
        x[i] = value[i] + h * (a41 * k1[i] + a42 * k2[i]);
      f(k4, size, time + c4 * h, x);

      // k5
      static constexpr double a51 = 15.0 / 16.0;
      static constexpr double a52 = -3.0 /  2.0;
      static constexpr double a53 =  3.0 /  4.0;
      static constexpr double a54 =  9.0 / 16.0;
      static constexpr double c5  =  3.0 /  4.0;
      for (std::size_t i = 0; i < size; i++)
        x[i] = value[i] + h * (a51 * k1[i] + a52 * k2[i] + a53 * k3[i] + a54 * k4[i]);
      f(k5, size, time + c5 * h, x);

      // k6
      static constexpr double a61 = -17.0 / 7.0;
      static constexpr double a62 =   4.0;
      static constexpr double a64 = -12.0 / 7.0;
      static constexpr double a65 =   8.0 / 7.0;
      static constexpr double c6  =   1.0;
      for (std::size_t i = 0; i < size; i++)
        x[i] = value[i] + h * (a61 * k1[i] + a62 * k2[i] + a64 * k4[i] + a65 * k5[i]);
      f(k6, size, time + c6 * h, x);

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
  struct butcher5v2_integrator: explicit_integrator_base<butcher5v2_integrator> {
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
      f(k1, size, time, value);

      // k2
      for (std::size_t i = 0; i < size; i++) {
        x[i]  = value[i] + a21 * h * k1[i];
        k6[i] = value[i] + a61 * h * k1[i];
        y[i]  = value[i] + b1 * h * k1[i];
      }
      f(k2, size, time + c2 * h, x);

      // k3
      for (std::size_t i = 0; i < size; i++) {
        x[i] = value[i] + h * (a31 * k1[i] + a32 * k2[i]);
        k6[i] += a62 * h * k2[i];
      }
      f(k3, size, time + c3 * h, x);

      // k4
      for (std::size_t i = 0; i < size; i++) {
        x[i] = a42 / a32 * x[i] + (1.0 - a42 / a32) * value[i] - ((a42 / a32) * a31) * h * k1[i] + a43 * h * k3[i];
        k6[i] += a63 * h * k3[i];
        y[i] += b3 * h * k3[i];
      }
      f(k4, size, time + c4 * h, x);

      // k5
      for (std::size_t i = 0; i < size; i++) {
        x[i] = value[i] + h * (a51 * k1[i] + a54 * k4[i]);
        k6[i] += a64 * h * k4[i];
        y[i] += b4 * h * k4[i];
      }
      f(k5, size, time + c5 * h, x);

      // k6
      for (std::size_t i = 0; i < size; i++) {
        x[i] = k6[i] + a65 * h * k5[i];
        y[i] += b5 * h * k5[i];
      }
      f(k6, size, time + c6 * h, x);

      // increment
      for (std::size_t i = 0; i < size; i++)
        value[i] = y[i] + b6 * h * k6[i];

      time += h;

    }
  };

  // http://www.330k.info/essay/Explicit-Runge-Kutta-Butcher-Tableau
  struct butcher5v3_integrator: explicit_integrator_base<butcher5v3_integrator> {
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
      f(k1, size, time, value);

      // k2
      static constexpr double a21 = -1.0 / 2.0;
      static constexpr double c2  = -1.0 / 2.0;
      for (std::size_t i = 0; i < size; i++)
        x[i] = value[i] + a21 * h * k1[i];
      f(k2, size, time + c2 * h, x);

      // k3
      static constexpr double a31 =  5.0 / 16.0;
      static constexpr double a32 = -1.0 / 16.0;
      static constexpr double c3  =  1.0 /  4.0;
      for (std::size_t i = 0; i < size; i++)
        x[i] = value[i] + h * (a31 * k1[i] + a32 * k2[i]);
      f(k3, size, time + c3 * h, x);

      // k4
      static constexpr double a41 = -3.0 / 4.0;
      static constexpr double a42 =  1.0 / 4.0;
      static constexpr double a43 =  1.0;
      static constexpr double c4  =  1.0 / 2.0;
      for (std::size_t i = 0; i < size; i++)
        x[i] = value[i] + h * (a41 * k1[i] + a42 * k2[i] + a43 * k3[i]);
      f(k4, size, time + c4 * h, x);

      // k5
      static constexpr double a51 = 3.0 / 16.0;
      static constexpr double a54 = 9.0 / 16.0;
      static constexpr double c5  = 3.0 /  4.0;
      for (std::size_t i = 0; i < size; i++)
        x[i] = value[i] + h * (a51 * k1[i] + a54 * k4[i]);
      f(k5, size, time + c5 * h, x);

      // k6
      static constexpr double a62 =  -1.0 / 7.0;
      static constexpr double a63 =  12.0 / 7.0;
      static constexpr double a64 = -12.0 / 7.0;
      static constexpr double a65 =   8.0 / 7.0;
      static constexpr double c6  = 1.0;
      for (std::size_t i = 0; i < size; i++)
        x[i] = value[i] + h * (a62 * k2[i] + a63 * k3[i] + a64 * k4[i] + a65 * k5[i]);
      f(k6, size, time + c6 * h, x);

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
  struct hammud6_integrator: explicit_integrator_base<hammud6_integrator> {
    static const int stage = 7;
    static const int order = 6;
    mutable working_buffer buffer;

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
      f(k1, size, time, value);

      // k2
      static constexpr double a21 = 4.0 / 7.0;
      static constexpr double c2  = 4.0 / 7.0;
      for (std::size_t i = 0; i < size; i++)
        x[i] = value[i] + a21 * h * k1[i];
      f(k2, size, time + c2 * h, x);

      // k3
      static constexpr double a31 = 115.0 / 112.0;
      static constexpr double a32 =  -5.0 /  16.0;
      static constexpr double c3  =   5.0 /   7.0;
      for (std::size_t i = 0; i < size; i++)
        x[i] = value[i] + h * (a31 * k1[i] + a32 * k2[i]);
      f(k3, size, time + c3 * h, x);

      // k4
      static constexpr double a41 = 589.0 / 630.0;
      static constexpr double a42 =   5.0 /  18.0;
      static constexpr double a43 = -16.0 /  45.0;
      static constexpr double c4  =   6.0 /   7.0;
      for (std::size_t i = 0; i < size; i++)
        x[i] = value[i] + h * (a41 * k1[i] + a42 * k2[i] + a43 * k3[i]);
      f(k4, size, time + c4 * h, x);

      // k5
      static constexpr double a51 = 229.0 / 1200.0 -  29.0 * sqrt5 / 6000.0;
      static constexpr double a52 = 119.0 /  240.0 - 187.0 * sqrt5 / 1200.0;
      static constexpr double a53 = -14.0 /   75.0 +  34.0 * sqrt5 /  375.0;
      static constexpr double a54 =                   -3.0 * sqrt5 /  100.0;
      static constexpr double c5  = (5.0 - sqrt5) / 10.0;
      for (std::size_t i = 0; i < size; i++)
        x[i] = value[i] + h * (a51 * k1[i] + a52 * k2[i] + a53 * k3[i] + a54 * k4[i]);
      f(k5, size, time + c5 * h, x);

      // k6
      static constexpr double a61 =  71.0 / 2400.0 - 587.0 * sqrt5 / 12000.0;
      static constexpr double a62 = 187.0 / 480.0  - 391.0 * sqrt5 /  2400.0;
      static constexpr double a63 = -38.0 / 75.0   +  26.0 * sqrt5 /   375.0;
      static constexpr double a64 =  27.0 / 80.0   -   3.0 * sqrt5 /   400.0;
      static constexpr double a65 =   1.0 / 4.0    +         sqrt5 /     4.0;
      static constexpr double c6  = (5.0 + sqrt5) / 10.0;
      for (std::size_t i = 0; i < size; i++)
        x[i] = value[i] + h * (a61 * k1[i] + a62 * k2[i] + a63 * k3[i] + a64 * k4[i] + a65 * k5[i]);
      f(k6, size, time + c6 * h, x);

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
      f(k7, size, time + c7 * h, x);

      // increment
      for (std::size_t i = 0; i < size; i++)
        value[i] += h * (b1 * k1[i] + b5 * k5[i] + b6 * k6[i] + b7 * k7[i]);

      time += h;
    }
  };

  // Butcher's 6th-order method
  // http://www.mymathlib.com/diffeq/runge-kutta/runge_kutta_butcher.html
  // http://www.mymathlib.com/c_source/diffeq/runge_kutta/runge_kutta_butcher.c
  struct butcher6_integrator: explicit_integrator_base<butcher6_integrator> {
    static const int stage = 7;
    static const int order = 6;
    mutable working_buffer buffer;

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

      static constexpr double b1 =  11.0 / 120.0;
      static constexpr double b3 =  81.0 / 120.0;
      static constexpr double b4 =  81.0 / 120.0;
      static constexpr double b5 = -32.0 / 120.0;
      static constexpr double b6 = -32.0 / 120.0;
      static constexpr double b7 =  11.0 / 120.0;

      // k1
      f(k1, size, time, value);

      // k2
      static constexpr double a21 = 1.0 / 3.0;
      static constexpr double c2  = 1.0 / 3.0;
      for (std::size_t i = 0; i < size; i++)
        x[i] = value[i] + a21 * h * k1[i];
      f(k2, size, time + c2 * h, x);

      // k3
      static constexpr double a32 = 2.0 / 3.0;
      static constexpr double c3  = 2.0 / 3.0;
      for (std::size_t i = 0; i < size; i++)
        x[i] = value[i] + a32 * h * k2[i];
      f(k3, size, time + c3 * h, x);

      // k4
      static constexpr double a41 =  1.0 / 12.0;
      static constexpr double a42 =  4.0 / 12.0;
      static constexpr double a43 = -1.0 / 12.0;
      static constexpr double c4  =  1.0 /  3.0;
      for (std::size_t i = 0; i < size; i++)
        x[i] = value[i] + h * (a41 * k1[i] + a42 * k2[i] + a43 * k3[i]);
      f(k4, size, time + c4 * h, x);

      // k5
      static constexpr double a51 = -1.0 / 16.0;
      static constexpr double a52 = 18.0 / 16.0;
      static constexpr double a53 = -3.0 / 16.0;
      static constexpr double a54 = -6.0 / 16.0;
      static constexpr double c5  =  1.0 /  2.0;
      for (std::size_t i = 0; i < size; i++)
        x[i] = value[i] + h * (a51 * k1[i] + a52 * k2[i] + a53 * k3[i] + a54 * k4[i]);
      f(k5, size, time + c5 * h, x);

      // k6
      static constexpr double a62 =  9.0 / 8.0;
      static constexpr double a63 = -3.0 / 8.0;
      static constexpr double a64 = -6.0 / 8.0;
      static constexpr double a65 =  4.0 / 8.0;
      static constexpr double c6  =  1.0 / 2.0;
      for (std::size_t i = 0; i < size; i++)
        x[i] = value[i] + h * (a62 * k2[i] + a63 * k3[i] + a64 * k4[i] + a65 * k5[i]);
      f(k6, size, time + c6 * h, x);

      // k7 <= k2
      static constexpr double a71 =   9.0 / 44.0;
      static constexpr double a72 = -36.0 / 44.0;
      static constexpr double a73 =  63.0 / 44.0;
      static constexpr double a74 =  72.0 / 44.0;
      static constexpr double a75 = -64.0 / 44.0;
      static constexpr double c7 = 1.0;
      for (std::size_t i = 0; i < size; i++)
        x[i] = value[i] + h * (a71 * k1[i] + a72 * k2[i] + a73 * k3[i] + a74 * k4[i] + a75 * k5[i]);
      f(k7, size, time + c7 * h, x);

      // increment
      for (std::size_t i = 0; i < size; i++)
        value[i] += h * (b1 * k1[i] + b3 * k3[i] + b4 * k4[i] + b5 * k5[i] + b6 * k6[i] + b7 * k7[i]);

      time += h;
    }
  };

  // http://www.330k.info/essay/Explicit-Runge-Kutta-Butcher-Tableau
  struct shanks7_integrator: explicit_integrator_base<shanks7_integrator> {
    static const int stage = 9;
    static const int order = 7;
    mutable working_buffer buffer;

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
      f(k1, size, time, value);

      // k2
      static constexpr double a21 = 2.0 / 9.0;
      static constexpr double c2  = 2.0 / 9.0;
      for (std::size_t i = 0; i < size; i++)
        x[i] = value[i] + a21 * h * k1[i];
      f(k2, size, time + c2 * h, x);

      // k3
      static constexpr double a31 = 1.0 / 12.0;
      static constexpr double a32 = 1.0 / 4.0;
      static constexpr double c3  = 1.0 / 3.0;
      for (std::size_t i = 0; i < size; i++)
        x[i] = value[i] + h * (a31 * k1[i] + a32 * k2[i]);
      f(k3, size, time + c3 * h, x);

      // k4
      static constexpr double a41 = 1.0 / 8.0;
      static constexpr double a43 = 3.0 / 8.0;
      static constexpr double c4  = 1.0 / 2.0;
      for (std::size_t i = 0; i < size; i++)
        x[i] = value[i] + h * (a41 * k1[i] + a43 * k3[i]);
      f(k4, size, time + c4 * h, x);

      // k5
      static constexpr double a51 = 23.0 / 216.0;
      static constexpr double a53 =  7.0 / 72.0;
      static constexpr double a54 = -1.0 / 27.0;
      static constexpr double c5  =  1.0 / 6.0;
      for (std::size_t i = 0; i < size; i++)
        x[i] = value[i] + h * (a51 * k1[i] + a53 * k3[i] + a54 * k4[i]);
      f(k5, size, time + c5 * h, x);

      // k6
      static constexpr double a61 = -4136.0 / 729.0;
      static constexpr double a63 = -4528.0 / 243.0;
      static constexpr double a64 =  5264.0 / 729.0;
      static constexpr double a65 =  1456.0 / 81.0;
      static constexpr double c6  =     8.0 / 9.0;
      for (std::size_t i = 0; i < size; i++)
        x[i] = value[i] + h * (a61 * k1[i] + a63 * k3[i] + a64 * k4[i] + a65 * k5[i]);
      f(k6, size, time + c6 * h, x);

      // k7
      static constexpr double a71 = 8087.0 / 11664.0;
      static constexpr double a73 =  484.0 / 243.0;
      static constexpr double a74 = -518.0 / 729.0;
      static constexpr double a75 = -658.0 / 351.0;
      static constexpr double a76 =    7.0 / 624.0;
      static constexpr double c7  =    1.0 / 9.0;
      for (std::size_t i = 0; i < size; i++)
        x[i] = value[i] + h * (a71 * k1[i] + a73 * k3[i] + a74 * k4[i] + a75 * k5[i] + a76 * k6[i]);
      f(k7, size, time + c7 * h, x);

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
      f(k8, size, time + c8 * h, x);

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
      f(k9, size, time + c9 * h, x);

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
  struct cooper_verner7_integrator: explicit_integrator_base<cooper_verner7_integrator> {
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
      f(k1, size, time, value);

      // k2
      static constexpr double a21 = (7.0 + 1.0 * sqrt21) / 42.0;
      static constexpr double c2  = (7.0 +       sqrt21) / 42.0;
      for (std::size_t i = 0; i < size; i++)
        x[i] = value[i] + a21 * h * k1[i];
      f(k2, size, time + c2 * h, x);

      // k3
      static constexpr double a32 = (7.0 + 1.0 * sqrt21) / 21.0;
      static constexpr double c3  = (7.0 +       sqrt21) / 21.0;
      for (std::size_t i = 0; i < size; i++)
        x[i] = value[i] + a32 * h * k2[i];
      f(k3, size, time + c3 * h, x);

      // k4
      static constexpr double a41 = ( 7.0 + 1.0 * sqrt21) / 56.0;
      static constexpr double a43 = (21.0 + 3.0 * sqrt21) / 56.0;
      static constexpr double c4  = ( 7.0 +       sqrt21) / 14.0;
      for (std::size_t i = 0; i < size; i++)
        x[i] = value[i] + h * (a41 * k1[i] + a43 * k3[i]);
      f(k4, size, time + c4 * h, x);

      // k5
      static constexpr double a51 = (  8.0 - 1.0 * sqrt21) / 16.0;
      static constexpr double a53 = (-21.0 + 6.0 * sqrt21) / 16.0;
      static constexpr double a54 = ( 21.0 - 5.0 * sqrt21) / 16.0;
      static constexpr double c5 = 1.0 / 2.0;
      for (std::size_t i = 0; i < size; i++)
        x[i] = value[i] + h * (a51 * k1[i] + a53 * k3[i] + a54 * k4[i]);
      f(k5, size, time + c5 * h, x);

      // k6
      static constexpr double a61 = (-1687.0 + 374.0 * sqrt21) / 196.0;
      static constexpr double a63 = (  969.0 - 210.0 * sqrt21) /  28.0;
      static constexpr double a64 = ( -381.0 +  83.0 * sqrt21) /  14.0;
      static constexpr double a65 = (   84.0 -  20.0 * sqrt21) /  49.0;
      static constexpr double c6  = (    7.0 -         sqrt21) /  14.0;
      for (std::size_t i = 0; i < size; i++)
        x[i] = value[i] + h * (a61 * k1[i] + a63 * k3[i] + a64 * k4[i] + a65 * k5[i]);
      f(k6, size, time + c6 * h, x);

      // k7
      static constexpr double a71 = (  583.0 - 131.0 * sqrt21) / 128.0;
      static constexpr double a73 = (-2373.0 + 501.0 * sqrt21) / 128.0;
      static constexpr double a74 = ( 4221.0 - 914.0 * sqrt21) / 288.0;
      static constexpr double a75 = (   -9.0 +   4.0 * sqrt21) /  18.0;
      static constexpr double a76 = (  189.0 +  35.0 * sqrt21) / 576.0;
      static constexpr double c7 = 1.0 / 2.0;
      for (std::size_t i = 0; i < size; i++)
        x[i] = value[i] + h * (a71 * k1[i] + a73 * k3[i] + a74 * k4[i] + a75 * k5[i] + a76 * k6[i]);
      f(k7, size, time + c7 * h, x);

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
      f(k8, size, time + c8 * h, x);

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
      f(k9, size, time + c9 * h, x);

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
  struct cooper_verner8_integrator: explicit_integrator_base<cooper_verner8_integrator> {
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
      f(k1, size, time, value);

      // k2
      static constexpr double a21 = 0.5;
      static constexpr double c20 = a21;
      for (std::size_t i = 0; i < size; i++) {
        delta[i] = b1 * h * k1[i];
        xnode[i] = value[i] + a21 * h * k1[i];
      }
      f(k2, size, time + c20 * h, xnode);

      // k3
      static constexpr double a31 = 0.25;
      static constexpr double a32 = 0.25;
      static constexpr double c30 = a31 + a32;
      for (std::size_t i = 0; i < size; i++)
        xnode[i] = value[i] + a31 * h * k1[i] + a32 * h * k2[i];
      f(k3, size, time + c30 * h, xnode);

      // k4 <= k2
      static constexpr double a41 = (1.0 /  7.0);
      static constexpr double a42 = (1.0 / 98.0) * (-7 - 3 * sqrt21);
      static constexpr double a43 = (1.0 / 49.0) * (21 + 5 * sqrt21);
      static constexpr double c40 = a41 + a42 + a43;
      for (std::size_t i = 0; i < size; i++)
        xnode[i] = value[i] + a41 * h * k1[i] + a42 * h * k2[i] + a43 * h * k3[i];
      f(k4, size, time + c40 * h, xnode);

      // k5
      static constexpr double a51 = (1.0 /  84.0) * (11 + 1 * sqrt21);
      static constexpr double a53 = (1.0 /  63.0) * (18 + 4 * sqrt21);
      static constexpr double a54 = (1.0 / 252.0) * (21 - 1 * sqrt21);
      static constexpr double c50 = a51 + a53 + a54;
      for (std::size_t i = 0; i < size; i++)
        xnode[i] = value[i] + a51 * h * k1[i] + a53 * h * k3[i] + a54 * h * k4[i];
      f(k5, size, time + c50 * h, xnode);

      // k6
      static constexpr double a61 = (1.0 /  48.0) * (   5 +  1 * sqrt21);
      static constexpr double a63 = (1.0 /  36.0) * (   9 +  1 * sqrt21);
      static constexpr double a64 = (1.0 / 360.0) * (-231 + 14 * sqrt21);
      static constexpr double a65 = (1.0 /  80.0) * (  63 -  7 * sqrt21);
      static constexpr double c60 = a61 + a63 + a64 + a65;
      for (std::size_t i = 0; i < size; i++)
        xnode[i] = value[i] + a61 * h * k1[i] + a63 * h * k3[i] + a64 * h * k4[i] + a65 * h * k5[i];
      f(k6, size, time + c60 * h, xnode);

      // k7 <= k3
      static constexpr double a71 = (1.0 /  42.0) * (  10 -   1 * sqrt21);
      static constexpr double a73 = (1.0 / 315.0) * (-432 +  92 * sqrt21);
      static constexpr double a74 = (1.0 /  90.0) * ( 633 - 145 * sqrt21);
      static constexpr double a75 = (1.0 /  70.0) * (-504 + 115 * sqrt21);
      static constexpr double a76 = (1.0 /  35.0) * (  63 -  13 * sqrt21);
      static constexpr double c70 = a71 + a73 + a74 + a75 + a76;
      for (std::size_t i = 0; i < size; i++)
        xnode[i] = value[i] + a71 * h * k1[i] + a73 * h * k3[i] + a74 * h * k4[i] + a75 * h * k5[i] + a76 * h * k6[i];
      f(k7, size, time + c70 * h, xnode);

      // k8 <= k4
      static constexpr double a81 =  1.0 /  14.0;
      static constexpr double a85 = (1.0 / 126.0) * (14 - 3 * sqrt21);
      static constexpr double a86 = (1.0 /  63.0) * (13 - 3 * sqrt21);
      static constexpr double a87 =  1.0 /   9.0;
      static constexpr double c80 = a81 + a85 + a86 + a87;
      for (std::size_t i = 0; i < size; i++)
        xnode[i] = value[i] + a81 * h * k1[i] + a85 * h * k5[i] + a86 * h * k6[i] + a87 * h * k7[i];
      f(k8, size, time + c80 * h, xnode);

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
      f(k9, size, time + c90 * h, xnode);

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
      f(kA, size, time + cA0 * h, xnode);

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
      f(kB, size, time + cB0 * h, xnode);

      // increment
      for (std::size_t i = 0; i < size; i++) {
        double const a = delta[i] + bB * h * kB[i];
        value[i] += a;
      }

      time += h;
    }
  };


  // Verner 9th order (埋込み型RKの一部らしい)
  // - https://www.330k.info/essay/explicit-runge-kutta-butcher-tableau/
  //   この方法はとても大きな誤差の係数を持つ。恐らく kC の係数が巨大な為。
  struct verner9_integrator: explicit_integrator_base<verner9_integrator> {
    static const int stage = 15;
    static const int order = 9;
    mutable working_buffer buffer;

    static constexpr double sqrt6 = std::sqrt(6.0);

    template<typename F>
    void operator()(double& time, double* ksh_restrict value, std::size_t size, F const& f, double h) const {
      buffer.ensure<double>(11 * size);
      double* ksh_restrict  xnode = buffer.ptr<double>();
      double* ksh_restrict  k1    = buffer.ptr<double>() + size * 1;
      double* ksh_restrict  k2    = buffer.ptr<double>() + size * 2;
      double* ksh_restrict& k3    = k2;
      double* ksh_restrict  k4    = buffer.ptr<double>() + size * 3;
      double* ksh_restrict& k5    = k3;
      double* ksh_restrict  k6    = buffer.ptr<double>() + size * 4;
      double* ksh_restrict& k7    = k4;
      double* ksh_restrict& k8    = k5;
      double* ksh_restrict  k9    = buffer.ptr<double>() + size * 5;
      double* ksh_restrict  kA    = buffer.ptr<double>() + size * 6;
      double* ksh_restrict  kB    = buffer.ptr<double>() + size * 7;
      double* ksh_restrict  kC    = buffer.ptr<double>() + size * 8;
      double* ksh_restrict  kD    = buffer.ptr<double>() + size * 9;
      // Note: kE は埋め込み型の誤差評価の為に使われるので跳んでいる。
      double* ksh_restrict  kF    = buffer.ptr<double>() + size * 10;
      double* ksh_restrict& kG    = k6;

      static constexpr double b1 =    23.0 /   525.0;
      static constexpr double b8 =   171.0 /  1400.0;
      static constexpr double b9 =    86.0 /   525.0;
      static constexpr double bA =    93.0 /   280.0;
      static constexpr double bB = -2048.0 /  6825.0;
      static constexpr double bC =    -3.0 / 18200.0;
      static constexpr double bD =    39.0 /   175.0;
      static constexpr double bF =     9.0 /    25.0;
      static constexpr double bG =   233.0 /  4200.0;

      // k1
      f(k1, size, time, value);

      // k2
      static constexpr double a21 = 1.0 / 12.0;
      static constexpr double c20 = a21;
      for (std::size_t i = 0; i < size; i++)
        xnode[i] = value[i] + a21 * h * k1[i];
      f(k2, size, time + c20 * h, xnode);

      // k3
      static constexpr double a31 = 1.0 / 27.0;
      static constexpr double a32 = 2.0 / 27.0;
      static constexpr double c30 = a31 + a32;
      for (std::size_t i = 0; i < size; i++)
        xnode[i] = value[i] + a31 * h * k1[i] + a32 * h * k2[i];
      f(k3, size, time + c30 * h, xnode);

      // k4
      static constexpr double a41 = 1.0 / 24.0;
      static constexpr double a43 = 1.0 / 8.0;
      static constexpr double c40 = a41 + a43;
      for (std::size_t i = 0; i < size; i++)
        xnode[i] = value[i] + a41 * h * k1[i] + a43 * h * k3[i];
      f(k4, size, time + c40 * h, xnode);

      // k5
      static constexpr double a51 = (1.0 / 375.0) * (   4.0 +  94.0 * sqrt6);
      static constexpr double a53 = (1.0 / 375.0) * (-282.0 - 252.0 * sqrt6);
      static constexpr double a54 = (1.0 / 375.0) * ( 328.0 + 208.0 * sqrt6);
      static constexpr double c50 = a51 + a53 + a54;
      for (std::size_t i = 0; i < size; i++)
        xnode[i] = value[i] + a51 * h * k1[i] + a53 * h * k3[i] + a54 * h * k4[i];
      f(k5, size, time + c50 * h, xnode);

      // k6
      static constexpr double a61 = (1.0 /  150.0) * (  9.0 -        sqrt6);
      static constexpr double a64 = (1.0 / 1425.0) * (312.0 + 32.0 * sqrt6);
      static constexpr double a65 = (1.0 /  570.0) * ( 69.0 + 29.0 * sqrt6);
      static constexpr double c60 = a61 + a64 + a65;
      for (std::size_t i = 0; i < size; i++)
        xnode[i] = value[i] + a61 * h * k1[i] + a64 * h * k4[i] + a65 * h * k5[i];
      f(k6, size, time + c60 * h, xnode);

      // k7
      static constexpr double a71 = (1.0 / 1250.0) * (   927.0 -  347.0 * sqrt6);
      static constexpr double a74 = (1.0 / 9375.0) * (-16248.0 + 7328.0 * sqrt6);
      static constexpr double a75 = (1.0 / 3750.0) * (  -489.0 +  179.0 * sqrt6);
      static constexpr double a76 = (1.0 / 9375.0) * ( 14268.0 - 5798.0 * sqrt6);
      static constexpr double c70 = a71 + a74 + a75 + a76;
      for (std::size_t i = 0; i < size; i++)
        xnode[i] = value[i] + a71 * h * k1[i] + a74 * h * k4[i] + a75 * h * k5[i] + a76 * h * k6[i];
      f(k7, size, time + c70 * h, xnode);

      // k8
      static constexpr double a81 =  2.0 / 27.0;
      static constexpr double a86 = (1.0 /  54.0) * (16.0 - sqrt6);
      static constexpr double a87 = (1.0 /  54.0) * (16.0 + sqrt6);
      static constexpr double c80 = a81 + a86 + a87;
      for (std::size_t i = 0; i < size; i++)
        xnode[i] = value[i] + a81 * h * k1[i] + a86 * h * k6[i] + a87 * h * k7[i];
      f(k8, size, time + c80 * h, xnode);

      // k9
      static constexpr double a91 =  19.0 / 256.0;
      static constexpr double a96 = ( 1.0 / 512.0) * (118.0 - 23.0 * sqrt6);
      static constexpr double a97 = ( 1.0 / 512.0) * (118.0 + 23.0 * sqrt6);
      static constexpr double a98 =  -9.0 / 256.0;
      static constexpr double c90 = a91 + a96 + a97 + a98;
      for (std::size_t i = 0; i < size; i++)
        xnode[i] = value[i] + a91 * h * k1[i] + a96 * h * k6[i] + a97 * h * k7[i] + a98 * h * k8[i];
      f(k9, size, time + c90 * h, xnode);

      // kA
      static constexpr double aA1 =  11.0 / 144.0;
      static constexpr double aA6 = ( 1.0 / 864.0) * (266.0 - sqrt6);
      static constexpr double aA7 = ( 1.0 / 864.0) * (266.0 + sqrt6);
      static constexpr double aA8 =  -1.0 /  16.0;
      static constexpr double aA9 =  -8.0 /  27.0;
      static constexpr double cA0 = aA1 + aA6 + aA7 + aA8 + aA9;
      for (std::size_t i = 0; i < size; i++)
        xnode[i] = value[i] + aA1 * h * k1[i] + aA6 * h * k6[i] + aA7 * h * k7[i] + aA8 * h * k8[i] + aA9 * h * k9[i];
      f(kA, size, time + cA0 * h, xnode);

      // kB
      static constexpr double aB1 = ( 5034.0 -  271.0 * sqrt6) / 61440.0;
      static constexpr double aB7 = ( 7859.0 - 1626.0 * sqrt6) / 10240.0;
      static constexpr double aB8 = (-2232.0 +  813.0 * sqrt6) / 20480.0;
      static constexpr double aB9 = ( -594.0 +  271.0 * sqrt6) /   960.0;
      static constexpr double aBA = (  657.0 -  813.0 * sqrt6) /  5120.0;
      static constexpr double cB0 = aB1 + aB7 + aB8 + aB9 + aBA;
      for (std::size_t i = 0; i < size; i++)
        xnode[i] = value[i] + aB1 * h * k1[i] + aB7 * h * k7[i] + aB8 * h * k8[i] + aB9 * h * k9[i] + aBA * h * kA[i];
      f(kB, size, time + cB0 * h, xnode);

      // kC
      static constexpr double aC1 = (   5996.0 -   3794.0 * sqrt6) /   405.0;
      static constexpr double aC6 = (  -4342.0 -    338.0 * sqrt6) /     9.0; // ~ -574
      static constexpr double aC7 = ( 154922.0 -  40458.0 * sqrt6) /   135.0; // ~ +413
      static constexpr double aC8 = (  -4176.0 +   3794.0 * sqrt6) /    45.0; // ~ +114
      static constexpr double aC9 = (-340864.0 + 242816.0 * sqrt6) /   405.0; // ~ +627
      static constexpr double aCA = (  26304.0 -  15176.0 * sqrt6) /    45.0; // ~ -241
      static constexpr double aCB = ( -26624.0                   ) /    81.0; // ~ -328
      static constexpr double cC0 = aC1 + aC6 + aC7 + aC8 + aC9 + aCA + aCB;
      for (std::size_t i = 0; i < size; i++)
        xnode[i] = value[i] + aC1 * h * k1[i] + aC6 * h * k6[i] + aC7 * h * k7[i] + aC8 * h * k8[i] + aC9 * h * k9[i] + aCA * h * kA[i] + aCB * h * kB[i];
      f(kC, size, time + cC0 * h, xnode);

      // kD
      static constexpr double aD1 = (   3793.0 +   2168.0 * sqrt6) / 103680.0;
      static constexpr double aD6 = (   4042.0 +   2263.0 * sqrt6) /  13824.0;
      static constexpr double aD7 = (-231278.0 +  40717.0 * sqrt6) /  69120.0;
      static constexpr double aD8 = (   7947.0 -   2168.0 * sqrt6) /  11520.0;
      static constexpr double aD9 = (   1048.0 -    542.0 * sqrt6) /    405.0;
      static constexpr double aDA = (  -1383.0 +    542.0 * sqrt6) /    720.0;
      static constexpr double aDB = (   2624.0                   ) /   1053.0;
      static constexpr double aDC = (      3.0                   ) /   1664.0;
      static constexpr double cD0 = aD1 + aD6 + aD7 + aD8 + aD9 + aDA + aDB + aDC;
      for (std::size_t i = 0; i < size; i++)
        xnode[i] = value[i] + aD1 * h * k1[i] + aD6 * h * k6[i] + aD7 * h * k7[i] + aD8 * h * k8[i] + aD9 * h * k9[i] + aDA * h * kA[i] + aDB * h * kB[i] + aDC * h * kC[i];
      f(kD, size, time + cD0 * h, xnode);

      // // kE
      // static constexpr double aE1 = (   -137.0                   ) /   1296.0;
      // static constexpr double aE6 = (   5642.0 -    337.0 * sqrt6) /    864.0;
      // static constexpr double aE7 = (   5642.0 +    337.0 * sqrt6) /    864.0;
      // static constexpr double aE8 = (   -299.0                   ) /     48.0;
      // static constexpr double aE9 = (    184.0                   ) /     81.0;
      // static constexpr double aEA = (    -44.0                   ) /      9.0;
      // static constexpr double aEB = (  -5120.0                   ) /   1053.0;
      // static constexpr double aEC = (    -11.0                   ) /    468.0;
      // static constexpr double aED = (     16.0                   ) /      9.0;
      // static constexpr double cE0 = aE1 + aE6 + aE7 + aE8 + aE9 + aEA + aEB + aEC + aED;
      // for (std::size_t i = 0; i < size; i++)
      //   xnode[i] = value[i] + aE1 * h * k1[i] + aE6 * h * k6[i] + aE7 * h * k7[i] + aE8 * h * k8[i] + aE9 * h * k9[i] + aEA * h * kA[i] + aEB * h * kB[i] + aEC * h * kC[i] + aED * h * kD[i];
      // f(kE, size, time + cE0 * h, xnode);

      // kF
      static constexpr double aF1 = (  33617.0 -   2168.0 * sqrt6) / 518400.0;
      static constexpr double aF6 = (  -3846.0 +     31.0 * sqrt6) /  13824.0;
      static constexpr double aF7 = ( 155338.0 -  52807.0 * sqrt6) / 345600.0;
      static constexpr double aF8 = ( -12537.0 +   2168.0 * sqrt6) /  57600.0;
      static constexpr double aF9 = (     92.0 +    542.0 * sqrt6) /   2025.0;
      static constexpr double aFA = (  -1797.0 -    542.0 * sqrt6) /   3600.0;
      static constexpr double aFB = (    320.0                   ) /    567.0;
      static constexpr double aFC = (     -1.0                   ) /   1920.0;
      static constexpr double aFD = (      4.0                   ) /    105.0;
      static constexpr double cF0 = aF1 + aF6 + aF7 + aF8 + aF9 + aFA + aFB + aFC + aFD;
      for (std::size_t i = 0; i < size; i++)
        xnode[i] = value[i] + aF1 * h * k1[i] + aF6 * h * k6[i] + aF7 * h * k7[i] + aF8 * h * k8[i] + aF9 * h * k9[i] + aFA * h * kA[i] + aFB * h * kB[i] + aFC * h * kC[i] + aFD * h * kD[i];
      f(kF, size, time + cF0 * h, xnode);

      // kG
      static constexpr double aG1 = ( -36487.0 -  30352.0 * sqrt6) / 279600.0;
      static constexpr double aG6 = ( -29666.0 -   4499.0 * sqrt6) /   7456.0;
      static constexpr double aG7 = (2779182.0 - 615973.0 * sqrt6) / 186400.0;
      static constexpr double aG8 = ( -94329.0 +  91056.0 * sqrt6) /  93200.0;
      static constexpr double aG9 = (-232192.0 + 121408.0 * sqrt6) /  17475.0;
      static constexpr double aGA = ( 101226.0 -  22764.0 * sqrt6) /   5825.0;
      static constexpr double aGB = (-169984.0                   ) /   9087.0;
      static constexpr double aGC = (    -87.0                   ) /  30290.0;
      static constexpr double aGD = (    492.0                   ) /   1165.0;
      static constexpr double aGF = (   1260.0                   ) /    233.0;
      static constexpr double cG0 = aG1 + aG6 + aG7 + aG8 + aG9 + aGA + aGB + aGC + aGD + aGF;
      for (std::size_t i = 0; i < size; i++)
        xnode[i] = value[i] + aG1 * h * k1[i] + aG6 * h * k6[i] + aG7 * h * k7[i] + aG8 * h * k8[i] + aG9 * h * k9[i] + aGA * h * kA[i] + aGB * h * kB[i] + aGC * h * kC[i] + aGD * h * kD[i] + aGF * h * kF[i];
      f(kG, size, time + cG0 * h, xnode);

      // increment
      for (std::size_t i = 0; i < size; i++)
        value[i] += h * (b1 * k1[i] + b8 * k8[i] + b9 * k9[i] + bA * kA[i] + bB * kB[i] + bC * kC[i] + bD * kD[i] + bF * kF[i] + bG * kG[i]);

      time += h;
    }
  };

  // 他
  // * http://www.mymathlib.com/diffeq/runge-kutta/runge_kutta_ralston_4.html Ralston's 4th
  // * ハイラーの Ralson(1962), Hull(1967) (u, v) = (0.4, 0.45) と同じ物か?
  // * http://www.mymathlib.com/diffeq/runge-kutta/runge_kutta_nystrom.html Nystrom's 5th
  // * https://gist.github.com/ketch/764e8dd4a91399eef5be このページには他にも沢山のステージの SSPRK が載っている。

}
}

#endif
