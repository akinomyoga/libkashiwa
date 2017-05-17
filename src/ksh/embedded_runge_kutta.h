// -*- mode:c++ -*-
#ifndef KASHIWA_EMBEDDED_RUNGE_KUTTA_H
#define KASHIWA_EMBEDDED_RUNGE_KUTTA_H
#include <cstddef>
#include <cfloat>
#include <mwg/except.h>
#include "def.h"
#include "buffer.h"

namespace kashiwa {
namespace runge_kutta {

  struct iequation_for_erk;

  struct stat_t {
    int nfcn   {0};
    int nstep  {0};
    int naccpt {0};
    int nrejct {0};

    iequation_for_erk* eq            {nullptr};
    double const*      previousValue {nullptr};
    std::size_t        previousSize  {0};
    double             previousTime  {0.0};
    double             previousStep  {0.0};
  };

  struct idense_output {
    virtual void get_values_at(
      double* ksh_restrict interpolated, double time,
      std::size_t const* icomp, std::size_t ncomp
    ) = 0;

    void get_values_at(double* ksh_restrict interpolated, double time) {
      return get_values_at(interpolated, time, nullptr, 0);
    }
  };

  struct iequation_for_erk {
    virtual void eval_f(double* ksh_restrict slope, double t, double const* ksh_restrict value) = 0;
    virtual void onstep() {}
    virtual void ondense(stat_t const&, idense_output&) {}
    virtual ~iequation_for_erk() {}
  };

  // DOP853: Dormand-Prince 8(5, 3)
  //
  //   The coefficients are taken from the supplementary material of the Heirer's book:
  //   from http://www.unige.ch/~hairer/prog/nonstiff/dop853.f
  //
  struct dop853_integrator {
    static const int stage = 12;
    static const int order = 8;
    mutable working_buffer buffer;

    // buffer の状態
    //   before [  ? | k1 |  ? |  ? |  ? |  ? |  ? |  ? |  ? |  ? ]
    //   after  [  x | k1 | k5 | kC | x6 | x7 | x8 | k9 | kA | kB ]
    //   但し x は次の step の値である。
    void _integrate8(
      double& time, double* ksh_restrict value, std::size_t size,
      iequation_for_erk& eq, double h,
      double atol, double rtol, double& _err, double& _stf
    ) const;

    struct param_t {
      double atol {1e-13};
      double rtol {1e-13};
      double beta {0.0};
      double fac1 {0.0};
      double fac2 {0.0};
      double safe {0.0};
      double hmax {0.0};
      double step {0.0};
      int    nstif {0};
      std::ptrdiff_t nmax {100000};
    };

    void _dense_output_initialize(
      working_buffer& interpBuffer, stat_t& stat, int* icomp, std::size_t nrd
    ) const;

    double _determine_initial_step(
      double time, double* ksh_restrict value, std::size_t size,
      iequation_for_erk& eq,
      int bwd, double atol, double rtol, double hmax
    ) const;

  private:
    template<typename F>
    struct eq_by_f:iequation_for_erk {
      F const& f;
      eq_by_f(F const& f):f(f) {}
      virtual void eval_f(double* ksh_restrict slope, double t, double const* ksh_restrict value) override {
        f(slope, t, value);
      }
    };

    template<typename F, typename CB>
    struct eq_by_f_and_cb:iequation_for_erk {
      F const& f;
      CB const& stepCallback;

      eq_by_f_and_cb(F const& f, CB const& stepCallback):f(f), stepCallback(stepCallback) {}
      virtual void eval_f(double* ksh_restrict slope, double t, double const* ksh_restrict value) override {
        f(slope, t, value);
      }
      virtual void onstep() {
        stepCallback();
      }
    };

  public:
    template<typename F>
    void operator()(double& time, double* ksh_restrict value, std::size_t size, F const& f, double h) const {
      eq_by_f<F> eq(f);

      buffer.ensure<double>(10 * size);

      double const atol = 1e-12;
      double const rtol = 1e-12;
      double err, stf;

      double* ksh_restrict  x  = buffer.ptr<double>();
      double* ksh_restrict  k1 = buffer.ptr<double>() + size * 1;
      eq.eval_f(k1, time, value);
      this->_integrate8(
        time, value, size, eq, h,
        atol, rtol, err, stf
      );

      for (std::size_t i = 0; i < size; i++)
        value[i] = x[i];
      time += h;
    }

    template<typename F>
    typename std::enable_if<!std::is_base_of<iequation_for_erk, F>::value, void>::type
    integrate(
      double& time, double* ksh_restrict value, std::size_t size, F const& f,
      double timeN, stat_t& stat, param_t const& params
    ) const {
      eq_by_f<F> eq(f);
      this->integrate(time, value, size, eq, timeN, stat, params);
    }

    template<typename F, typename CB>
    void integrate(
      double& time, double* ksh_restrict value, std::size_t size, F const& f,
      double timeN, stat_t& stat, param_t const& params, CB const& stepCallback
    ) const {
      eq_by_f_and_cb<F, CB> eq(f, stepCallback);
      this->integrate(time, value, size, eq, timeN, stat, params);
    }

  public:
    void integrate(
      double& time, double* ksh_restrict value, std::size_t size,
      iequation_for_erk& eq,
      double timeN, stat_t& stat, param_t const& params
    ) const;

    struct dense_output;
  };

}
}
#endif
