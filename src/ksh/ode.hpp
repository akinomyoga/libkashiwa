// -*- mode:c++ -*-
#ifndef kashiwa_ode_hpp
#define kashiwa_ode_hpp
#include <cstddef>
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

    void clear() {
      nfcn   = 0;
      nstep  = 0;
      naccpt = 0;
      nrejct = 0;
    }
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
    virtual ~iequation_for_erk() {}
    virtual void eval_derivative(double* ksh_restrict slope, std::size_t size, double t, double const* ksh_restrict value) = 0;
    virtual void onstep(double t, double const* ksh_restrict value) {}
    virtual void ondense(stat_t const&, idense_output&) {}
  public:
    virtual bool is_stopping() { return false; }
  };

  template<typename CRTP>
  class explicit_integrator_base {
    typedef CRTP integrator_type;

  public:
    void integrate_with_constant_step(
      double& time, double* ksh_restrict value, std::size_t size,
      iequation_for_erk& eq,
      double time_end, stat_t& stat, double time_step
    ) const {
      if (time == time_end) return;
      double const sgn = time < time_end ? 1 : -1;

      bool last_step = false;
      while (!last_step) {
        if (time_step > sgn * (time_end - time) * (1.0 - 1e-8)) {
          last_step = true;
          time_step = sgn * (time_end - time);
        }

        static_cast<integrator_type const&>(*this)(
          time, value, size, [&eq](double* ksh_restrict slope, std::size_t size, double t, double const* ksh_restrict value) {
            eq.eval_derivative(slope, size, t, value);
          }, sgn * time_step);
        stat.nfcn += integrator_type::stage;
        stat.nstep++;
        stat.naccpt++;
        eq.onstep(time, value);
      }
    }

    void integrate_with_equal_step(
      double& time, double* ksh_restrict value, std::size_t size,
      iequation_for_erk& eq,
      double time_end, stat_t& stat, std::size_t nstep
    ) const {
      for (std::size_t i = 0; i < nstep; i++) {
        double const time_step = (time_end - time) / (nstep - i);
        static_cast<integrator_type const&>(*this)(
          time, value, size, [&eq](double* ksh_restrict slope, std::size_t size, double t, double const* ksh_restrict value) {
            eq.eval_derivative(slope, size, t, value);
          }, time_step);
        stat.nfcn += integrator_type::stage;
        stat.nstep++;
        stat.naccpt++;
        eq.onstep(time, value);
      }
    }
  };

}
}
#endif
