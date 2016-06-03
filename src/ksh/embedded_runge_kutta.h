// -*- mode:c++ -*-
#pragma once
#ifndef KASHIWA_RK_EMBEDDED_RUNGE_KUTTA_H
#define KASHIWA_RK_EMBEDDED_RUNGE_KUTTA_H
#include <cstddef>
#include <cfloat>
#include "buffer.h"

namespace kashiwa{
namespace rk16{

  struct iequation_for_erk{
    virtual void eval_f(double* __restrict__ slope,double t,double const* __restrict__ value) = 0;
    virtual void onstep(){}
    virtual ~iequation_for_erk(){}
  };

  // DOP853: Dormand-Prince 8(5,3)
  //   from http://www.unige.ch/~hairer/prog/nonstiff/dop853.f (ハイラーの本の附録)
  struct dop853_integrator{
    static const int stage = 12;
    static const int order = 8;
    mutable working_buffer buffer;

    // buffer の状態
    //   before [  ? | k1 |  ? |  ? |  ? |  ? |  ? |  ? |  ? |  ? ]
    //   after  [  x | k1 | k5 | kC | x6 | x7 | x8 | k9 | kA | kB ]
    void _integrate8(
      double& time,double* __restrict__ value,std::size_t size,
      iequation_for_erk& eq,double h,
      double atol,double rtol,double& _err,double& _stf
    ) const;

    struct stat_t{
      int nfcn   { 0 };
      int nstep  { 0 };
      int naccpt { 0 };
      int nrejct { 0 };
    };

    struct param_t{
      double atol {1e-13};
      double rtol {1e-13};
      double beta {0.0};
      double fac1 {0.0};
      double fac2 {0.0};
      double safe {0.0};
      double hmax {0.0};
      double step {0.0};
      int    nstif{0};
      std::ptrdiff_t nmax {100000};
    };

    // buffer の状態
    //   before [  x | k1 | kD | kC | k6 | k7 | k8 | k9 | kA | kB ]
    //   after  [  x | k1 | kD | kC | xE | xF | xG | kE | kF | kG ]
    void _dense_output_initialize(
      double const& time,double const* __restrict__ value,std::size_t size,
      iequation_for_erk& eq,double h,stat_t& stat,int* icomp,std::size_t nrd,double* __restrict__ cont
    ) const;

    double _determine_initial_step(
      double time,double* __restrict__ value,std::size_t size,
      iequation_for_erk& eq,
      int bwd,double atol,double rtol,double hmax
    ) const;

  private:
    template<typename F>
    struct eq_by_f:iequation_for_erk{
      F const& f;
      eq_by_f(F const& f):f(f){}
      virtual void eval_f(double* __restrict__ slope,double t,double const* __restrict__ value) override{
        f(slope,t,value);
      }
    };

    template<typename F,typename CB>
    struct eq_by_f_and_cb:iequation_for_erk{
      F const& f;
      CB const& stepCallback;

      eq_by_f_and_cb(F const& f,CB const& stepCallback):f(f),stepCallback(stepCallback){}
      virtual void eval_f(double* __restrict__ slope,double t,double const* __restrict__ value) override{
        f(slope,t,value);
      }
      virtual void onstep(){
        stepCallback();
      }
    };

  public:
    template<typename F>
    void operator()(double& time,double* __restrict__ value,std::size_t size,F const& f,double h) const{
      eq_by_f<F> eq(f);

      buffer.ensure(10*size);

      double const atol = 1e-12;
      double const rtol = 1e-12;
      double err,stf;

      double* __restrict__  x = buffer.ptr();
      double* __restrict__  k1 = buffer.ptr()+size*1;
      eq.eval_f(k1,time,value);
      this->_integrate8(
        time,value,size,eq,h,
        atol,rtol,err,stf
      );

      for(std::size_t i=0;i<size;i++)
        value[i] = x[i];
      time+=h;
    }

    template<typename F>
    typename std::enable_if<!std::is_base_of<iequation_for_erk,F>::value,void>::type
    integrate(
      double& time,double* __restrict__ value,std::size_t size,F const& f,
      double timeN,stat_t& stat,param_t const& params
    ) const{
      eq_by_f<F> eq(f);
      this->integrate(time,value,size,eq,timeN,stat,params);
    }

    template<typename F,typename CB>
    void integrate(
      double& time,double* __restrict__ value,std::size_t size,F const& f,
      double timeN,stat_t& stat,param_t const& params,CB const& stepCallback
    ) const{
      eq_by_f_and_cb<F,CB> eq(f,stepCallback);
      this->integrate(time,value,size,eq,timeN,stat,params);
    }

  private:
    std::size_t denseVersion {0};
    std::size_t previousSize;
    double      previousTime;
    double      previousStep;

  public:
    void integrate(
      double& time,double* __restrict__ value,std::size_t size,
      iequation_for_erk& eq,
      double timeN,stat_t& stat,param_t const& params
    ) const;

    // struct dense_type{
    //   dop853_integrator const& integ;
    //   std::vector<std::size_t> icomp;
    //   working_buffer buffer;

    //   std::size_t denseVersion {-1};
    //   double previousTime;
    //   double previousStep;
    // public:
    //   dense_type(dop853_integrator const& integ):integ(integ){}
    //   dense_type(dop853_integrator const& integ,std::size_t const* icomp,std::size_t ncomp)
    //     :integ(integ),icomp(icomp,icomp+ncomp){}

    //   void get_value_at(double* __restrict__ value,double time) const{
    //     int const ncomp = icomp.size()? icomp.size(): integ.presviousSize;

    //     // void _dense_output_initialize(
    //     //   double time,double const* __restrict__ value,std::size_t size,
    //     //   F const& f,double h,stat_t& stat,int* icomp,std::size_t nrd,double* __restrict__ cont
    //     // );

    //     integ._dense_output_initialize(time,value,size,f,h,stat,&icomp[0],ncomp);

    //     double const s = (time-integ.previousTime)/integ.previousStep;
    //     double const t = 1.0-s;
    //     mwg_check(0.0<=s&&s<=1.0, "time out of range.");
    //     double const con0 = con;
    //     double const con1 = con+ncomp;
    //     double const con2 = con+ncomp*2;
    //     double const con3 = con+ncomp*3;
    //     double const con4 = con+ncomp*4;

    //     for(std::size_t i=0; i<ncomp; i++){
    //       int const j = icomp.size()? icomp[i]: i;

    //       double a=con[i+7*ncomp];
    //       a=a*s+con[i+6*ncomp];
    //       a=a*t+con[i+5*ncomp];
    //       a=a*s+con[i+4*ncomp];
    //       a=a*t+con[i+3*ncomp];
    //       a=a*s+con[i+2*ncomp];
    //       a=a*t+con[i+1*ncomp];
    //       a=a*s+con[i+0*ncomp];
    //       value[i]=a;
    //     }
    //   }
    // };

  };

}
}
#endif
