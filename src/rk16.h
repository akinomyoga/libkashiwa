// -*- mode:c++ -*-
#pragma once
#ifndef KASHIWA_RK16_H
#define KASHIWA_RK16_H
#include <cstdlib>
#include <cmath>
#include <vector>

namespace kashiwa{
namespace rk16{

  struct working_buffer{
    std::vector<double> data;
  public:
    void ensure(std::size_t minimalSize){
      if(data.size()<minimalSize)
        data.resize(minimalSize);
    }

    double      * ptr()      {return &data[0];}
    double const* ptr() const{return &data[0];}
  };


  struct euler_integrator{
    static const int stage = 1;
    static const int order = 1;

    mutable working_buffer buffer;

    template<typename F>
    void operator()(double& time,double* value,std::size_t size,F const& f,double h) const{
      buffer.ensure(size);
      double* __restrict__ knode = buffer.ptr();
      f(knode,time,value);
      for(std::size_t i=0;i<size;i++)
        value[i] += h*knode[i];

      time+=h;
    }
  };

  struct midpoint_integrator{
    static const int stage = 2;
    static const int order = 2;
    mutable working_buffer buffer;

    template<typename F>
    void operator()(double& time,double* value,std::size_t size,F const& f,double h) const{
      buffer.ensure(2*size);
      double* __restrict__ knode = buffer.ptr();
      double* __restrict__ xnode = buffer.ptr()+size;

      f(knode,time,value);
      for(std::size_t i=0;i<size;i++)
        xnode[i] = value[i] + 0.5*h*knode[i];

      f(knode,time+0.5*h,xnode);
      for(std::size_t i=0;i<size;i++)
        value[i] += h*knode[i];

      time+=h;
    }
  };

  struct heun_integrator{
    static const int stage = 2;
    static const int order = 2;
    mutable working_buffer buffer;

    template<typename F>
    void operator()(double& time,double* value,std::size_t size,F const& f,double h) const{
      buffer.ensure(2*size);
      double* __restrict__ k = buffer.ptr();
      double* __restrict__ x = buffer.ptr()+size;

      f(k,time,value);
      for(std::size_t i=0;i<size;i++){
        x[i] = value[i] + h*k[i];
        value[i] += (1.0/2.0)*h*k[i];
      }

      f(k,time+0.5*h,x);
      for(std::size_t i=0;i<size;i++)
        value[i] += (1.0/2.0)*h*k[i];

      time+=h;
    }
  };

  struct ralston_integrator{
    static const int stage = 2;
    static const int order = 2;
    mutable working_buffer buffer;

    template<typename F>
    void operator()(double& time,double* value,std::size_t size,F const& f,double h) const{
      buffer.ensure(2*size);
      double* __restrict__ k = buffer.ptr();
      double* __restrict__ x = buffer.ptr()+size;

      f(k,time,value);
      for(std::size_t i=0;i<size;i++){
        x[i] = value[i] + (2.0/3.0)*h*k[i];
        value[i] += (1.0/4.0)*h*k[i];
      }

      f(k,time+0.5*h,x);
      for(std::size_t i=0;i<size;i++)
        value[i] += (3.0/4.0)*h*k[i];

      time+=h;
    }
  };

  // RK4 (classical Runge-Kutta method)

  struct rk4_integrator{
    static const int stage = 4;
    static const int order = 4;

    mutable working_buffer buffer;

    template<typename F>
    void operator()(double& time,double* value,std::size_t size,F const& f,double h) const{
      buffer.ensure(3*size);
      double* __restrict__ knode = buffer.ptr();
      double* __restrict__ xnode = buffer.ptr()+size;
      double* __restrict__ delta = buffer.ptr()+size*2;

      // k1
      f(knode,time,value);
      for(std::size_t i=0;i<size;i++){
        delta[i] = (1.0/6.0)*h*knode[i];
        xnode[i] = value[i] + 0.5*h*knode[i];
      }

      // k2
      f(knode,time+0.5*h,xnode);
      for(std::size_t i=0;i<size;i++){
        delta[i] += (2.0/6.0)*h*knode[i];
        xnode[i] = value[i] + 0.5*h*knode[i];
      }

      // k3
      f(knode,time+0.5*h,xnode);
      for(std::size_t i=0;i<size;i++){
        delta[i] += (2.0/6.0)*h*knode[i];
        xnode[i] = value[i] + h*knode[i];
      }

      // k4
      f(knode,time+h,xnode);
      for(std::size_t i=0;i<size;i++){
        double const a = delta[i] + (1.0/6.0)*h*knode[i];
        value[i] += a;
      }

      time+=h;
    }
  };

  // Kutta 3/8-rule
  struct kutta_3_8_integrator{
    static const int stage = 4;
    static const int order = 4;
    mutable working_buffer buffer;

    template<typename F>
    void operator()(double& time,double* value,std::size_t size,F const& f,double h) const{
      buffer.ensure(3*size);
      double* __restrict__ k  = buffer.ptr();
      double* __restrict__ xi = buffer.ptr()+size;
      double* __restrict__ x4 = buffer.ptr()+size*2;

      // k1
      f(k,time,value);
      for(std::size_t i=0;i<size;i++){
        xi[i] = value[i] + (1.0/3.0)*h*k[i];
        x4[i] = value[i] + h*k[i];
        value[i] += (1.0/8.0)*h*k[i];
      }

      // k2
      f(k,time+0.5*h,xi);
      for(std::size_t i=0;i<size;i++){
        xi[i] = 2.0*xi[i] - x4[i] + h*k[i];
        x4[i] -= h*k[i];
        value[i] += (3.0/8.0)*h*k[i];
      }

      // k3
      f(k,time+0.5*h,xi);
      for(std::size_t i=0;i<size;i++){
        x4[i] += h*k[i];
        value[i] += (3.0/8.0)*h*k[i];
      }

      // k4
      f(k,time+h,x4);
      for(std::size_t i=0;i<size;i++)
        value[i] += (1.0/8.0)*h*k[i];

      time+=h;
    }
  };

  // Cooper-Verner 8th order
  //   [cv8.1] E.ハイラー, 三井 斌友, 『常微分方程式の数値解法 I 』, 丸善出版 (2012/7/17).
  //   [cv8.2] http://www.330k.info/essay/Explicit-Runge-Kutta-Butcher-Tableau
  struct cooper_verner8_integrator{
    static const int stage = 11;
    static const int order = 8;
    mutable working_buffer buffer;

    static constexpr double sqrt21 = std::sqrt(21.0); // Ref [cv8.2] では -sqrt(21.0). どちらでも OK.

    template<typename F>
    void operator()(double& time,double* __restrict__ value,std::size_t size,F const& f,double h) const{
      buffer.ensure(8*size);
      double* __restrict__  delta = buffer.ptr();
      double* __restrict__  xnode = buffer.ptr()+size;
      double* __restrict__  k1    = buffer.ptr()+size*2;
      double* __restrict__  k2    = buffer.ptr()+size*3;
      double* __restrict__  k3    = buffer.ptr()+size*4;
      double* __restrict__& k4    = k2;
      double* __restrict__  k5    = buffer.ptr()+size*5;
      double* __restrict__  k6    = buffer.ptr()+size*6;
      double* __restrict__& k7    = k3;
      double* __restrict__& k8    = k4;
      double* __restrict__  k9    = buffer.ptr()+size*7;
      double* __restrict__& kA    = k1;
      double* __restrict__& kB    = k5;

      static constexpr double b1 =  1.0/ 20.0;
      static constexpr double b8 = 49.0/180.0;
      static constexpr double b9 = 16.0/ 45.0;
      static constexpr double bA = 49.0/180.0;
      static constexpr double bB =  1.0/ 20.0;

      // k1
      f(k1,time,value);

      // k2
      static constexpr double a21 = 0.5;
      static constexpr double c20 = a21;
      for(std::size_t i=0;i<size;i++){
        delta[i] = b1*h*k1[i];
        xnode[i] = value[i] + a21*h*k1[i];
      }
      f(k2,time+c20*h,xnode);

      // k3
      static constexpr double a31 = 0.25;
      static constexpr double a32 = 0.25;
      static constexpr double c30 = a31 + a32;
      for(std::size_t i=0;i<size;i++)
        xnode[i] = value[i] + a31*h*k1[i] + a32*h*k2[i];
      f(k3,time+c30*h,xnode);

      // k4 <= k2
      static constexpr double a41 = (1.0/ 7.0);
      static constexpr double a42 = (1.0/98.0)*(-7-3*sqrt21);
      static constexpr double a43 = (1.0/49.0)*(21+5*sqrt21);
      static constexpr double c40 = a41 + a42 + a43;
      for(std::size_t i=0;i<size;i++)
        xnode[i] = value[i] + a41*h*k1[i] + a42*h*k2[i] + a43*h*k3[i];
      f(k4,time+c40*h,xnode);

      // k5
      static constexpr double a51 = (1.0/ 84.0)*(11+1*sqrt21);
      static constexpr double a53 = (1.0/ 63.0)*(18+4*sqrt21);
      static constexpr double a54 = (1.0/252.0)*(21-1*sqrt21);
      static constexpr double c50 = a51 + a53 + a54;
      for(std::size_t i=0;i<size;i++)
        xnode[i] = value[i] + a51*h*k1[i] + a53*h*k3[i] + a54*h*k4[i];
      f(k5,time+c50*h,xnode);

      // k6
      static constexpr double a61 = (1.0/ 48.0)*(   5 +1*sqrt21);
      static constexpr double a63 = (1.0/ 36.0)*(   9 +1*sqrt21);
      static constexpr double a64 = (1.0/360.0)*(-231+14*sqrt21);
      static constexpr double a65 = (1.0/ 80.0)*(  63 -7*sqrt21);
      static constexpr double c60 = a61 + a63 + a64 + a65;
      for(std::size_t i=0;i<size;i++)
        xnode[i] = value[i] + a61*h*k1[i] + a63*h*k3[i] + a64*h*k4[i] + a65*h*k5[i];
      f(k6,time+c60*h,xnode);

      // k7 <= k3
      static constexpr double a71 = (1.0/ 42.0)*(  10-  1*sqrt21);
      static constexpr double a73 = (1.0/315.0)*(-432+ 92*sqrt21);
      static constexpr double a74 = (1.0/ 90.0)*( 633-145*sqrt21);
      static constexpr double a75 = (1.0/ 70.0)*(-504+115*sqrt21);
      static constexpr double a76 = (1.0/ 35.0)*(  63- 13*sqrt21);
      static constexpr double c70 = a71 + a73 + a74 + a75 + a76;
      for(std::size_t i=0;i<size;i++)
        xnode[i] = value[i] + a71*h*k1[i] + a73*h*k3[i] + a74*h*k4[i] + a75*h*k5[i] + a76*h*k6[i];
      f(k7,time+c70*h,xnode);

      // k8 <= k4
      static constexpr double a81 =  1.0/ 14.0;
      static constexpr double a85 = (1.0/126.0)*(14-3*sqrt21);
      static constexpr double a86 = (1.0/ 63.0)*(13-3*sqrt21);
      static constexpr double a87 =  1.0/  9.0;
      static constexpr double c80 = a81 + a85 + a86 + a87;
      for(std::size_t i=0;i<size;i++)
        xnode[i] = value[i] + a81*h*k1[i] + a85*h*k5[i] + a86*h*k6[i] + a87*h*k7[i];
      f(k8,time+c80*h,xnode);

      // k9
      static constexpr double a91 =   1.0/  32.0;
      static constexpr double a95 = ( 1.0/ 576.0)*(  91.0-21.0*sqrt21);
      static constexpr double a96 =  11.0/  72.0;
      static constexpr double a97 = ( 1.0/1152.0)*(-385.0-75.0*sqrt21);
      static constexpr double a98 = ( 1.0/ 128.0)*(  63.0+13.0*sqrt21);
      static constexpr double c90 = a91 + a95 + a96 + a97 + a98;
      for(std::size_t i=0;i<size;i++){
        delta[i] += b8*h*k8[i];
        xnode[i] = value[i] + a91*h*k1[i] + a95*h*k5[i] + a96*h*k6[i] + a97*h*k7[i] + a98*h*k8[i];
      }
      f(k9,time+c90*h,xnode);

      // kA <= k1
      static constexpr double aA1 =  1.0/  14.0;
      static constexpr double aA5 =  1.0/   9.0;
      static constexpr double aA6 = (1.0/2205.0)*(-733.0-147.0*sqrt21);
      static constexpr double aA7 = (1.0/ 504.0)*( 515.0+111.0*sqrt21);
      static constexpr double aA8 = (1.0/  56.0)*(- 51.0- 11.0*sqrt21);
      static constexpr double aA9 = (1.0/ 245.0)*( 132.0+ 28.0*sqrt21);
      static constexpr double cA0 = aA1 + aA5 + aA6 + aA7 + aA8 + aA9;
      for(std::size_t i=0;i<size;i++){
        delta[i] += b9*h*k9[i];
        xnode[i] = value[i] + aA1*h*k1[i] + aA5*h*k5[i] + aA6*h*k6[i] + aA7*h*k7[i] + aA8*h*k8[i] + aA9*h*k9[i];
      }
      f(kA,time+cA0*h,xnode);

      // kB <= k5
      static constexpr double aB5 = (1.0/18.0)*(- 42.0+ 7.0*sqrt21);
      static constexpr double aB6 = (1.0/45.0)*(- 18.0+28.0*sqrt21);
      static constexpr double aB7 = (1.0/72.0)*(-273.0-53.0*sqrt21);
      static constexpr double aB8 = (1.0/72.0)*( 301.0+53.0*sqrt21);
      static constexpr double aB9 = (1.0/45.0)*(  28.0-28.0*sqrt21);
      static constexpr double aBA = (1.0/18.0)*(  49.0- 7.0*sqrt21);
      static constexpr double cB0 = aB5 + aB6 + aB7 + aB8 + aB9 + aBA;
      for(std::size_t i=0;i<size;i++){
        delta[i] += bA*h*kA[i];
        xnode[i] = value[i] + aB5*h*k5[i] + aB6*h*k6[i] + aB7*h*k7[i] + aB8*h*k8[i] + aB9*h*k9[i] + aBA*h*kA[i];
      }
      f(kB,time+cB0*h,xnode);

      // increment
      for(std::size_t i=0;i<size;i++){
        double const a = delta[i] + bB*h*kB[i];
        value[i] += a;
      }

      time+=h;
    }
  };

}
}

#endif
