// -*- mode:c++ -*-
#pragma once
#ifndef KASHIWA_RKTEST_H
#define KASHIWA_RKTEST_H
#include <cstdio>
#include <vector>
#include <cstddef>
#include <cmath>
#include <algorithm>
#include <utility>
#include <mwg/except.h>

template<typename F>
double binary_search_function(double lowerBound,double upperBound,double value,double tolerance,F func){
  double yl = func(lowerBound)-value;
  double yu = func(upperBound)-value;
  if(yl*yu>=0){
    return std::abs(yl)<=std::abs(yu)?lowerBound: upperBound;
  }else if(yl>0){
    using namespace std;
    swap(lowerBound, upperBound);
  }

  for(int i=0;i<54;i++){
    double const middle = 0.5*(lowerBound+upperBound);
    if(std::abs(lowerBound - upperBound)<=tolerance)
      return middle;

    double const ym = func(middle)-value;
    (ym<=0?lowerBound:upperBound)=middle;
  }

  return 0.5*(lowerBound+upperBound);
}

//-----------------------------------------------------------------------------
// 方程式
//
// dx/dt = 1/x
//   解 √(2(t-t0))
//   初期条件を X0 とすると、√(-2 t0) = X0
//   解 √(2t + X0^2)

// static const double initialCondition = 1.0;
// static const double finalTime = 1.0;
// static void f(double* __restrict__ slope,double t,double const* __restrict__ value){
//   slope[0] = 1.0/value[0];
// }
// static double exactSolution(double t){
//   return std::sqrt(2*t+initialCondition*initialCondition);
// }


// 方程式
//
// dx/dt = (tan(y) + 1)/2.
//   解: t - t0 = x + ln(sin(x)+cos(x)). (但し、t in [0, t*).)
//   初期条件を x(t=0) = 0 とすると A = 1.
//   特異点 t* = t(x=pi/2) = pi/2 = 1.5707963267948966
//

static const double initialCondition = 1.0;
static const double initialTime = 1.0 + std::log(std::sin(1.0)+std::cos(1.0)); // = 1.3233676675153825
static const double finalTime = initialTime + 0.2; // ~ 1.52
static void f(double* __restrict__ slope,double t,double const* __restrict__ value){
  slope[0] = (std::tan(value[0])+1.0)/2;
}
static double exactSolution(double t){
  return binary_search_function(
    0.0,0.5*M_PI,t,0.0,[](double x){
      return x+std::log(std::cos(x)+std::sin(x));
    }
  );
}

//-----------------------------------------------------------------------------
// test ERK

template<typename Integrator>
void test_method(std::FILE* file,Integrator const& integ){
  int const nstage = Integrator::stage;

  int    previousNStep = 0;
  double previousError = 0;
  bool   roundingDominant = false;
  double maxOrder = 0.0;
  std::vector<double> orders;

  double const sol = exactSolution(finalTime);

  std::size_t nvalMax =
    Integrator::order>=5?0x1000:
    Integrator::order>=4?0x10000:
    0x100000;
  for(std::size_t nval=1;nval<nvalMax;nval*=2){
    double time = initialTime;
    double value[1] = { initialCondition };

    if(nval*2 <= nstage) continue;

    int    const nstep = (nval+nstage-1)/nstage;
    double const h = (finalTime-time)/nstep;
    for(int i=0;i<nstep;i++)
      integ(time,value,1,f,h);

    double const err = sol-value[0];
    std::fprintf(file,"%u %g %g %g %g\n",nstep*nstage,time,value[0],sol,err);

    // check order
    if(!roundingDominant){
      if(previousNStep){
        double const order = -std::log(std::abs(err/previousError))/std::log((double)nstep/(double)previousNStep);
        if(order<0)
          roundingDominant = true;
        else{
          orders.push_back(order);
          if(order>maxOrder)
            maxOrder = order;
        }
      }
      previousNStep = nstep;
      previousError = err;
    }
  }

  // check order
  std::vector<double>::const_iterator b = orders.begin();
  std::vector<double>::const_iterator e = orders.end();
  while(e!=b&&e[-1]<maxOrder-1.0) e--;
  while(b!=e&&*b<maxOrder-1.0) b++;
  double const estimatedOrder = std::accumulate(b, e, 0.0)/(e - b);
  mwg_assert_nothrow(
    std::round(estimatedOrder)>=Integrator::order,
    "RK order mismatch: expected = %g, estimated = %g (%d-%d), max = %g",
    (double)Integrator::order,
    estimatedOrder, b - orders.begin(), e - orders.begin(),
    maxOrder);

}
#endif
