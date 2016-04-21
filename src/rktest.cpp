#include <cstdio>
#include <cmath>
#include <vector>
#include <mwg/except.h>
#include <mwg/xprintf.h>
#include "rk16.h"

//-----------------------------------------------------------------------------
// 方程式
//
// dx/dt = 1/x
//   解 √(2(t-t0))
//   初期条件を X0 とすると、√(-2 t0) = X0
//   解 √(2t + X0^2)

void f(double* __restrict__ slope,double t,double const* __restrict__ value){
  slope[0] = 1.0/value[0];
}

static const double initialCondition = 1.0;
double exactSolution(double t){
  return std::sqrt(2*t+initialCondition*initialCondition);
}

//-----------------------------------------------------------------------------
// RK integrator

//-----------------------------------------------------------------------------

template<typename Integrator>
void test_method(std::FILE* file,Integrator const& integ){
  double const finalTime = 1.0;
  int const nstage = Integrator::stage;

  int    previousNStep = 0;
  double previousError = 0;
  bool   roundingDominant = false;
  double maxOrder = 0.0;
  std::vector<double> orders;

  for(std::size_t nval=1;nval<0x100000;nval*=2){
    double time = 0.0;
    double value[1] = { initialCondition };

    if(nval*2 <= nstage) continue;

    int    const nstep = (nval+nstage-1)/nstage;
    double const h = (finalTime-time)/nstep;
    for(int i=0;i<nstep;i++)
      integ(time,value,1,f,h);

    double const sol = exactSolution(time);
    double const err = sol-value[0];
    std::fprintf(file,"%zu %g %g %g %g\n",nstep*nstage,time,value[0],sol,err);

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
  mwg_assert(
    std::round(estimatedOrder)==Integrator::order,
    "RK order mismatch: expected = %g, estimated = %g (%d-%d), max = %g\n",
    (double)Integrator::order,
    estimatedOrder, b - orders.begin(), e - orders.begin(),
    maxOrder);

}

int main(){
  std::FILE* file;

  mwg_printd("euler");
  file = std::fopen("out/rk/rkeuler.txt","wb");
  test_method(file,kashiwa::rk16::euler_integrator());
  std::fclose(file);

  mwg_printd("mid");
  file = std::fopen("out/rk/rkmid.txt","wb");
  test_method(file,kashiwa::rk16::midpoint_integrator());
  std::fclose(file);

  mwg_printd("rk4");
  file = std::fopen("out/rk/rkrk4.txt","wb");
  test_method(file,kashiwa::rk16::rk4_integrator());
  std::fclose(file);

  mwg_printd("cv8");
  file = std::fopen("out/rk/rkcv8.txt","wb");
  test_method(file,kashiwa::rk16::cooper_verner8_integrator());
  std::fclose(file);

  return 0;
}
