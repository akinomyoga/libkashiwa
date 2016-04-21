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

  for(std::size_t nval=1;nval<0x100000;nval*=2){
    double time = 0.0;
    double value[1] = { initialCondition };

    if(nval*2 <= nstage) continue;

    int    const nstep = (nval+nstage-1)/nstage;
    double const h = (finalTime-time)/nstep;
    for(int i=0;i<nstep;i++)
      integ(time,value,1,f,h);

    double const sol=exactSolution(time);
    std::fprintf(file,"%zu %g %g %g %g\n",nstep*nstage,time,value[0],sol,sol-value[0]);
  }
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
