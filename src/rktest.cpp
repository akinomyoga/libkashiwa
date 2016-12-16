#include <cstdio>
#include <mwg/except.h>
#include "rktest.h"
#include "ksh/embedded_runge_kutta.h"

void test_erk1();

template<typename Integrator>
void test_embedded_integrate(std::FILE* file,Integrator const& integ){
  double const sol = exactSolution(finalTime);

  for(double tol=1e-15;tol<1e-4;tol*=10){
    kashiwa::runge_kutta::stat_t stat;
    typename Integrator::param_t params;
    params.atol = tol;
    params.rtol = tol;
    //params.step = 0.1;

    double time = initialTime;
    double value[1] = { initialCondition };
    integ.integrate(time,value,1,f,finalTime,stat,params);

    double const err = sol-value[0];
    std::fprintf(file,"%u %g %g %g %g\n",stat.nfcn,time,value[0],sol,err);
    std::printf("DOP853: nfe = %d, nac:nrj = %d:%d\n",stat.nfcn,stat.naccpt,stat.nrejct);
  }
}

void test_erk2(){
  std::FILE* file;

  mwg_printd("dop853");
  file = std::fopen("../out/rk/dop853.txt","wb");
  test_method(file,kashiwa::runge_kutta::dop853_integrator());
  std::fclose(file);

  mwg_printd("dop853");
  file = std::fopen("../out/rk/mdop853.txt","wb");
  test_embedded_integrate(file,kashiwa::runge_kutta::dop853_integrator());
  std::fclose(file);
}

int main(){
  test_erk1();
  test_erk2();
  return 0;
}
