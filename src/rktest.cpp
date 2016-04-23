#include <cstdio>
#include <cmath>
#include <vector>
#include <mwg/except.h>
#include <mwg/xprintf.h>
#include "ksh/erk.h"

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
// void f(double* __restrict__ slope,double t,double const* __restrict__ value){
//   slope[0] = 1.0/value[0];
// }
// double exactSolution(double t){
//   return std::sqrt(2*t+initialCondition*initialCondition);
// }


// 方程式
//
// dx/dt = (tan(y) + 1)/2.
//   解: t - t0 = x + ln(sin(x)+cos(x)). (但し、t in [0, t*).)
//   初期条件を x(t=0) = 0 とすると A = 1.
//   特異点 t* = pi/4 + (1/2)ln(2) = 1.131971753677421
//

static const double initialCondition = 0.0;
static const double finalTime = 1.13;
void f(double* __restrict__ slope,double t,double const* __restrict__ value){
  slope[0] = (std::tan(value[0])+1.0)/2;
}
double exactSolution(double t){
  return binary_search_function(
    0.0,0.25*M_PI,t,0.0,[](double x){
      return x+std::log(std::cos(x)+std::sin(x));
    }
  );
}


//-----------------------------------------------------------------------------
// RK integrator

//-----------------------------------------------------------------------------

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
    double time = 0.0;
    double value[1] = { initialCondition };

    if(nval*2 <= nstage) continue;

    int    const nstep = (nval+nstage-1)/nstage;
    double const h = (finalTime-time)/nstep;
    for(int i=0;i<nstep;i++)
      integ(time,value,1,f,h);

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
  mwg_assert_nothrow(
    std::round(estimatedOrder)>=Integrator::order,
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

  //
  // 2次公式
  //

  mwg_printd("mid");
  file = std::fopen("out/rk/rkmid.txt","wb");
  test_method(file,kashiwa::rk16::midpoint_integrator());
  std::fclose(file);

  mwg_printd("ralston");
  file = std::fopen("out/rk/rkral.txt","wb");
  test_method(file,kashiwa::rk16::ralston_integrator());
  std::fclose(file);

  mwg_printd("heun");
  file = std::fopen("out/rk/rkheun.txt","wb");
  test_method(file,kashiwa::rk16::heun_integrator());
  std::fclose(file);

  //
  // 3次公式
  //

  mwg_printd("runge3");
  file = std::fopen("out/rk/runge3.txt","wb");
  test_method(file,kashiwa::rk16::runge3_integrator());
  std::fclose(file);

  mwg_printd("heun3");
  file = std::fopen("out/rk/rkheun3.txt","wb");
  test_method(file,kashiwa::rk16::heun3_integrator());
  std::fclose(file);

  mwg_printd("ralston3");
  file = std::fopen("out/rk/rkral3.txt","wb");
  test_method(file,kashiwa::rk16::ralston3_integrator());
  std::fclose(file);

  mwg_printd("kutta3");
  file = std::fopen("out/rk/rkk3.txt","wb");
  test_method(file,kashiwa::rk16::kutta3_integrator());
  std::fclose(file);

  //
  // 4次公式
  //

  mwg_printd("gill");
  file = std::fopen("out/rk/rkgill.txt","wb");
  test_method(file,kashiwa::rk16::gill_integrator());
  std::fclose(file);

  mwg_printd("rk4");
  file = std::fopen("out/rk/rkrk4.txt","wb");
  test_method(file,kashiwa::rk16::rk4_integrator());
  std::fclose(file);

  mwg_printd("kutta38");
  file = std::fopen("out/rk/rkk38.txt","wb");
  test_method(file,kashiwa::rk16::kutta_3_8_integrator());
  std::fclose(file);

  //
  // 5次公式
  //

  mwg_printd("b5v1");
  file = std::fopen("out/rk/b5v1.txt","wb");
  test_method(file,kashiwa::rk16::butcher5v1_integrator());
  std::fclose(file);

  mwg_printd("b5v2");
  file = std::fopen("out/rk/b5v2.txt","wb");
  test_method(file,kashiwa::rk16::butcher5v2_integrator());
  std::fclose(file);

  mwg_printd("b5v3");
  file = std::fopen("out/rk/b5v3.txt","wb");
  test_method(file,kashiwa::rk16::butcher5v3_integrator());
  std::fclose(file);

  //
  // 高次公式
  //

  mwg_printd("hammud6");
  file = std::fopen("out/rk/hammud6.txt","wb");
  test_method(file,kashiwa::rk16::hammud6_integrator());
  std::fclose(file);

  mwg_printd("shanks7");
  file = std::fopen("out/rk/shanks7.txt","wb");
  test_method(file,kashiwa::rk16::shanks7_integrator());
  std::fclose(file);

  mwg_printd("cv7");
  file = std::fopen("out/rk/rkcv7.txt","wb");
  test_method(file,kashiwa::rk16::cooper_verner7_integrator());
  std::fclose(file);

  mwg_printd("cv8");
  file = std::fopen("out/rk/rkcv8.txt","wb");
  test_method(file,kashiwa::rk16::cooper_verner8_integrator());
  std::fclose(file);

  return 0;
}
