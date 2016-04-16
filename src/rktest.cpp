#include <cstdio>
#include <cmath>

//-----------------------------------------------------------------------------
// 方程式
//
// dx/dt = 1/x
//   解 √(2(t-t0))
//   初期条件を X0 とすると、√(-2 t0) = X0
//   解 √(2t + X0^2)
double f(double t,double value){
  return 1.0/value;
}
static const double initialCondition = 1.0;
double exactSolution(double t){
  return std::sqrt(2*t+initialCondition*initialCondition);
}

//-----------------------------------------------------------------------------

int main(){
  double value = initialCondition;
  double h = 0.001;
  double time = 0.0;

  for(int i=0;i<100;i++){
    // classical RK4
    double const k1 = h*f(time      ,value);
    double const k2 = h*f(time+0.5*h,value + 0.5*k1);
    double const k3 = h*f(time+0.5*h,value + 0.5*k2);
    double const k4 = h*f(time+h    ,value + k3);
    value+=(1.0/6.0)*(k1+2.0*(k2+k3)+k4);

    // // euler
    // double const k1 = h*f(time      ,value);
    // value+=k1;

    time+=h;

    double const sol=exactSolution(time);
    std::printf("%g %g %g %g\n",time,value,sol,sol-value);
  }

  return 0;
}
