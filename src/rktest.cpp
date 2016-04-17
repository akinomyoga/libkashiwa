#include <cstdio>
#include <cmath>
#include <vector>
#include <mwg/except.h>
#include <mwg/xprintf.h>

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

// RK4 (classical Runge-Kutta method)

struct rk4_integrator{
  mutable working_buffer buffer;

  template<typename F>
  void operator()(double& time,double* value,std::size_t size,F const& f,double h) const{
    buffer.ensure(3*size);
    double* knode = buffer.ptr();
    double* xnode = buffer.ptr()+size;
    double* delta = buffer.ptr()+size*2;

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

struct euler_integrator{
  mutable working_buffer buffer;

  template<typename F>
  void operator()(double& time,double* value,std::size_t size,F const& f,double h) const{
    buffer.ensure(size);
    double* knode = buffer.ptr();
    f(knode,time,value);
    for(std::size_t i=0;i<size;i++)
      value[i] += h*knode[i];

    time+=h;
  }
};

}
}

//-----------------------------------------------------------------------------

template<typename Integrator>
void test_method(std::FILE* file,Integrator const& integ,int nstage){
  double const finalTime = 1.0;

  for(std::size_t nval=1;nval<0x10000;nval*=2){
    double time = 0.0;
    double value[1] = { initialCondition };

    int    const nstep = (nval+nstage-1)/nstage;
    double const h = (finalTime-time)/nstep;
    for(int i=0;i<nstep;i++)
      integ(time,value,1,f,h);

    double const sol=exactSolution(time);
    std::fprintf(file,"%zu %g %g %g %g\n",nval,time,value[0],sol,sol-value[0]);
  }
}

int main(){
  std::FILE* file;

  mwg_printd("rk4");
  file = std::fopen("out/rk/rkrk4.txt","wb");
  test_method(file,kashiwa::rk16::rk4_integrator(), 4);
  std::fclose(file);

  mwg_printd("euler");
  file = std::fopen("out/rk/rkeuler.txt","wb");
  test_method(file,kashiwa::rk16::euler_integrator(), 1);
  std::fclose(file);

  return 0;
}
