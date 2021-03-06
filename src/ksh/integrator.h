// -*- C++ -*-
#ifndef KASHIWA_INTEGRATOR_H
#define KASHIWA_INTEGRATOR_H
#ifdef _MSC_VER
# define _USE_MATH_DEFINES
#endif
#include <cmath>
#include <algorithm>
namespace kashiwa{
//NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
namespace integrator_detail{
  struct point_weight_pair{
    double t;
    double w;
  };

  template<int I>
  struct GaussLegendre{
    static const int order=I;
    static const int data_size=(I+1)/2;
    static point_weight_pair data[data_size];
  };
}
//NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN

  template<typename F>
  double IntegrateByTrapezoid(const double lower,const double upper,const int iN,const F& f){
    double const dx=(upper-lower)/iN;

    double s=0.5*(f(lower)+f(upper));
    for(int i=1;i<iN;i++)
      s+=f(lower+i*dx);

    return s*dx;
  }

  template<int I,typename F>
  double IntegrateByGaussLegendre(const double lower,const double upper,const F& f){
    typedef integrator_detail::GaussLegendre<I> traits_t;
    double const center=0.5*(upper+lower);
    double const dxdt  =0.5*(upper-lower);

    double s=0;
    for(int i=0;i<I/2;i++){
      double const t=traits_t::data[i].t;
      double const w=traits_t::data[i].w;
      s+=w*(f(center+dxdt*t)+f(center-dxdt*t));
    }

    return s*dxdt;
  }

  template<int I>
  class GaussLegendreIntegrator{
    double const xmin;
    double const xmax;
    double const center;
    double const dxdt;
    integrator_detail:: point_weight_pair data[(I+1)/2];
  public:
    GaussLegendreIntegrator(double xmin,double xmax)
      :xmin(xmin),xmax(xmax),
       center(0.5*(xmax+xmin)),
       dxdt(0.5*(xmax-xmin))
    {
      typedef integrator_detail::GaussLegendre<I> traits_t;

      for(int i=0;i<(I+1)/2;i++){
        double const t=traits_t::data[i].t;
        double const w=traits_t::data[i].w;
        this->data[i].t=dxdt*t;
        this->data[i].w=dxdt*w;
      }
    }
    template<typename F>
    double Integrate(const F& f){
      double s=0;
      for(int i=0;i<I/2;i++){
        double const t=this->data[i].t;
        double const w=this->data[i].w;
        s+=w*(f(center+t)+f(center-t));
      }

      return s;
    }
  };

  template<int I, typename F>
  void gauss_legendre_quadrature(int N, double* result, double xmin, double xmax, F f) {
    typedef integrator_detail::GaussLegendre<I> traits_t;
    double const center = 0.5 * (xmax + xmin);
    double const dxdt   = 0.5 * (xmax - xmin);

    std::fill(result, result + N, 0.0);

    for (int i = 0; i < I / 2; i++) {
      double const t = traits_t::data[i].t;
      double const w = traits_t::data[i].w;

      double value1[N];
      f(value1, center + dxdt * t);
      double value2[N];
      f(value2, center - dxdt * t);

      for (int k = 0; k < N; k++)
        result[k] += w * (value1[k] + value2[k]);
    }

    for (int k = 0; k < N; k++)
      result[k] *= dxdt;
  }

  template<typename F>
  void gauss_chebyshev_quadrature(int Nd, double* result, int Ni, double xmin, double xmax, F f) {
    double const center = 0.5 * (xmax + xmin);
    double const dxdt   = 0.5 * (xmax - xmin);
    double const dh = M_PI / (2.0 * Ni);

    if (Ni & 1) {
      f(result, center);
    } else {
      std::fill(result, result + Nd, 0.0);
    }

    for (int i = 0, iN = Ni / 2; i < iN; i++) {
      double const t = std::cos((2 * i + 1) * dh);
      double const w = std::sqrt(1.0 - t * t);

      double value1[Nd];
      f(value1, center + dxdt * t);
      double value2[Nd];
      f(value2, center - dxdt * t);

      for (int k = 0; k < Nd; k++)
        result[k] += w * (value1[k] + value2[k]);
    }

    double const factor = dxdt * M_PI / Ni;
    for (int k = 0; k < Nd; k++)
      result[k] *= factor;
  }

//NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
}
#endif
