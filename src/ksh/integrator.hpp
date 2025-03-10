// -*- mode: c++ -*-
#ifndef kashiwa_integrator_hpp
#define kashiwa_integrator_hpp
#include <cstddef>
#include <cmath>
#include <complex>
#include <algorithm>
namespace kashiwa {
namespace integrator_detail {
  struct point_weight_pair {
    double t;
    double w;
  };

  template<std::size_t I>
  struct GaussLegendre {
    static const std::size_t order = I;
    static const std::size_t data_size = (I + 1) / 2;
    static point_weight_pair data[data_size];
  };

  template<std::size_t I>
  struct GaussLaguerre {
    static const std::size_t order = I;
    static const std::size_t data_size = I;
    static point_weight_pair data[I];
  };
}

  template<typename F>
  double IntegrateByTrapezoid(const double lower, const double upper, const int iN, const F& f) {
    double const dx = (upper - lower) / iN;

    double s = 0.5 * (f(lower) + f(upper));
    for (int i = 1; i < iN; i++)
      s += f(lower + i * dx);

    return s * dx;
  }

  template<std::size_t I, typename F>
  double IntegrateByGaussLegendre(const double lower, const double upper, const F& f) {
    typedef integrator_detail::GaussLegendre<I> traits_t;
    double const center = 0.5 * (upper + lower);
    double const dxdt   = 0.5 * (upper - lower);

    double s = 0.0;
    for (std::size_t i = 0; i < I / 2; i++) {
      double const t = traits_t::data[i].t;
      double const w = traits_t::data[i].w;
      s += w * (f(center + dxdt * t) + f(center - dxdt * t));
    }

    return s * dxdt;
  }

  template<std::size_t I, typename F>
  double IntegrateByGaussLaguerre(const double lower, const double scale, F f) {
    typedef integrator_detail::GaussLaguerre<I> traits_t;
    double s = 0.0;
    for (std::size_t i = 0; i < traits_t::data_size; i++) {
      double const t = traits_t::data[i].t;
      double const w = traits_t::data[i].w;
      s += w * f(lower + scale * t);
    }
    return s * scale;
  }

  template<std::size_t I>
  class GaussLegendreIntegrator {
    double const xmin;
    double const xmax;
    double const center;
    double const dxdt;
    integrator_detail:: point_weight_pair data[(I + 1) / 2];
  public:
    GaussLegendreIntegrator(double xmin,double xmax):
      xmin(xmin), xmax(xmax),
      center(0.5 * (xmax + xmin)),
      dxdt(0.5 * (xmax - xmin))
    {
      typedef integrator_detail::GaussLegendre<I> traits_t;

      for (int i = 0; i < (I + 1) / 2; i++) {
        double const t = traits_t::data[i].t;
        double const w = traits_t::data[i].w;
        this->data[i].t = dxdt * t;
        this->data[i].w = dxdt * w;
      }
    }
    template<typename F>
    double Integrate(const F& f) {
      double s = 0;
      for (int i = 0; i < I / 2; i++) {
        double const t = this->data[i].t;
        double const w = this->data[i].w;
        s += w * (f(center + t) + f(center - t));
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

  template<int I, typename F>
  void gauss_legendre_quadrature(int N, std::complex<double>* result, double xmin, double xmax, F f) {
    gauss_legendre_quadrature<I>(N * 2, reinterpret_cast<double*>(result), xmin, xmax, [&] (double* integrand, double x) {
      f(reinterpret_cast<std::complex<double>*>(integrand), x);
    });
  }
  template<int I, typename F>
  void gauss_chebyshev_quadrature(int N, std::complex<double>* result, double xmin, double xmax, F f) {
    gauss_chebyshev_quadrature<I>(N * 2, reinterpret_cast<double*>(result), xmin, xmax, [&] (double* integrand, double x) {
      f(reinterpret_cast<std::complex<double>*>(integrand), x);
    });
  }

}
#endif
