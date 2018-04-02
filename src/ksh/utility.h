// -*- mode:c++ -*-
#ifndef KASHIWA_UTILITY_H
#define KASHIWA_UTILITY_H
namespace kashiwa {

  template<typename T>
  T const& clamp(T const& value, T const& lowerBound, T const& upperBound) {
    return value < lowerBound? lowerBound: value > upperBound? upperBound: value;
  }

  template<typename F>
  double binary_search_function(double lowerBound, double upperBound, double value, double tolerance, F func) {
    double yl = func(lowerBound) - value;
    double yu = func(upperBound) - value;
    if (yl * yu >= 0) {
      return std::abs(yl) <= std::abs(yu)? lowerBound: upperBound;
    } else if (yl > 0) {
      using namespace std;
      swap(lowerBound, upperBound);
    }

    for (int i = 0; i < 54; i++) {
      double const middle = 0.5 * (lowerBound + upperBound);
      if (std::abs(lowerBound - upperBound) <= tolerance)
        return middle;

      double const ym = func(middle) - value;
      (ym <= 0? lowerBound: upperBound) = middle;
    }

    return 0.5 * (lowerBound + upperBound);
  }

}
#endif
