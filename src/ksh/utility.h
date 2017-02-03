// -*- mode:c++ -*-
#ifndef KASHIWA_UTILITY_H
#define KASHIWA_UTILITY_H
namespace kashiwa {

  template<typename T>
  T const& clamp(T const& value, T const& lowerBound, T const& upperBound) {
    return value < lowerBound? lowerBound: value > upperBound? upperBound: value;
  }

}
#endif
