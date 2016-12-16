// -*- mode:c++ -*-
#ifndef KASHIWA_DEF_H
#define KASHIWA_DEF_H

namespace kashiwa {
  typedef unsigned char byte;

  template<typename T>
  T const& clamp(T const& value, T const& lowerBound, T const& upperBound) {
    return value < lowerBound? lowerBound: value > upperBound? upperBound: value;
  }
}

#endif
