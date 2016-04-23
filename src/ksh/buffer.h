// -*- mode:c++ -*-
#pragma once
#ifndef KASHIWA_BUFFER_H
#define KASHIWA_BUFFER_H
#include <cstddef>
#include <vector>

namespace kashiwa{
  struct working_buffer{
    std::vector<double> data;
  public:
    void ensure(std::size_t minimalSize){
      if(data.size()<minimalSize)
        data.resize(minimalSize,0.0);
    }

    double      * ptr()      {return &data[0];}
    double const* ptr() const{return &data[0];}
  };

  template<typename T>
  T const& clamp(T const& value,T const& lowerBound,T const& upperBound){
    return value<lowerBound?lowerBound: value>upperBound?upperBound: value;
  }
}
#endif
