// -*- mode:c++ -*-
#ifndef KASHIWA_BUFFER_H
#define KASHIWA_BUFFER_H
#include <cstddef>
#include <vector>

namespace kashiwa {

  struct working_buffer {
    std::vector<char> data;

  public:
    template<typename T = char>
    void ensure(std::size_t minimalSize) {
      if(data.size() < sizeof(T) * minimalSize)
        data.resize(sizeof(T) * minimalSize, 0.0);
    }

    template<typename T = char> T      * ptr()       {return reinterpret_cast<T*>(&data[0]);}
    template<typename T = char> T const* ptr() const {return reinterpret_cast<T const*>(&data[0]);}

    std::size_t size() const {return this->data.size();}
  };

}
#endif
