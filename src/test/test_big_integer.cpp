#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <ksh/big_integer.h>
#include <mwg/except.h>

template<typename S, typename C, C M>
void dump(kashiwa::big_integer<S, C, M> const& value) {
  std::printf("sign = %d, data = [ ", value.sign);
  for (std::size_t i = value.data.size(); i--; ) {
    std::printf("%llx", (unsigned long long) value.data[i]);
    if (i != 0) std::printf(", ");
  }
  std::printf(" ]\n");
}

void test1() {
  kashiwa::big_integer<std::uint32_t, std::uint64_t, 0x100> a = 1234;
  dump(a);
  for (int i = 0; i <= 20; i++) {
    a += a;
    a++;
    dump(a);
  }
}

void test2() {
  //kashiwa::big_integer<std::uint32_t, std::uint64_t, 0x100> a = 1;
  kashiwa::big_integer<unsigned, unsigned, 10> a = 1;
  for (int n = 1; n <= 20; n++) {
    a = a * n;
    std::cout << n << "! = " << a << std::endl;
  }
}

int main() {
  test1();
  test2();
  return 0;
}
