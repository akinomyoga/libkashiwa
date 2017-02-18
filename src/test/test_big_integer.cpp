#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <ksh/big_integer.h>
#include <mwg/except.h>

void test() {
  using bigint = kashiwa::bigint;
  mwg_check((bigint {1234} == 1234));
  mwg_check((bigint {1234} == bigint {1200} + bigint {34}));
}

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

// http://d.hatena.ne.jp/ku-ma-me/20080116/p1
void test3() {
  // Note: 基数10^nで計算すると遅い (2017-02-05)
  //auto result = pow(kashiwa::big_integer<std::uint32_t, std::uint64_t, 1000000000> {5}, 100000);

  int const n = 10000;
  // int const n = 30000;
  // int const n = 100000;
  // int const n = 300000;
  // int const n = 1000000;
  auto result = pow(kashiwa::bigint {5}, n);
  std::cout << "5 ** " << n << " = " << result << std::endl;
}

void test_div() {
  int const n = 100;
  auto a = pow(kashiwa::bigint {5}, n);
  auto b = a * a;
  mwg_check(b == pow(kashiwa::bigint {5}, 2 * n));
  mwg_check(b / a == a);
  mwg_check(b % a == 0);
  mwg_check(a / (uint32_t) 25 == pow(kashiwa::bigint {5}, n - 2));
  mwg_check(a % (uint32_t) 25 == 0);
}

int main() {
  test();
  // test1();
  // test2();
  test3();
  test_div();
  return 0;
}
