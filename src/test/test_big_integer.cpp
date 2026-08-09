#include <cinttypes>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <vector>
#include <ksh/big_integer.hpp>
#include <mwg/except.h>

namespace early_tests {
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
    using int_t = kashiwa::bigint;

    int const n = 100;

    auto a = pow(int_t{5}, n);
    mwg_check(a / (uint32_t) 25 == pow(int_t {5}, n - 2));
    mwg_check(a % (uint32_t) 25 == 0);

    auto b = a * a;
    mwg_check(b == pow(int_t {5}, 2 * n));
    mwg_check(b / a == a);
    mwg_check(b % a == 0);

    auto uint32_mod = pow(int_t{2}, 32);
    auto uint64_mod = pow(int_t{2}, 64);

    // 2017-05-18 bug: 繰り下がりで桁数が減少した時 resize していなかった
    mwg_check(uint64_mod - 1 + 1 == uint64_mod);

    // 2017-05-18 bug: arhs (商の概算用の法) の overflow があった
    mwg_check(int_t{39118700544} % int_t{37497702955} == 1620997589);

    // 2017-05-18 bug: arhs で除ききれない remainder が残っていた
    mwg_check(2 * uint64_mod / (uint64_mod - 1) == 2);
    mwg_check(a / a == 1);
    mwg_check(a % a == 0);
  }

  void run_tests() {
    test();
    test1();
    test2();
    test3();
    test_div();
  }
}

//-----------------------------------------------------------------------------

namespace kashiwa {
  std::vector<std::int32_t> primes_up_to(std::int32_t B);
}

namespace factorize_tests {

  void make_table() {
    std::vector<std::uint64_t> mod255_table((std::size_t) 4, (std::uint64_t) 0);
    std::vector<std::uint64_t> mod256_table((std::size_t) 4, (std::uint64_t) 0);
    std::vector<std::uint64_t> mod257_table((std::size_t) 5, (std::uint64_t) 0);
    for (std::uint32_t i = 0; i < 257; i++) {
      std::uint32_t const mod255 = i * i % 255;
      std::uint32_t const mod256 = i * i % 256;
      std::uint32_t const mod257 = i * i % 257;
      mod255_table[mod255 / 64] |= (std::uint64_t) 1 << (mod255 % 64);
      mod256_table[mod256 / 64] |= (std::uint64_t) 1 << (mod256 % 64);
      mod257_table[mod257 / 64] |= (std::uint64_t) 1 << (mod257 % 64);
    }
    std::printf("static std::uint64_t mod255_exclude_table[] = {");
    for (std::uint64_t x: mod255_table) std::printf("%#018" PRIx64 ", ", ~x);
    std::printf("};\n");
    std::printf("static std::uint64_t mod256_exclude_table[] = {");
    for (std::uint64_t x: mod256_table) std::printf("%#018" PRIx64 ", ", ~x);
    std::printf("};\n");
    std::printf("static std::uint64_t mod257_exclude_table[] = {");
    for (std::uint64_t x: mod257_table) std::printf("%#018" PRIx64 ", ", ~x);
    std::printf("};\n");
  }

  void check_primes_up_to() {
    // for (auto a: primes_up_to(1000000)) std::cout << a << "\n";

    constexpr int count_expect = 78498;
    constexpr std::uint64_t sum_expect = 37550402023;
    int count = 0;
    std::uint64_t sum = 0;
    for (auto a: kashiwa::primes_up_to(1000000)) {
      count++;
      sum += a;
    }
    std::cout << "count = " << count << " (expect: " << count_expect << ")" << std::endl;
    std::cout << "sum = " << sum << " (expect: " << sum_expect << ")" << std::endl;
    mwg_check(count == count_expect);
    mwg_check(sum == sum_expect);
  }
}

//-----------------------------------------------------------------------------

int main() {
  early_tests::run_tests();
  factorize_tests::make_table();
  factorize_tests::check_primes_up_to();
  return 0;
}
