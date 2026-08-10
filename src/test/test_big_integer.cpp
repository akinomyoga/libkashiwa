#include <cinttypes>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <random>
#include <sstream>
#include <string>
#include <vector>
#include <mwg/except.h>
#include <ksh/big_integer.hpp>

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

  void make_table_for_squares() {
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

namespace basic_tests {

  typedef kashiwa::big_integer<> int_t;
  const std::string literal_prefix = "0x";
  constexpr std::size_t literal_radix = 16;
  constexpr std::size_t elem_digits = 8;

  char digit(int number) {
    if (0 <= number && number <= 9) {
      return number + '0';
    } else if (10 <= number && number <= 25) {
      return "abcdefghijklmnopqrstuvwxyz"[number - 10];
    } else {
      mwg_check(0, "BUG");
    }
  }

  std::string make_literal_random(std::size_t digits, unsigned radix) {
    if (digits == 0) return "0";

    static std::mt19937 rng(1337);

    std::uniform_int_distribution<int> urandf(1, radix - 1);
    std::uniform_int_distribution<int> urand(0, radix - 1);

    int_t x = (int_t::element_type) urandf(rng);
    for (std::size_t i = 1; i < digits; ++i)
      x = x * radix + urand(rng);
    std::ostringstream s;
    s << x;
    return s.str();
  }

  // 9999...9999
  std::string make_literal_nines(std::size_t digits, unsigned radix) {
    std::ostringstream s;
    s << kashiwa::ipow((int_t) radix, digits) - 1;
    return s.str();
  }

  // 10000...0000
  std::string make_literal_one_zeros(std::size_t digits, unsigned radix) {
    if (digits == 0) return "0";
    std::ostringstream s;
    s << kashiwa::ipow((int_t) radix, digits - 1);
    return s.str();
  }

  bool arith_save_cases(const char* filename) {
    std::ofstream out(filename);
    if (!out) {
      std::cerr << "ファイルを開けませんでした。" << std::endl;
      return false;
    }

    auto output_with_signs = [&out] (std::string const& a, const char* op, std::string const& b) {
      out << a << " " << op << " " << b << "\n";
      if (b != "0")
        out << a << " " << op << " -" << b << "\n";
      if (a != "0")
        out << "-" << a << " " << op << " " << b << "\n";
      if (a != "0" && b != "0")
        out << "-" << a << " " << op << " -" << b << "\n";
    };
    auto output_with_signs_and_opts = [&out, &output_with_signs] (std::string const& a, std::string const& b) {
      static const char* op_other[3] = {"+", "-", "*"};
      static const char* op_div[2] = {"/", "%"};
      for (const char* op: op_other) {
        output_with_signs(a, op, b);
        if (b != a)
          output_with_signs(b, op, a);
      }
      for (const char* op: op_div) {
        if (b != "0")
          output_with_signs(a, op, b);
        if (b != a && a != "0")
          output_with_signs(b, op, a);
      }
    };

    std::vector<std::pair<std::size_t, std::size_t>> size_pairs = {
      {100, 100}, {100, 99}, {100, 2}, {1000, 1}, {50, 80},
    };
    for (const auto& pair : size_pairs) {
      for (std::size_t i = 0; i < 100; i++) {
        output_with_signs_and_opts(
          make_literal_random(pair.first, literal_radix),
          make_literal_random(pair.second, literal_radix));
      }
    }

    std::vector<std::size_t> test_digits = {5, 10, 50, 100};
    for (std::size_t digits: test_digits) {
      for (std::size_t i = 0; i < 100; i++) {
        std::string const rand1 = make_literal_random(digits, literal_radix);
        std::string const rand2 = make_literal_random(digits, literal_radix);
        std::string const nines = make_literal_nines(digits, literal_radix);
        std::string const zeros = make_literal_one_zeros(digits, literal_radix);
        output_with_signs_and_opts(rand1, rand2);
        output_with_signs_and_opts(rand1, nines);
        output_with_signs_and_opts(rand1, zeros);
        output_with_signs_and_opts(rand2, nines);
        output_with_signs_and_opts(rand2, zeros);
        output_with_signs_and_opts(nines, zeros);
      }
    }

    std::vector<std::size_t> border_digits = {
      elem_digits - 1,
      elem_digits,
      elem_digits + 1,
      elem_digits * 2,
      elem_digits * 2 + 1
    };
    for (std::size_t len_a : border_digits) {
      for (std::size_t len_b : border_digits) {
        for (std::size_t i = 0; i < 100; i++) {
          output_with_signs_and_opts(
            make_literal_random(len_a, literal_radix),
            make_literal_random(len_b, literal_radix));
        }
      }
    }

    const std::string values[7] = {
      "0", "1", "2",
      make_literal_nines(elem_digits, literal_radix),
      make_literal_one_zeros(elem_digits, literal_radix),
      make_literal_nines(elem_digits * 2, literal_radix),
      make_literal_one_zeros(elem_digits * 2, literal_radix)};
    for (std::size_t i = 0; i < std::size(values); i++) {
      for (std::size_t j = i; j < std::size(values); j++)
        output_with_signs_and_opts(values[i], values[j]);
      for (std::size_t len: test_digits)
        for (std::size_t j = 0; j < 100; j++)
          output_with_signs_and_opts(values[i], make_literal_random(len, literal_radix));
      for (std::size_t len: border_digits)
        for (std::size_t j = 0; j < 100; j++)
          output_with_signs_and_opts(values[i], make_literal_random(len, literal_radix));
    }

    return true;
  }

  void arith_run(const char* input_filename, const char* output_filename) {
    std::ifstream istr(input_filename);
    if (!istr) {
      std::cerr << "BUG(arith_run): failed to open the file '" << input_filename << "' to read." << std::endl;
      std::exit(3);
    }

    std::ofstream ostr(output_filename);
    if (!ostr) {
      std::cerr << "error(arith_run): failed to open the file '" << input_filename << "' to write." << std::endl;
      std::exit(1);
    }

    std::string line;
    while (std::getline(istr, line)) {
      std::istringstream sstr(line);
      std::string astr, op, bstr;
      sstr >> astr >> op >> bstr;
      int_t a(astr.c_str());
      int_t b(bstr.c_str());
      int_t r1, r2, r3, r4;
      switch (op[0]) {
      case '+':
        {
          r1 = a + b;
          r2 = a; r2 += b;
          int_t r3t = a + b; r3 = r3t;
          int_t r4t = a; r4t += b; r4 = r4t;
        }
        break;
      case '-':
        {
          r1 = a - b;
          r2 = a; r2 -= b;
          int_t r3t = a - b; r3 = r3t;
          int_t r4t = a; r4t -= b; r4 = r4t;
        }
        break;
      case '*':
        {
          r1 = a * b;
          r2 = a; r2 *= b;
          int_t r3t = a * b; r3 = r3t;
          int_t r4t = a; r4t *= b; r4 = r4t;
        }
        break;
      case '/':
        {
          r1 = a / b;
          r2 = a; r2 /= b;
          int_t r3t = a / b; r3 = r3t;
          int_t r4t = a; r4t /= b; r4 = r4t;
        }
        break;
      case '%':
        {
          r1 = a % b;
          r2 = a; r2 %= b;
          int_t r3t = a % b; r3 = r3t;
          int_t r4t = a; r4t %= b; r4 = r4t;
        }
        break;
      default:
        mwg_assert(0, "unknwon op = '%s'", op.c_str());
        break;
      }
      mwg_check(r1 == r2);
      mwg_check(r3 == r4);
      mwg_check(r1 == r3);
      ostr << r1 << std::endl;
    }
  }

  void arith_check() {
    std::error_code ec;
    std::filesystem::create_directories("../out/test", ec);
    if (ec) {
      std::cerr << "failed to create the directory ../out/test: " << ec.message() << std::endl;
      std::exit(1);
    }

    arith_save_cases("../out/test/big_integer.input.txt");
    std::system("BC_LINE_LENGTH=0 bc < ../out/test/big_integer.input.txt > ../out/test/big_integer.expect.txt");
    arith_run("../out/test/big_integer.input.txt", "../out/test/big_integer.result.txt");
    std::system("cd ../out/test; diff -bwu big_integer.{expect,result}.txt > big_integer.diff || echo \"see $PWD/big_integer.diff\"");
    // arith_run("../out/test/big_integer.input.1.txt", "../out/test/big_integer.result.1.txt");
  }

}

//-----------------------------------------------------------------------------

int main() {
  early_tests::run_tests();
  factorize_tests::make_table_for_squares();
  factorize_tests::check_primes_up_to();
  basic_tests::arith_check();
  return 0;
}
