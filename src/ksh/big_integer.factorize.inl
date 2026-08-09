// -*- mode: c++ -*-
#ifndef kashiwa_biginteger_factorize_inl
#define kashiwa_biginteger_factorize_inl
#include <cstdlib>
#include <cstdint>
#include <algorithm>
#include <limits>
#include <map>
#include <random>
#include <vector>
#include <mwg/except.h>
#include "big_integer.hpp"

namespace kashiwa {

  std::vector<std::int32_t> primes_up_to(std::int32_t B);

namespace {

  constexpr std::uint32_t small_prime_threshold = 10000;

  template<typename E, typename C, C M>
  void factorize_large(std::map<big_integer<E, C, M>, unsigned>& factors, big_integer<E, C, M>& value);

  template<typename E, typename C, C M>
  void factorize_product_gcd(std::map<big_integer<E, C, M>, unsigned>& factors, big_integer<E, C, M>& value, E next_prime_max) {
    auto _count_prime = [] (E x, E p) -> unsigned {
      unsigned power = 0;
      for (;;) {
        int const q = x / p;
        int const r = x % p;
        if (r != 0) break;
        power++;
        x = q;
      }
      return power;
    };

    std::vector<E> primes;
    std::vector<E> primes_new;
    primes.reserve(9);
    primes.emplace_back(2);
    primes.emplace_back(3);
    E prime_product = (E) 2 * 3;

    // The following variable contains a prime candidate that will be added to
    // "primes" next.  $p_k$ is a sequence of prime candidates. The first two
    // terms are special: $p_0 = 2, p_1 = 3$.  The other terms are given by
    // $p_{2k} = 6k - 1$ and $p_{2k+1} = 6k +1$.  The initial value is $k=2$ and
    // $p_2 = 5$.
    std::size_t next_prime_index = 2;
    E next_prime = 5;

    while (next_prime <= next_prime_max) {
      while (prime_product <= std::numeric_limits<E>::max() / next_prime) {
        primes.emplace_back(next_prime);
        prime_product *= next_prime;
        next_prime += next_prime_index++ % 2 == 0 ? 2 : 4;
      }
      if (value < primes[0] * primes[0]) break;

      for (;;) {
        E const g = kashiwa::gcd(value, prime_product);
        if (g <= 1) break;
        value /= g;

        unsigned max_power = 0;
        primes_new.clear();
        primes_new.reserve(primes.size());
        prime_product = 1;
        for (E p: primes) {
          unsigned const power = _count_prime(g, p);
          if (power == 0) continue;
          factors[p] += power;

          if (power > max_power) {
            max_power = power;
            primes_new.clear();
            prime_product = 1;
          }
          primes_new.emplace_back(p);
          prime_product *= p;
        }
        primes.swap(primes_new);
        mwg_assert(max_power >= 1, "max_power=%d", (int) max_power);

        E const prime_product_1 = prime_product;
        while (prime_product < std::numeric_limits<E>::max() / prime_product_1)
          prime_product *= prime_product_1;
      }

      primes.clear();
      prime_product = 1;
    }
  }

  template<typename E, typename C, C M>
  bool is_possibly_square(big_integer<E, C, M> const& value) {
    using int_t = big_integer<E, C, M>;
    static const std::uint64_t mod255_exclude_table[] = {0xef75ffebbdd67dec, 0xfcf77befbfcdef9a, 0xdbffb5bffb7cfe7f, 0xffbefbdde7ffcfe7, };
    static const std::uint64_t mod256_exclude_table[] = {0xfdfdfdedfdfcfdec, 0xfdfdfdedfdfdfdec, 0xfdfdfdedfdfcfded, 0xfdfdfdedfdfdfded, };
    static const std::uint64_t mod257_exclude_table[] = {0x81e9abe2191854e8, 0xe0897ee36ced7ca6, 0, };

    mwg_assert(value.data.size() >= 1);

    // mod 256
    std::uint32_t const mod256 = value.data[0] & 0xFF;
    if (mod256_exclude_table[mod256 / 64] & ((std::uint64_t) 1 << (mod256 % 64))) return false;

    // mod 2^32-1 ⇔ mod (2^16+1)(2^8+1)(2^8-1) ⇔ mod (65537*257*255)
    mwg_check(value.data.size() <= std::numeric_limits<std::uint64_t>::max() / int_t::element_max, "too large");
    std::uint64_t sum32 = 0;
    for (std::uint32_t d: value.data) sum32 += d;

    // mod255
    std::uint32_t const mod255 = sum32 % 255;
    if (mod255_exclude_table[mod255 / 64] & ((std::uint64_t) 1 << (mod255 % 64))) return false;

    // mod257: The table size of mod257 can be reduced to half because the
    // table is known to be self conjugate (i.e., table[i] = table[257 - i]).
    // This is because -1 is a perfect square with mod 257 (i.e., there exists
    // an integer $z$ that satisfies $z \equiv -1$), so if $n$ s.t. $n^2 \equiv
    // i$ exists, $n' = zn$ s.t. $n'^2 = (zn)^2 \equiv M - i$ exists.  However,
    // this doesn't happen with the other case with mod255 and mod256.
    std::uint32_t mod257 = sum32 % 257;
    if (mod257 > 128) mod257 = 257 - mod257;
    if (mod257_exclude_table[mod257 / 64] & ((std::uint64_t) 1 << (mod257 % 64))) return false;
    return true;
  }

  template<typename E, typename C, C M>
  bool factorize_square(std::map<big_integer<E, C, M>, unsigned>& factors, big_integer<E, C, M>& value) {
    using int_t = big_integer<E, C, M>;

    if (!is_possibly_square(value)) return false;
    mwg_assert(value.data.size() >= 1);

    std::size_t const o = (value.data.size() - 1) / 2;

    double a = value.data[2 * o];
    if ((value.data.size() & 1) == 0) {
      a += int_t::modulo * value.data[2 * o + 1];
    }
    if (o >= 1)
      a += value.data[2 * o - 1] / int_t::modulo;
    a = std::sqrt(a);

    int_t x;
    x.sign = 1;
    x.data.resize(o + 1, (E) 0);
    x.data[o] = (E) a;
    if (o >= 1)
      x.data[o - 1] = std::min<E>((a - x.data[o]) * int_t::modulo, int_t::element_max);

    int_t last_x;
    while (true) {
      int_t next_x = (x + value / x) / (E) 2;
      if (next_x == x || next_x == last_x) break;
      last_x = std::move(x);
      x = std::move(next_x);//
    }

    int const c = kashiwa::compare(x * x, value);
    if (c != 0) {
      x += c > 0 ? -1 : 1;
      if (kashiwa::compare(x * x, value) != 0) return false;
    }

    //std::cout << "value=" << value << ", x=" << x << std::endl;
    value = 1;
    std::map<int_t, unsigned> xfactors;
    factorize_large(xfactors, x);
    for (auto const& pair: xfactors)
      factors[pair.first] += 2 * pair.second;
    return true;
  }

  template<typename E, typename C, C M>
  E urand_mask(big_integer<E, C, M> const& n) {
    if (n.data.empty()) return 0;
    E const back = n.data.back();
    C mask = 1;
    while (mask <= back) mask <<= 1;
    mask--;
    return (E) mask;
  }

  template<typename E, typename C, C M>
  void urand_gen(big_integer<E, C, M>& r, big_integer<E, C, M> const& n, E mask) {
    static std::random_device rd;
    static std::mt19937 rng(rd());
    r.data.reserve(n.data.size());
    do {
      r.data.clear();
      std::generate_n(std::back_inserter(r.data), n.data.size(), std::ref(rng));
      r.data.back() &= mask;
      while (!r.data.empty() && r.data.back() == 0) r.data.pop_back();
      r.sign = r.data.empty() ? 0 : 1;
    } while (!(r < n));
  }

  template<typename E, typename C, C M>
  bool is_miller_rabin_prime(big_integer<E, C, M> const& value) {
    using int_t = big_integer<E, C, M>;

    static const E bases[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
    mwg_assert(value >= 4 && value > bases[sizeof bases / sizeof bases[0] - 1]);

    int_t d = value - 1;
    std::uint64_t s = 0;
    while (d % 2 == 0) {
      d /= 2;
      s++;
    }

    int_t const value_minus_one = value - 1;
    int_t x; // Declare outside to reuse the buffer.
    auto _check = [&value, &d, s, &value_minus_one, &x] (auto const& a) {
      x = kashiwa::ipow_mod<int_t>(a, d, value);

      // If $a^d \equiv \pm 1 (\mod N)$, OK.
      if (x == 1 || x == value_minus_one) return true;

      // Requires at least one $r$ s.t. $(a^d)^(2^r) \equiv -1 (\mod value)$
      for (std::uint64_t r = 1; r < s; ++r) {
        x = (x * x) % value;
        if (x == value_minus_one) return true;
        if (x == 1) return false;
      }
      return false;
    };

    if (value.data.size() <= 5) {
      for (E a: bases)
        if (!_check(a))
          return false;
    } else {
      // An appropriate number of random tests is estimated by the following
      // formula.  This is determined based on my observation on a table of
      // pairs of $(n_\mathrm{bits}, n_\mathrm{test})$.
      int ncheck = 3 + 100 / value.data.size();

      E const mask = urand_mask(value);
      int_t a;
      while (ncheck) {
        // generate a random $a \in [2, value -2]$
        urand_gen(a, value_minus_one, mask);
        if (a <= 1) continue;

        if (!_check(a)) return false;
        ncheck--;
      }
    }

    // This is almost certainly a prime.
    return true;
  }

  template<typename E, typename C, C M>
  bool factorize_pollard_rho(std::map<big_integer<E, C, M>, unsigned>& factors, big_integer<E, C, M>& value) {
    using int_t = big_integer<E, C, M>;
    constexpr int max_steps = 100000;
    constexpr int gcd_steps = 200;
    constexpr E max_orbits = 9;

    mwg_assert(value % 2 != 0);

    // stepping by f(t) = (t*t + c) % value
    auto _step_f = [&value] (int_t& t, E c) {
      t *= t;
      t += c;
      t %= value;
    };

    int_t x, y, g;
    auto _next_g = [&_step_f] (int_t& x, int_t& y, int_t& g, E c) {
      _step_f(x, c);
      _step_f(y, c);
      _step_f(y, c);
      g = x - y;
      if (g.sign < 0) kashiwa::chneg(g);
    };

    int_t x0, y0, product;
    for (E c = 1; c <= max_orbits; ++c) {
      x = 2;
      y = 2;
      g = 1;

      x0 = x;
      y0 = y;
      product = 1;
      for (int step = 0; step < max_steps; step++) {
        _next_g(x, y, g, c);
        product *= g;
        product %= value;

        // We check the termination only for every `gcd_steps`.
        if ((step + 1) % gcd_steps != 0 && step + 1 != max_steps) continue;

        g = kashiwa::gcd(product, value);
        if (g == 1) {
          x0 = x;
          y0 = y;
          product = 1;
          continue;
        }

        if (g == value) {
          // We try to go back to the first step where "gcd(g, n)" became
          // different from 1.
          x = x0;
          y = y0;
          do {
            _next_g(x, y, g, c);
            g = kashiwa::gcd(g, value);
          } while (g == 1);

          // If `g` is still `value` in the first step where "gcd(g, n) !=
          // 1", we conclude that this orbit reaches the infinite point by
          // the trivial factor "value".  We discard this orbit.
          if (g == value) break;
        }

        mwg_assert(value % g == 0, "BUG: this should not happen");
        value /= g;
        factorize_large(factors, g);
        return true;
      }
    }

    return false;
  }

  // We currently use 20000 for B.  The Google AI suggested 20000..50000 for
  // search in a few minutes.  If we want to search more aggressively, we may
  // increase this number to a very large value.  In such a case, the type
  // holding `B`, `ec_mul_t` should be changed to `std::uint64_t` or a bigger
  // one, though I'm not sure if the computation time would be realistic.
  typedef std::int32_t ec_mul_t;

  template<typename int_t>
  struct ec_point { int_t x, y; };

  // This function calculates p1 += p2.  When we encounter a non-trivial factor
  // in the middle of the process, this function cancels the operation and
  // returns the factor.  Otherwise, this function completes the addition and
  // returns 0.
  template<typename E, typename C, C M>
  big_integer<E, C, M> ec_chadd(ec_point<big_integer<E, C, M>>& p1, ec_point<big_integer<E, C, M>> const& p2, const big_integer<E, C, M>& a, const big_integer<E, C, M>& n) {
    using int_t = big_integer<E, C, M>;

    int_t num, den;
    if (&p1 == &p2 || (p1.x == p2.x && p1.y != p2.y)) {
      num = 3 * p1.x * p1.x + a;
      den = 2 * p1.y;
    } else {
      num = p2.y - p1.y;
      den = p2.x - p1.x;
    }
    kashiwa::chmod(num, n);
    kashiwa::chmod(den, n);

    int_t const den_inv = kashiwa::inv_mod(den, n);
    if (den_inv == 0) return kashiwa::gcd(den, n);

    num *= den_inv;
    kashiwa::chmod(num, n);

    int_t const xnew = kashiwa::mod(num * num - p1.x - p2.x, n);
    p1.y = kashiwa::mod(num * (p1.x - xnew) - p1.y, n);
    p1.x = std::move(xnew);
    return 0;
  }

  template<typename E, typename C, C M>
  big_integer<E, C, M> ec_chmul(ec_point<big_integer<E, C, M>>& p, ec_mul_t mul, big_integer<E, C, M> const& a, big_integer<E, C, M> const& n) {
    using int_t = big_integer<E, C, M>;

    if (mul <= 0) return 1;

    while (mul % 2 == 0) {
      if (int_t const g = ec_chadd(p, p, a, n)) return g;
      mul /= 2;
    }

    ec_point<int_t> base = p;
    while ((mul /= 2) > 0) {
      if (int_t const g = ec_chadd(base, base, a, n)) return g;
      if (mul % 2) {
        if (int_t const g = ec_chadd(p, base, a, n)) return g;
      }
    }
    return 0;
  }

  template<typename E, typename C, C M>
  bool factorize_elliptic_curve(std::map<big_integer<E, C, M>, unsigned>& factors, big_integer<E, C, M>& value) {
    using int_t = big_integer<E, C, M>;
    constexpr ec_mul_t B = 20000;

    static const std::vector<ec_mul_t> primes = primes_up_to(B);

    E const mask = urand_mask(value);
    for (int trial = 0; trial < 50; ++trial) {
      // Initialize $p$ and parameter $a$
      ec_point<int_t> p;
      urand_gen(p.x, value, mask);
      urand_gen(p.y, value, mask);
      int_t a;
      urand_gen(a, value, mask);
      for (E prime: primes) {
        ec_mul_t mul = prime;
        while (mul * prime <= B) mul *= prime;
        int_t g = ec_chmul(p, mul, a, value);
        if (g > 1 && g < value) {
          // When a non-trivial factor is found
          mwg_assert(value % g == 0, "BUG: this should not happen");
          value /= g;
          factorize_large(factors, g);
          return true;
        }
        if (g == value) break;
      }
    }

    return false;
  }

  template<typename E, typename C, C M>
  void factorize_large(std::map<big_integer<E, C, M>, unsigned>& factors, big_integer<E, C, M>& value) {
    using int_t = big_integer<E, C, M>;
    static std::vector<int_t> prime_cache;
    std::size_t icache = 0;

    // We are supposed to have already excluded the prime factors $p (\le
    // n_\mathrm{thresh})$, so if the value is less than $n_\mathrm{thresh}^2$,
    // it implies that value is a prime.
    while (value >= small_prime_threshold * small_prime_threshold) {
      // We scan the known big primes first.  The elements before `icache` has
      // already been tested for the present `value`, so we can start the test
      // from prime_cache[icache].  It should be noted that new primes may be
      // appeded in the nested calls of `factorize_large`.
      if (icache < prime_cache.size()) {
        for (; icache < prime_cache.size(); icache++) {
          int_t q;
          if (kashiwa::div(value, prime_cache[icache], q) == 0) {
            value = q;
            factors[prime_cache[icache]]++;
            break;
          }
        }

        // The following condition implies that a factor is found and the above
        // loop was canceled, so let us check the lower bound of `value` again.
        if (icache < prime_cache.size()) continue;
      }

      // When "factorize_square" succeeds, value is supposed to become 1, so
      // the factorization is complete.
      if (factorize_square(factors, value)) return;

      // When the Miller-Rabin test considers it to be a prime number, it is
      // probably a prime number, so we stop here.
      if (is_miller_rabin_prime(value)) break;

      // When any of the following methods succeed to find a factor, we restart
      // on the remaining "value" from the beginning.
      if (factorize_pollard_rho(factors, value)) continue;
      if (factorize_elliptic_curve(factors, value)) continue;

      // When none succeeds any more, we give up.
      break;
    }

    if (value != 1) {
      if (std::find(prime_cache.begin() + icache, prime_cache.end(), value) == prime_cache.end())
        prime_cache.emplace_back(value);

      factors[value]++;
      value = 1;
    }
  }
}

template<typename E, typename C, C M>
std::map<big_integer<E, C, M>, unsigned> big_integer<E, C, M>::factorize(big_integer<E, C, M> value) {
  using int_t = big_integer<E, C, M>;

  std::map<int_t, unsigned> ret;
  if (value < 0) {
    ret.emplace((int_t) -1, 1);
    kashiwa::chneg(value);
  }

  if (value <= 3) {
    if (value != 1)
      ret.emplace((int_t) value, 1);
    return ret;
  }

  factorize_product_gcd(ret, value, small_prime_threshold);
  factorize_large(ret, value);
  return ret;
}

}

#endif
