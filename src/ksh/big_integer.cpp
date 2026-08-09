#include <cstdint>
#include <cstdlib>
#include <vector>
#include "arithmetic.hpp"
#include "big_integer.hpp"
#include "big_integer.factorize.inl"

namespace kashiwa {

  std::vector<std::int32_t> primes_up_to(std::int32_t B) {
    mwg_check(B >= 3);

    std::vector<std::int32_t> ret;

    // We first fill 2 and 3 and the prime candidates in the form $p = 6k \pm
    // 1$ (which are not divisable by 2 or 3).
    ret.reserve(2 + (B + 1)/ 6 + (B - 1) / 6);
    ret.emplace_back(2);
    ret.emplace_back(3);
    std::int32_t p = 1;
    for (;;) {
      if ((p += 4) > B) break;
      ret.emplace_back(p);
      if ((p += 2) > B) break;
      ret.emplace_back(p);
    }

    for (std::size_t i = 2; i < ret.size(); i++) {
      std::int32_t p = ret[i];
      if (p == 0) continue;
      if (p > B / p) break;

      // For the $(6k - 1)$-type candidate number, we want to find one $k$
      // satisfying $6k-1 \equiv 0 \mod p$.  This accomplished by $k = 6^{-1}$
      // in $\mod p$.
      std::int32_t k0 = kashiwa::inv_mod((std::int32_t) 6, p);

      // We want to start search from p^2.  When $p^2 = 6n ... 6n + 5$, the
      // minimum candidate $6k-1 (\ge p^2)$ is $6n + 5$, so $k = n + 1$.  Then,
      // we find the minimum $k_1 \equiv k_0$ that is larger than $n$.
      std::int32_t k1 = k0 + p * std::max((p * p / 6 - k0 + p) / p, (std::int32_t) 0);
      std::int32_t p1 = 6 * k1 - 1;
      for (; p1 <= B; k1 += p, p1 += 6 * p) ret[k1 * 2] = 0;

      // The new k0 is the solution for $(6k+1)$-type candidate: i.e., $6 k_0 +
      // 1 \equiv p$.  Then, when $p^2 = 6n + 2 ... 6n + 7$, the minimum
      // candidate $6k + 1 (\ge p^2)$ is $6n + 7$, so $k = n + 1$.  We find the
      // minimum $k_1 \equiv k_0$ that is larger than $n$.
      k0 = (p - k0) % p;
      k1 = k0 + p * std::max(((p * p  - 2) / 6 - k0 + p) / p, (std::int32_t) 0);
      p1 = 6 * k1 + 1;
      for (; p1 <= B; k1 += p, p1 += 6 * p) ret[k1 * 2 + 1] = 0;
    }

    ret.erase(
      std::remove_if(ret.begin(), ret.end(), [] (std::int32_t a) { return a == 0; }),
      ret.end());
    return ret;
  }

  // template instantiation
  template std::map<big_integer<>, unsigned> big_integer<>::factorize(big_integer<> value);

}
