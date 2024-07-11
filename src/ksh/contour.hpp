#ifndef kashiwa_contour_hpp
#define kashiwa_contour_hpp
#include <cstddef>

#include <complex>
#include <stdexcept>
#include <vector>

// * XXX--領域が単連結でない場合には内側の boundary を正しく抽出する事ができない。
//
// * XXX--現在の実装だと、一点から4以上の線が出ている場合に無限ループになった
//   り取り尽くしていない点が発生したりする気がする。ちゃんと確認する必要がある。
//
// * ToDo: 現在 limitx, limity によって生じた boudary も出力しているが、これを
//   省略する機能?

namespace kashiwa::contour {

class binning {
  std::size_t m_size = 1u;
  double m_min = 0.0;
  double m_max = 1.0;
  double m_dx = 1.0;
public:
  std::size_t size() const { return m_size; }
  double min() const { return m_min; }
  double max() const { return m_max; }
  double mesh() const { return m_dx; }
public:
  void initialize(double const min_value, double const max_value, std::size_t const count) {
    m_size = count;
    m_min = min_value;
    m_max = max_value;
    m_dx = (max_value - min_value) / count;
  }
  void initialize(double const min_value, std::size_t const count, double const dx) {
    m_size = count;
    m_min = min_value;
    m_max = min_value + count * dx;
    m_dx = dx;
  }
  binning() {}
  binning(double const min_value, double const max_value, std::size_t const count) {
    initialize(min_value, max_value, count);
  }
  binning(double const min_value, std::size_t const count, double const dx) {
    initialize(min_value, count, dx);
  }
  double at(double const index) const {
    return m_min + index * m_dx;
  }
  double operator[](std::size_t const index) const {
    return at(index);
  }
  double mid_point(std::size_t index) const {
    return at(index + 0.5);
  }
};

struct contour_search_params {
  binning binx;
  binning biny;

  std::pair<double, double> limitx = {
    std::numeric_limits<double>::quiet_NaN(),
    std::numeric_limits<double>::quiet_NaN()
  };
  std::pair<double, double> limity = {
    std::numeric_limits<double>::quiet_NaN(),
    std::numeric_limits<double>::quiet_NaN()
  };
};

template<typename T, typename F>
T binary_search(T lower, T upper, F predicate, int count = 10) {
  bool const up = predicate(lower);
  bool const lp = predicate(upper);
  if (up == lp)
    throw std::runtime_error("binary_search: wrong boundary condition");

  while (count--) {
    if (T m =(T)0.5 * (lower + upper); predicate(m) == up) {
      lower = m;
    } else {
      upper = m;
    }
  }

  return (T)0.5 * (lower + upper);
}

template<typename F>
struct contour_tracer {
  contour_search_params m_params;
  F m_func;
  contour_tracer(contour_search_params const& params, F func):
    m_params(params), m_func(func)
  {
    if (std::isnan(m_params.limitx.first))
      m_params.limitx.first = 1.5 * m_params.binx.min() - 0.5 * m_params.binx.max();
    if (std::isnan(m_params.limitx.second))
      m_params.limitx.second = 1.5 * m_params.binx.max() - 0.5 * m_params.binx.min();
    if (std::isnan(m_params.limity.first))
      m_params.limity.first = 1.5 * m_params.biny.min() - 0.5 * m_params.biny.max();
    if (std::isnan(m_params.limity.second))
      m_params.limity.second = 1.5 * m_params.biny.max() - 0.5 * m_params.biny.min();
  }

  int evaluate(double const x, double const y) {
    if (x < m_params.limitx.first || x > m_params.limitx.second) return -1;
    if (y < m_params.limity.first || y > m_params.limity.second) return -1;
    return m_func(x, y);
  }

  // -1: already checked
  // 0: region not to check
  // 1: region to check
  std::vector<int> m_table;
  void generate_map() {
    std::size_t const ixN = m_params.binx.size();
    std::size_t const iyN = m_params.biny.size();
    m_table.resize((ixN + 1) * (iyN + 1));
    for (std::size_t ix = 0; ix <= ixN; ix++) {
      double const x = m_params.binx[ix];
      for (std::size_t iy = 0; iy <= iyN; iy++) {
        double const y = m_params.biny[iy];
        m_table[ix * (iyN + 1) + iy] = evaluate(x, y);
      }
    }
  }

  std::vector<std::complex<double>> m_region_points;
  void find_connected_regions() {
    int const ixN = m_params.binx.size();
    int const iyN = m_params.biny.size();
    double const dx = m_params.binx.mesh();
    double const dy = m_params.biny.mesh();
    std::vector<int> table = m_table; // copy

    std::vector<std::pair<int, int>> list1, list2;
    for (int ix = 0; ix <= ixN; ix++) {
      for (int iy = 0; iy <= iyN; iy++) {
        int const color = table[ix * (iyN + 1) + iy];
        if (color <= 0) continue;

        int ixmax = -1, iymax = -1;
        auto propagate_adjacent = [&table, &list2, &ixmax, &iymax, ixN, iyN, color] (int jx, int jy) {
          if (jx < 0 || ixN < jx) return;
          if (jy < 0 || iyN < jy) return;
          if (table[jx * (iyN + 1) + jy] != color) return;
          if (jx > ixmax) ixmax = jx, iymax = jy;
          table[jx * (iyN + 1) + jy] = -1;
          list2.emplace_back(jx, jy);
        };

        // Note: We here use BFS to avoid stack overflow. With DFS, as many
        // recursion level as the area is required.
        list1.clear();
        list1.push_back(std::make_pair(ix, iy));
        propagate_adjacent(ix, iy);
        while (list1.size()) {
          list2.clear();
          for (auto [ix1, iy1]: list1) {
            propagate_adjacent(ix1 + 1, iy1);
            propagate_adjacent(ix1 - 1, iy1);
            propagate_adjacent(ix1, iy1 + 1);
            propagate_adjacent(ix1, iy1 - 1);

            // 斜めも
            propagate_adjacent(ix1 + 1, iy1 + 1);
            propagate_adjacent(ix1 + 1, iy1 - 1);
            propagate_adjacent(ix1 - 1, iy1 + 1);
            propagate_adjacent(ix1 - 1, iy1 - 1);
          }
          list1.swap(list2);
        }

        double const x = m_params.binx[ixmax];
        double const y = m_params.biny[iymax];
        m_region_points.emplace_back(x, y);
      }
    }
  }

  std::vector<std::vector<std::complex<double>>> get_contours() {
    std::vector<std::vector<std::complex<double>>> ret;
    constexpr int scope_resolution = 100;
    constexpr double scope_offset = 1e-5;

    double const dx = m_params.binx.mesh();
    double const ds = 0.5 * dx;

    int color = 0;
    auto _in_region = [this, &color] (std::complex<double> value) -> bool {
      return evaluate(value.real(), value.imag()) == color;
    };

    for (std::complex<double> z0: m_region_points) {
      try {
        color = evaluate(z0.real(), z0.imag());

        // Identify the smallest x where (x, Im(z0)) is in the region.  Since
        // we clip the region by m_params.limitx, we should be able to find
        // such a point at least until we reach the right edge of the clipping
        // box, m_params.limitx.second.
        std::complex<double> zr = z0 + 2.0 * dx;
        while (zr.real() < m_params.limitx.second + 2.0 * dx && _in_region(zr))
          zr += 2.0 * dx;
        z0 = binary_search(z0, zr, _in_region, 32);

        std::complex<double> zt = z0 - ds;
        while (!_in_region(zt))
          zt = 0.5 * (z0 + binary_search(zt, z0 + (zt - z0) * 1e-6, _in_region, 20));

        double theta = binary_search(0.0, M_PI, [z0, zt, &_in_region] (double theta) {
          return _in_region(z0 + (zt - z0) * std::polar(1.0, theta));
        }, 20);
        std::complex<double> z1 = z0;
        std::complex<double> z2 = z0 + (zt - z0) * std::polar(1.0, theta);

        std::vector<std::complex<double>> contour;
        contour.push_back(z1);
        contour.push_back(z2);

        for (;;) {
          std::complex<double> n = (z2 - z1) / std::abs(z2 - z1);

          double scope;
          for (int i = 1; i <= scope_resolution; i++) {
            scope = (M_PI - scope_offset) * i / scope_resolution;
            if (_in_region(z2 + n * std::polar(ds, -scope))
              != _in_region(z2 + n * std::polar(ds, scope)))
              break;
          }

          theta = binary_search(-scope, scope, [z2, n, ds, &_in_region] (double theta) {
            return _in_region(z2 + n * std::polar(ds, theta));
          }, 20);
          z1 = z2;
          z2 = z1 + n * std::polar(ds, theta);
          contour.push_back(z2);

          if (std::abs(z2 - z0) < 2.0 * ds && contour.size() >= 5) {
            contour.push_back(z0);
            break;
          }
          if (contour.size() >= 10000) break;
        }

        if (contour.size())
          ret.emplace_back(std::move(contour));
      } catch (std::runtime_error& e) {
        std::fprintf(stderr, "contour:region(%g,%g): failed in binary search\n", z0.real(), z0.imag());
      }
    }
    return ret;
  }
  std::size_t print_contours(std::FILE* file) {
    auto contours = get_contours();
    for (auto contour: contours) {
      for (auto z: contour)
        std::fprintf(file, "%g %g\n", z.real(), z.imag());
      std::fprintf(file, "\n");
    }
    return contours.size();
  }
  void save_contours(const char* filename) {
    std::FILE* file = std::fopen(filename, "wb");
    if (!file) {
      std::fprintf(stderr, "countour_generator(save_contours): failed to open the output file '%s'\n");
      throw std::runtime_error("countour_generator(save_contours): failed to open the output file");
    }
    std::size_t const count = print_contours(file);
    std::fclose(file);
    std::printf("%s: saved %zu contours\n", filename, count);
  }
};

template<typename F>
void save_contours(const char* filename, contour_search_params const& params, F func) {
  contour_tracer gen(params, func);
  gen.generate_map();
  gen.find_connected_regions();
  gen.save_contours(filename);
}

}
#endif
