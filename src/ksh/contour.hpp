#ifndef kashiwa_contour_hpp
#define kashiwa_contour_hpp
#include <cstddef>
#include <cstdint>

#include <complex>
#include <stdexcept>
#include <vector>

#include "def.hpp"
#include "utility.hpp"

// tracer_inchworm:
//
// * XXX--領域が単連結でない場合には内側の boundary を正しく抽出する事ができない。
//
// * XXX--現在の実装だと、一点から4以上の線が出ている場合に無限ループになったり
//   取り尽くしていない点が発生したりする気がする。ちゃんと確認する必要がある。

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

  // From value to index
  double locate(double value) const {
    return (value - m_min) / m_dx;
  }
};

struct contour_search_params {
  binning binx;
  binning biny;
  bool closepath = false;

  enum tracer_type {
    tracer_inchworm,
    tracer_grid,
  };
  tracer_type tracer = tracer_grid;

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
        m_table[ix * (iyN + 1) + iy] = m_func(x, y);
      }
    }
  }

  std::vector<std::complex<double>> m_region_points;
  void update_connected_regions() {
    int const ixN = m_params.binx.size();
    int const iyN = m_params.biny.size();
    double const dx = m_params.binx.mesh();
    double const dy = m_params.biny.mesh();
    std::vector<int> table = m_table; // copy

    m_region_points.clear();

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

public:
  /*?lwiki
   * @param[in] std::complex<double> const& z;
   *   The current position.
   * @param[in] std::complex<double> const& n;
   *   The outward direction for the previous line element (from the current
   *   position `z` to the previous position).
   * @param[in] double ds;
   *   The search radius, or equivalently, the length of the line element.
   * @param[in] Predicate const& in_region;
   *   The predicate to test whether the specified  point is in the region.
   */
  template<typename Predicate>
  std::pair<double, double> determine_angle_scope(std::complex<double> const& z, std::complex<double> const& n, double const ds, Predicate const& in_region) const {
#if 0
    // D0020: To properly handle critical points at which multiple contours
    // cross, we first scan [0, 2\pi) with some sampling points. これにより確か
    // に臨界点上での振る舞いは改善したが、代わりに別の場所で角度の二分探索に失
    // 敗する。
    {
      constexpr int angle_sample_count = 12;
      double const dtheta = 2.0 * M_PI / angle_sample_count;

      // We first determine a small angle $\theta$ outside the region.  We
      // should be able to find such an angle in the clockwise direction
      // because we trace the contour so that it encloses the region on the
      // left-hand side.
      int count = 0;
      double theta = 0.5 * dtheta;
      while (in_region(z + n * std::polar(ds, theta))) {
        theta *= 0.5;
        if (++count >= 20)
          throw std::runtime_error("binary_search: wrong boundary direction");
      }

      // We then increase $\theta$ by $d\theta$ and find the first $\theta$ in
      // the region.
      for (int i = 1; i < angle_sample_count; i++) {
        theta += dtheta;
        if (in_region(z + n * std::polar(ds, theta)))
          return {theta - dtheta, theta};
      }
    }
#endif

    constexpr int scope_resolution = 100;
    constexpr double scope_offset = 1e-5;
    // The above implementations should basically succeed if enabled, but in
    // case it fails due to the coarse sampling, we fallback to the older
    // approach, which selects the contour in the same direction as the
    // previous line element as much as possible.
    double scope;
    for (int i = 1; i <= scope_resolution; i++) {
      scope = (M_PI - scope_offset) * i / scope_resolution;
      if (in_region(z - n * std::polar(ds, -scope))
        != in_region(z - n * std::polar(ds, scope)))
        break;
    }

    return {M_PI - scope, M_PI + scope};
  }

  std::vector<std::vector<std::complex<double>>> trace_contours_inchworm() {
    this->update_connected_regions();

    std::vector<std::vector<std::complex<double>>> ret;

    double const dx = m_params.binx.mesh();
    double const ds = 0.5 * dx;

    int color = 0;
    auto _in_region = [this, &color] (std::complex<double> value) -> bool {
      double const x = value.real(), y = value.imag();
      if (x < m_params.limitx.first || x > m_params.limitx.second) return false;
      if (y < m_params.limity.first || y > m_params.limity.second) return false;
      return m_func(x, y) == color;
    };

    for (std::complex<double> z0: m_region_points) {
      try {
        color = m_func(z0.real(), z0.imag());

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

        double theta = binary_search(-M_PI, 0.0, [z0, zt, &_in_region] (double theta) {
          return _in_region(z0 + (zt - z0) * std::polar(1.0, theta));
        }, 20);
        std::complex<double> z1 = z0;
        std::complex<double> z2 = z0 + (zt - z0) * std::polar(1.0, theta);

        std::vector<std::complex<double>> contour;
        contour.push_back(z1);
        contour.push_back(z2);

        for (;;) {
          std::complex<double> n = (z1 - z2) / std::abs(z2 - z1);
          std::pair<double, double> const scope = this->determine_angle_scope(z2, n, ds, _in_region);
          theta = binary_search(scope.first, scope.second, [z2, n, ds, &_in_region] (double theta) {
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

public:
  enum grid_contour_point_type {
    limit_corner,
    link_up,
    link_down,
    link_right,
    link_left,
  };

  enum grid_link_mark_flags {
    grid_link_mark_vertical = 0x1,
    grid_link_mark_horizontal = 0x2,
  };

  void grid_initialize_link_mark(std::vector<std::uint32_t>& link_mark) const {
    int const ixN = m_params.binx.size();
    int const iyN = m_params.biny.size();
    link_mark.assign(ixN * iyN, 0);
    for (int ix = 0; ix < ixN; ix++) {
      for (int iy = 0; iy < iyN; iy++) {
        int const color0 = m_table[ix * (iyN + 1) + iy];
        int const color1 = m_table[(ix + 1) * (iyN + 1) + iy];
        int const color2 = m_table[ix * (iyN + 1) + iy + 1];
        if (color0 != color1 && (color0 > 0 || color1 > 0))
          link_mark[ix * iyN + iy] |= grid_link_mark_horizontal;
        if (color0 != color2 && (color0 > 0 || color2 > 0))
          link_mark[ix * iyN + iy] |= grid_link_mark_vertical;
      }
    }
  }

  void grid_unmark_link(std::vector<std::uint32_t>& link_mark, grid_contour_point_type type, int ix, int iy) const {
    int const iyN = m_params.biny.size();
    switch (type) {
    case link_down:
      iy--;
      ksh_fallthrough;
      /*FALL-THROUGH*/
    case link_up:
      link_mark[ix * iyN + iy] &= ~grid_link_mark_vertical;
      break;
    case link_left:
      ix--;
      ksh_fallthrough;
      /*FALL-THROUGH*/
    case link_right:
      link_mark[ix * iyN + iy] &= ~grid_link_mark_horizontal;
      break;
    }
  }

  void grid_trace_contour(std::vector<std::vector<std::complex<double>>>& out, std::vector<std::uint32_t>& link_mark, grid_contour_point_type type, int ix, int iy, bool include_first) {
    double const dx = m_params.binx.mesh();
    int const ixN = m_params.binx.size();
    int const iyN = m_params.biny.size();
    double const xthresh = dx * dx * 1e-4;

    std::vector<std::complex<double>> contour;
    auto _add_contour = [&out, &contour] {
      if (contour.size()) {
        out.emplace_back(std::move(contour));
        contour.clear();
      }
    };

    int const color = m_table[ix * (iyN + 1) + iy];
    auto _in_region = [this] (std::complex<double> value) -> bool {
      return this->m_func(value.real(), value.imag());
    };
    auto _index2contained = [this, iyN, color] (int ix, int iy) -> bool {
      return this->m_table[ix * (iyN + 1) + iy] == color;
    };
    auto _index2complex = [this] (int ix, int iy) -> std::complex<double> {
      return std::complex<double>(m_params.binx[ix], m_params.biny[iy]);
    };

    std::complex<double> p1, p2;
    std::complex<double> pinit, plast;
    bool pinit_initialized = false;

    if (include_first) {
      p1 = _index2complex(ix, iy);
      switch (type) {
      case link_right: p2 = _index2complex(ix + 1, iy); break;
      case link_left: p2 = _index2complex(ix - 1, iy); break;
      case link_up: p2 = _index2complex(ix, iy + 1); break;
      case link_down: p2 = _index2complex(ix, iy - 1); break;
      }
      plast = binary_search(p1, p2, _in_region, 20);
      contour.emplace_back(plast);
      grid_unmark_link(link_mark, type, ix, iy);
      pinit = plast;
      pinit_initialized = true;
    }

    for (int iloop = 0; iloop < 10000; iloop++) {
      switch (type) {
      case link_right:
        if (iy == iyN) {
          if (!m_params.closepath) _add_contour();
          while (ix > 0 && _index2contained(ix - 1, iy)) ix--;
          if (ix == 0) {
            type = link_up;
            goto corner;
          }
          type = link_left;
          p1 = _index2complex(ix, iy);
          p2 = _index2complex(ix - 1, iy);
        } else if (_index2contained(ix + 1, iy + 1)) {
          type = link_down;
          p1 = _index2complex(++ix, ++iy);
        } else if (_index2contained(ix, iy + 1)) {
          type = link_right;
          p1 = _index2complex(ix, ++iy);
          p2 = _index2complex(ix + 1, iy);
        } else {
          type = link_up;
          p2 = _index2complex(ix, iy + 1);
        }
        goto link;
      case link_up:
        if (ix == 0) {
          if (!m_params.closepath) _add_contour();
          while (iy > 0 && _index2contained(ix, iy - 1)) iy--;
          if (iy == 0) {
            type = link_left;
            goto corner;
          }
          type = link_down;
          p1 = _index2complex(ix, iy);
          p2 = _index2complex(ix, iy - 1);
        } else if (_index2contained(ix - 1, iy + 1)) {
          type = link_right;
          p1 = _index2complex(--ix, ++iy);
        } else if (_index2contained(ix - 1, iy)) {
          type = link_up;
          p1 = _index2complex(--ix, iy);
          p2 = _index2complex(ix, iy + 1);
        } else {
          type = link_left;
          p2 = _index2complex(ix - 1, iy);
        }
        goto link;
      case link_left:
        if (iy == 0) {
          if (!m_params.closepath) _add_contour();
          while (ix < ixN && _index2contained(ix + 1, iy)) ix++;
          if (ix == ixN) {
            type = link_down;
            goto corner;
          }
          type = link_right;
          p1 = _index2complex(ix, iy);
          p2 = _index2complex(ix + 1, iy);
        } else if (_index2contained(ix - 1, iy - 1)) {
          type = link_up;
          p1 = _index2complex(--ix, --iy);
        } else if (_index2contained(ix, iy - 1)) {
          type = link_left;
          p1 = _index2complex(ix, --iy);
          p2 = _index2complex(ix - 1, iy);
        } else {
          type = link_down;
          p2 = _index2complex(ix, iy - 1);
        }
        goto link;
      case link_down:
        if (ix == ixN) {
          if (!m_params.closepath) _add_contour();
          while (iy < iyN && _index2contained(ix, iy + 1)) iy++;
          if (iy == iyN) {
            type = link_right;
            goto corner;
          }
          type = link_up;
          p1 = _index2complex(ix, iy);
          p2 = _index2complex(ix, iy + 1);
        } else if (_index2contained(ix + 1, iy - 1)) {
          type = link_left;
          p1 = _index2complex(++ix, --iy);
        } else if (_index2contained(ix + 1, iy)) {
          type = link_down;
          p1 = _index2complex(++ix, iy);
          p2 = _index2complex(ix, iy - 1);
        } else {
          type = link_right;
          p2 = _index2complex(ix + 1, iy);
        }
        goto link;
      corner:
        plast = _index2complex(ix, iy);
        if (m_params.closepath)
          contour.emplace_back(plast);
        break;
      link:
        plast = binary_search(p1, p2, _in_region, 20);
        contour.emplace_back(plast);
        grid_unmark_link(link_mark, type, ix, iy);
        break;
      }

      if (!pinit_initialized) {
        pinit = plast;
        pinit_initialized = true;
      } else if (std::norm(plast - pinit) < xthresh) {
        break;
      }
    }

    _add_contour();
  }

  std::vector<std::vector<std::complex<double>>> trace_contours_grid() {
    int const ixN = m_params.binx.size();
    int const iyN = m_params.biny.size();

    std::vector<std::uint32_t> link_mark;
    grid_initialize_link_mark(link_mark);

    std::vector<std::vector<std::complex<double>>> ret;

    for (int ix = 0; ix < ixN; ix++) {
      for (int iy = 0; iy < iyN; iy++) {
        if (link_mark[ix * iyN + iy] & grid_link_mark_horizontal) {
          if (m_table[ix * (iyN + 1) + iy] > 0)
            grid_trace_contour(ret, link_mark, link_right, ix, iy, true);
          else
            grid_trace_contour(ret, link_mark, link_left, ix + 1, iy, true);
        }
        if (link_mark[ix * iyN + iy] & grid_link_mark_vertical) {
          if (m_table[ix * (iyN + 1) + iy] > 0)
            grid_trace_contour(ret, link_mark, link_up, ix, iy, true);
          else
            grid_trace_contour(ret, link_mark, link_down, ix, iy + 1, true);
        }
      }
    }

    return ret;
  }

public:
  std::vector<std::vector<std::complex<double>>> trace_contours() {
    switch (m_params.tracer) {
    case contour_search_params::tracer_grid:
      return  trace_contours_grid();
    case contour_search_params::tracer_inchworm:
      return  trace_contours_inchworm();
    default:
      throw std::logic_error("FATAL: unrecognized tracer_type");
    }
  }

  std::size_t print_contours(std::FILE* file) {
    auto contours = trace_contours();
    for (auto contour: contours) {
      for (auto z: contour)
        std::fprintf(file, "%g %g\n", z.real(), z.imag());
      std::fprintf(file, "\n");
    }
    return contours.size();
  }
  void save_contours(const char* filename) {
    std::string filename_part = filename;
    filename_part += ".part";
    std::FILE* file = std::fopen(filename_part.c_str(), "wb");
    if (!file) {
      std::fprintf(stderr, "countour_generator(save_contours): failed to open the output file '%s'\n");
      throw std::runtime_error("countour_generator(save_contours): failed to open the output file");
    }
    std::size_t const count = print_contours(file);
    std::fclose(file);
    std::remove(filename);
    std::rename(filename_part.c_str(), filename);
    std::printf("%s: saved %zu contours\n", filename, count);
  }
};

template<typename F>
void save_contours(const char* filename, contour_search_params const& params, F func) {
  contour_tracer gen(params, func);
  gen.generate_map();
  gen.save_contours(filename);
}

}
#endif
