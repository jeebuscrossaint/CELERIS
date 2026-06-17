#include "celeris/analysis/throughfocus.hpp"

#include "celeris/core.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

#ifdef CELERIS_USE_CUDA_KERNELS
#include "celeris/cuda/propagate.hpp"
#endif

namespace celeris {

ThroughFocus analyze_through_focus(const MetalensDesign& lens,
                                   const UnitCellLibrary& lib,
                                   double focal_length_um, double wavelength_um,
                                   double diameter_um) {
    const int n = lens.n_cells;
    const double p = lens.period_um;
    const double center = (n - 1) / 2.0;
    const double R = diameter_um / 2.0;
    const double k = 2.0 * pi / wavelength_um;

    struct Cell { double x, y; cdouble t; };
    std::vector<Cell> cells;
    for (int iy = 0; iy < n; ++iy)
        for (int ix = 0; ix < n; ++ix) {
            const double x = (ix - center) * p, y = (iy - center) * p;
            if (std::sqrt(x * x + y * y) > R) continue;
            const double fill = lens.fill_map[static_cast<std::size_t>(iy) * n + ix];
            cells.push_back({x, y, lib.transmission_for_fill(fill)});
        }

    auto on_axis = [&](double z) {
        cdouble E{0, 0};
        for (const Cell& c : cells) {
            double r = std::sqrt(c.x * c.x + c.y * c.y + z * z);
            E += c.t * std::polar(1.0 / r, k * r);
        }
        return std::norm(E);
    };

    // Scan a window around the nominal focus; widen it for low NA (long DOF).
    const double na = std::sin(std::atan(R / focal_length_um));
    const double half = std::max(2.0 * wavelength_um / (na * na + 1e-6), 0.25 * focal_length_um);
    const double zlo = std::max(0.1 * focal_length_um, focal_length_um - half);
    const double zhi = focal_length_um + half;
    const int N = 161;

    ThroughFocus tf;
    double peak = 0.0;
    std::vector<double> I(N);
    for (int i = 0; i < N; ++i) {
        double z = zlo + (zhi - zlo) * i / (N - 1);
        I[i] = on_axis(z);
        peak = std::max(peak, I[i]);
        tf.z_um.push_back(static_cast<float>(z));
    }
    if (peak <= 0) peak = 1;

    int ipk = 0;
    for (int i = 0; i < N; ++i) {
        tf.intensity.push_back(static_cast<float>(I[i] / peak));
        if (I[i] > I[ipk]) ipk = i;
    }
    tf.z_peak_um = tf.z_um[ipk];

    // Depth of focus = axial full width at half maximum of the peak.
    double half_lvl = peak / 2.0;
    int lo = -1, hi = -1;
    for (int i = 0; i < N; ++i)
        if (I[i] >= half_lvl) { if (lo < 0) lo = i; hi = i; }
    tf.dof_um = (hi > lo && lo >= 0) ? (tf.z_um[hi] - tf.z_um[lo]) : 0.0;

    // Longitudinal caustic (x-z intensity slice through the focus). The lateral
    // window is a few diffraction-limited spot widths; the axial window matches
    // the on-axis scan. Computed on the GPU when available (this is exactly the
    // dense propagation it excels at); a coarse CPU grid otherwise.
    const double dl = wavelength_um * focal_length_um / std::max(diameter_um, 1e-6);
    const double Wx = std::max(8.0 * dl, 3.0);
    tf.caustic_xmin = -Wx; tf.caustic_xmax = Wx;
    tf.caustic_zmin = zlo;  tf.caustic_zmax = zhi;

    std::vector<double> px(cells.size()), py(cells.size());
    std::vector<cdouble> pt(cells.size());
    for (std::size_t c = 0; c < cells.size(); ++c) {
        px[c] = cells[c].x; py[c] = cells[c].y; pt[c] = cells[c].t;
    }

    int cnx = 0, cnz = 0;
    std::vector<double> cmap;
#ifdef CELERIS_USE_CUDA_KERNELS
    {
        cnx = 201; cnz = 201;
        cmap.assign(static_cast<std::size_t>(cnx) * cnz, 0.0);
        if (!cuda::propagate_zx(px.data(), py.data(), pt.data(),
                                static_cast<int>(pt.size()), tf.caustic_xmin,
                                tf.caustic_xmax, cnx, tf.caustic_zmin,
                                tf.caustic_zmax, cnz, k, cmap.data())) {
            cnx = cnz = 0; cmap.clear();
        }
    }
#endif
    if (cnx == 0) {
        // CPU fallback: keep it small so it stays interactive without a GPU.
        cnx = 81; cnz = 81;
        cmap.assign(static_cast<std::size_t>(cnx) * cnz, 0.0);
        for (int iz = 0; iz < cnz; ++iz) {
            double z = tf.caustic_zmin + (tf.caustic_zmax - tf.caustic_zmin) * iz / (cnz - 1);
            for (int ix = 0; ix < cnx; ++ix) {
                double fx = tf.caustic_xmin + (tf.caustic_xmax - tf.caustic_xmin) * ix / (cnx - 1);
                cdouble E{0, 0};
                for (const Cell& c : cells) {
                    double r = std::sqrt((fx - c.x) * (fx - c.x) + c.y * c.y + z * z);
                    E += c.t * std::polar(1.0 / r, k * r);
                }
                cmap[static_cast<std::size_t>(iz) * cnx + ix] = std::norm(E);
            }
        }
    }

    double cpeak = 0.0;
    for (double v : cmap) cpeak = std::max(cpeak, v);
    if (cpeak <= 0) cpeak = 1.0;
    tf.caustic_nx = cnx; tf.caustic_nz = cnz;
    tf.caustic.resize(cmap.size());
    for (std::size_t i = 0; i < cmap.size(); ++i)
        tf.caustic[i] = static_cast<float>(cmap[i] / cpeak);

    return tf;
}

} // namespace celeris
