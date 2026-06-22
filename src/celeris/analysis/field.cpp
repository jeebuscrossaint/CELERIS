#include "celeris/analysis/field.hpp"

#include "celeris/analysis/focal.hpp"  // propagate_pillars (GPU/CPU far-field)
#include "celeris/core.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

namespace celeris {

std::vector<FieldPoint> analyze_field_of_view(const MetalensDesign& lens,
                                              const UnitCellLibrary& lib,
                                              double focal_length_um,
                                              double wavelength_um,
                                              double diameter_um,
                                              const std::vector<double>& angles_deg) {
    const double p = lens.period_um;
    const double center = (lens.n_cells - 1) / 2.0;
    const double R_ap = diameter_um / 2.0;
    const double k = 2.0 * pi / wavelength_um;

    struct Cell { double x, y; cdouble t; };
    std::vector<Cell> cells;
    for (int iy = 0; iy < lens.n_cells; ++iy)
        for (int ix = 0; ix < lens.n_cells; ++ix) {
            const double x = (ix - center) * p, y = (iy - center) * p;
            if (std::sqrt(x * x + y * y) > R_ap) continue;
            const double fill = lens.fill_map[static_cast<std::size_t>(iy) * lens.n_cells + ix];
            cells.push_back({x, y, lib.transmission_for_fill(fill)});
        }

    // Focal-plane field for an incident tilt theta (in x-z): the illumination
    // adds a linear phase exp(i k sin(theta) x) across the aperture.
    auto peak_near = [&](double theta, double xc) {
        const double s = std::sin(theta);
        const double dl = wavelength_um * focal_length_um / diameter_um;
        // The spot lands near f*tan(theta) for a hyperbolic lens (chief ray) but
        // near f*sin(theta) for a recentered parabolic (quadratic) lens; at large
        // angles those differ by more than a few spot widths. Widen the half-window
        // to span both predictions so the search captures the peak for EITHER lens
        // type (at small angles f*tan ~ f*sin, so this collapses to the old window).
        const double x_sin = focal_length_um * s;
        const double x_tan = focal_length_um * std::tan(theta);
        const double W = std::max(3.0 * dl, 2.0) + std::abs(x_tan - x_sin);
        const int NG = 81;
        double best = 0.0, best_x = xc;
        for (int j = 0; j < NG; ++j) {
            double fy = -W + 2.0 * W * j / (NG - 1);
            for (int i = 0; i < NG; ++i) {
                double fx = xc - W + 2.0 * W * i / (NG - 1);
                cdouble E{0.0, 0.0};
                for (const Cell& c : cells) {
                    double R = std::sqrt((fx - c.x) * (fx - c.x) +
                                         (fy - c.y) * (fy - c.y) +
                                         focal_length_um * focal_length_um);
                    E += c.t * std::polar(1.0, k * s * c.x) * std::polar(1.0 / R, k * R);
                }
                double I = std::norm(E);
                if (I > best) { best = I; best_x = fx; }
            }
        }
        return std::pair<double, double>{best, best_x};
    };

    double on_axis_peak = peak_near(0.0, 0.0).first;

    std::vector<FieldPoint> out;
    for (double ang : angles_deg) {
        double theta = ang * pi / 180.0;
        double xc = focal_length_um * std::tan(theta);  // paraxial spot position
        auto [pk, px] = peak_near(theta, xc);
        out.push_back({ang, on_axis_peak > 0 ? pk / on_axis_peak : 0.0, px});
    }
    return out;
}

std::vector<FieldPoint> analyze_wide_fov(const MetalensDesign& lens,
                                         const UnitCellLibrary& lib,
                                         double focal_length_um,
                                         double wavelength_um,
                                         double lens_diameter_um,
                                         double stop_diameter_um,
                                         double stop_distance_um,
                                         const std::vector<double>& angles_deg) {
    const double p = lens.period_um;
    const double center = (lens.n_cells - 1) / 2.0;
    const double R_lens = lens_diameter_um / 2.0;
    const double R_stop = stop_diameter_um / 2.0;
    const double k = 2.0 * pi / wavelength_um;
    // The stop is the limiting aperture, so resolution (spot scale) ~ lambda*f/stop.
    const double dl_stop = wavelength_um * focal_length_um / stop_diameter_um;

    struct Cell { double x, y; cdouble t; };
    std::vector<Cell> cells;
    for (int iy = 0; iy < lens.n_cells; ++iy)
        for (int ix = 0; ix < lens.n_cells; ++ix) {
            const double x = (ix - center) * p, y = (iy - center) * p;
            if (std::sqrt(x * x + y * y) > R_lens) continue;
            const double fill = lens.fill_map[static_cast<std::size_t>(iy) * lens.n_cells + ix];
            cells.push_back({x, y, lib.transmission_for_fill(fill)});
        }

    // Peak focal intensity for a field angle theta. Only the cells under the
    // DECENTERED stop patch (disk of radius R_stop centered at x_patch) are
    // illuminated; the tilt adds exp(i k sin(theta) x) across that patch.
    auto peak = [&](double theta, double x_patch) {
        const double s = std::sin(theta);
        // Search window: the spot lands between f*sin(theta) (recentered parabola)
        // and f*tan(theta) (chief ray / distortion); cover both with margin.
        const double x_sin = focal_length_um * s;
        const double x_tan = focal_length_um * std::tan(theta);
        const double xc = 0.5 * (x_sin + x_tan);
        const double W = std::max(3.0 * dl_stop, 2.0) + 0.5 * std::abs(x_tan - x_sin);
        const int NG = 81;
        const double R_stop2 = R_stop * R_stop;
        double best = 0.0, best_x = xc;
        for (int j = 0; j < NG; ++j) {
            double fy = -W + 2.0 * W * j / (NG - 1);
            for (int i = 0; i < NG; ++i) {
                double fx = xc - W + 2.0 * W * i / (NG - 1);
                cdouble E{0.0, 0.0};
                for (const Cell& c : cells) {
                    const double dx = c.x - x_patch, dy = c.y;
                    if (dx * dx + dy * dy > R_stop2) continue;  // stop patch mask
                    double R = std::sqrt((fx - c.x) * (fx - c.x) +
                                         (fy - c.y) * (fy - c.y) +
                                         focal_length_um * focal_length_um);
                    E += c.t * std::polar(1.0, k * s * c.x) * std::polar(1.0 / R, k * R);
                }
                double I = std::norm(E);
                if (I > best) { best = I; best_x = fx; }
            }
        }
        return std::pair<double, double>{best, best_x};
    };

    double on_axis_peak = peak(0.0, 0.0).first;

    std::vector<FieldPoint> out;
    for (double ang : angles_deg) {
        double theta = ang * pi / 180.0;
        double x_patch = stop_distance_um * std::tan(theta);  // walk-off on the lens
        auto [pk, px] = peak(theta, x_patch);
        out.push_back({ang, on_axis_peak > 0 ? pk / on_axis_peak : 0.0, px});
    }
    return out;
}

FieldGrid analyze_field_grid(const MetalensDesign& lens,
                             const UnitCellLibrary& lib,
                             double focal_length_um, double wavelength_um,
                             double diameter_um, double max_angle_deg,
                             int n_half, int psf_n) {
    const double p = lens.period_um;
    const double center = (lens.n_cells - 1) / 2.0;
    const double R_ap = diameter_um / 2.0;
    const double k = 2.0 * pi / wavelength_um;

    // The aperture cells (positions + designed transmission), built once; each
    // field point reuses them with its own tilt phase.
    std::vector<double> px, py;
    std::vector<cdouble> t0;
    for (int iy = 0; iy < lens.n_cells; ++iy)
        for (int ix = 0; ix < lens.n_cells; ++ix) {
            const double x = (ix - center) * p, y = (iy - center) * p;
            if (std::sqrt(x * x + y * y) > R_ap) continue;
            const double fill = lens.fill_map[static_cast<std::size_t>(iy) * lens.n_cells + ix];
            px.push_back(x);
            py.push_back(y);
            t0.push_back(lib.transmission_for_fill(fill));
        }
    const std::size_t npil = px.size();

    // Window: a few diffraction-limited spot widths, large enough to hold the
    // coma flare that grows off-axis. Sampling resolves the spot for FWHM.
    const double dl = wavelength_um * focal_length_um / diameter_um;
    const double W = std::max(6.0 * dl, 4.0);
    const double step = 2.0 * W / (psf_n - 1);

    // Peak intensity + tangential/sagittal FWHM of one tilted-illumination PSF.
    struct Spot { double peak, fwhm_x, fwhm_y; };
    std::vector<cdouble> tt(npil);
    auto spot_at = [&](double tx, double ty, double cx, double cy) -> Spot {
        const double sx = std::sin(tx), sy = std::sin(ty);
        for (std::size_t c = 0; c < npil; ++c)
            tt[c] = t0[c] * std::polar(1.0, k * (sx * px[c] + sy * py[c]));
        PsfMap m = propagate_pillars(px, py, tt, cx, cy, focal_length_um,
                                     wavelength_um, psf_n, W);
        // Peak pixel.
        int pix = 0, piy = 0;
        double best = -1.0;
        for (int j = 0; j < psf_n; ++j)
            for (int i = 0; i < psf_n; ++i) {
                double v = m.intensity[static_cast<std::size_t>(j) * psf_n + i];
                if (v > best) { best = v; pix = i; piy = j; }
            }
        const double half = best / 2.0;
        auto fwhm_along = [&](bool horizontal) {
            int lo = -1, hi = -1;
            for (int s = 0; s < psf_n; ++s) {
                std::size_t idx = horizontal
                    ? static_cast<std::size_t>(piy) * psf_n + s
                    : static_cast<std::size_t>(s) * psf_n + pix;
                if (m.intensity[idx] >= half) { if (lo < 0) lo = s; hi = s; }
            }
            return (hi > lo && lo >= 0) ? (hi - lo) * step : 0.0;
        };
        return {best, fwhm_along(true), fwhm_along(false)};
    };

    Spot on_axis = spot_at(0.0, 0.0, 0.0, 0.0);
    const double ref = on_axis.peak > 0 ? on_axis.peak : 1.0;

    FieldGrid g;
    g.n = 2 * n_half + 1;
    g.max_angle_deg = max_angle_deg;
    g.points.reserve(static_cast<std::size_t>(g.n) * g.n);
    for (int jy = -n_half; jy <= n_half; ++jy) {
        const double ang_y = n_half > 0 ? max_angle_deg * jy / n_half : 0.0;
        const double ty = ang_y * pi / 180.0;
        for (int jx = -n_half; jx <= n_half; ++jx) {
            const double ang_x = n_half > 0 ? max_angle_deg * jx / n_half : 0.0;
            const double tx = ang_x * pi / 180.0;
            const double cx = focal_length_um * std::tan(tx);
            const double cy = focal_length_um * std::tan(ty);
            Spot s = spot_at(tx, ty, cx, cy);
            g.points.push_back({ang_x, ang_y, s.peak / ref, s.fwhm_x, s.fwhm_y,
                                cx, cy});
        }
    }
    return g;
}

} // namespace celeris
