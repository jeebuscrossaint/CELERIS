#include "celeris/design/achromatic.hpp"

#include "celeris/core.hpp"

#include <algorithm>
#include <cmath>
#include <future>
#include <thread>

namespace celeris {
namespace {

// Speed of light in vacuum, microns per femtosecond. With wavelengths in um this
// makes angular frequency omega = 2*pi*c/lambda come out in rad/fs, so a group
// delay d(phi)/d(omega) is in femtoseconds -- the natural unit for these designs.
constexpr double c_um_per_fs = 0.299792458;

// Smallest signed difference between two angles, in (-pi, pi].
double angle_diff(double a, double b) {
    double d = a - b;
    while (d > pi) d -= 2.0 * pi;
    while (d <= -pi) d += 2.0 * pi;
    return d;
}

// Unwrap a phase sequence in place (remove 2*pi jumps between consecutive samples)
// so it can be differentiated. The input order must be monotone in the variable
// the derivative is taken against.
void unwrap(std::vector<double>& p) {
    for (std::size_t i = 1; i < p.size(); ++i) {
        double d = p[i] - p[i - 1];
        while (d > pi) { p[i] -= 2.0 * pi; d -= 2.0 * pi; }
        while (d < -pi) { p[i] += 2.0 * pi; d += 2.0 * pi; }
    }
}

// Least-squares slope of y vs x.
double lstsq_slope(const std::vector<double>& x, const std::vector<double>& y) {
    const std::size_t n = x.size();
    if (n < 2) return 0.0;
    double mx = 0, my = 0;
    for (std::size_t i = 0; i < n; ++i) { mx += x[i]; my += y[i]; }
    mx /= n; my /= n;
    double num = 0, den = 0;
    for (std::size_t i = 0; i < n; ++i) {
        num += (x[i] - mx) * (y[i] - my);
        den += (x[i] - mx) * (x[i] - mx);
    }
    return den != 0.0 ? num / den : 0.0;
}

} // namespace

DispersiveLibrary build_dispersive_library(
    const Material& pillar, const Material& background, const Material& incident,
    const Material& substrate, double period_um,
    const std::vector<double>& band_wavelengths_um, double center_wavelength_um,
    double fill_min, double fill_max, int n_fills, double thick_lo, double thick_hi,
    int n_heights, int M) {
    DispersiveLibrary lib;
    lib.wavelengths_um = band_wavelengths_um;
    lib.center_wavelength_um = center_wavelength_um;
    lib.period_um = period_um;
    lib.n_fill = n_fills;
    lib.n_height = std::max(1, n_heights);
    const int nb = static_cast<int>(band_wavelengths_um.size());

    // Center sample = the band wavelength nearest the requested center.
    lib.center_index = 0;
    double best_d = std::abs(band_wavelengths_um[0] - center_wavelength_um);
    for (int j = 1; j < nb; ++j) {
        double d = std::abs(band_wavelengths_um[j] - center_wavelength_um);
        if (d < best_d) { best_d = d; lib.center_index = j; }
    }

    // Angular frequency at each band sample (rad/fs), for the group-delay slope.
    std::vector<double> omega(nb);
    for (int j = 0; j < nb; ++j)
        omega[j] = 2.0 * pi * c_um_per_fs / band_wavelengths_um[j];

    const int n_atoms = n_fills * lib.n_height;
    lib.atoms.resize(n_atoms);

    // Each (fill, height) atom is independent: solve it at every band wavelength,
    // then derive the center phase and the group delay. Parallel over atoms.
    auto solve_one = [&](int idx) {
        const int ih = idx / n_fills;
        const int iff = idx % n_fills;
        double f = fill_min + (fill_max - fill_min) * iff / (n_fills - 1);
        double h = lib.n_height == 1
                       ? thick_lo
                       : thick_lo + (thick_hi - thick_lo) * ih / (lib.n_height - 1);
        DispersiveAtom a;
        a.fill = f;
        a.thickness_um = h;
        a.phase_rad.resize(nb);
        a.amplitude.resize(nb);
        double amp_sum = 0.0;
        for (int j = 0; j < nb; ++j) {
            Rcwa2DStack cell{period_um, period_um,
                             {RectCell2D{pillar, background, f, f, h}}};
            auto r = solve_rcwa_2d(incident, cell, substrate, band_wavelengths_um[j],
                                   0.0, 0.0, /*Ex0=*/1.0, /*Ey0=*/0.0, M, M);
            a.phase_rad[j] = std::arg(r.tx0);
            a.amplitude[j] = std::abs(r.tx0);
            amp_sum += a.amplitude[j];
        }
        a.phase0_rad = a.phase_rad[lib.center_index];
        a.mean_amplitude = amp_sum / nb;
        // Group delay = slope of unwrapped phase vs omega (the dispersion lever).
        std::vector<double> ph = a.phase_rad;
        unwrap(ph);
        a.group_delay_fs = lstsq_slope(omega, ph);
        lib.atoms[idx] = std::move(a);
    };

    unsigned hw = std::max(1u, std::thread::hardware_concurrency());
    int workers = std::min<int>(static_cast<int>(hw), n_atoms);
    std::vector<std::future<void>> jobs;
    for (int w = 0; w < workers; ++w)
        jobs.push_back(std::async(std::launch::async, [&, w] {
            for (int i = w; i < n_atoms; i += workers) solve_one(i);
        }));
    for (auto& j : jobs) j.get();

    lib.gd_min_fs = lib.gd_max_fs = lib.atoms.empty() ? 0.0 : lib.atoms[0].group_delay_fs;
    for (const auto& a : lib.atoms) {
        lib.gd_min_fs = std::min(lib.gd_min_fs, a.group_delay_fs);
        lib.gd_max_fs = std::max(lib.gd_max_fs, a.group_delay_fs);
    }
    return lib;
}

AchromaticDesign design_achromatic_metalens(const DispersiveLibrary& lib,
                                            double focal_length_um,
                                            double diameter_um, double gd_weight,
                                            double amplitude_weight) {
    const double p = lib.period_um;
    const double lambda0 = lib.center_wavelength_um;
    const double k0 = 2.0 * pi / lambda0;
    const double f = focal_length_um;
    const double R = diameter_um / 2.0;
    const int n = std::max(1, static_cast<int>(std::round(diameter_um / p)));
    const double center = (n - 1) / 2.0;

    AchromaticDesign d;
    d.n_cells = n;
    d.period_um = p;
    d.center_wavelength_um = lambda0;
    d.focal_length_um = f;
    d.diameter_um = diameter_um;
    d.atom_index.resize(static_cast<std::size_t>(n) * n);
    d.fill_map.resize(static_cast<std::size_t>(n) * n);
    d.thickness_map.resize(static_cast<std::size_t>(n) * n);

    // The required group-delay span over the aperture (the achromatic budget) and
    // the offset that centers the target range within what the library can supply
    // -- exploiting the freedom to add a radius-independent reference group delay.
    const double required_span = (std::sqrt(R * R + f * f) - f) / c_um_per_fs;
    const double lib_mid = 0.5 * (lib.gd_min_fs + lib.gd_max_fs);
    const double gd_ref = lib_mid + 0.5 * required_span;  // raw target midpoint = -span/2

    // Convert a group-delay error to the phase error it causes at the band EDGE
    // (omega - omega0 at the extreme), which is what actually blurs focus there.
    double domega_half = 0.0;
    if (lib.wavelengths_um.size() >= 2) {
        double w_lo = 2.0 * pi * c_um_per_fs / lib.wavelengths_um.front();
        double w_hi = 2.0 * pi * c_um_per_fs / lib.wavelengths_um.back();
        domega_half = 0.5 * std::abs(w_lo - w_hi);
    }

    double sq_phase = 0.0, sq_gd = 0.0, amp_sum = 0.0;
    int count = 0;
    bool any = false;
    double hmin = 0, hmax = 0;
    for (int iy = 0; iy < n; ++iy) {
        for (int ix = 0; ix < n; ++ix) {
            const double x = (ix - center) * p, y = (iy - center) * p;
            const double r = std::sqrt(x * x + y * y);
            const double sag = std::sqrt(r * r + f * f) - f;
            const double target_phase = -k0 * sag;                 // base phase (mod 2pi)
            const double target_gd = -sag / c_um_per_fs + gd_ref;  // group delay [fs]

            // Two-objective atom pick. The group-delay error is mapped to the phase
            // error it causes at the band edge (dgd * domega_half), so both terms
            // are in radians and gd_weight ~ 1 weighs center vs band-edge focus.
            int best = 0;
            double best_cost = 1e300;
            for (int q = 0; q < static_cast<int>(lib.atoms.size()); ++q) {
                const DispersiveAtom& a = lib.atoms[q];
                double dphi = angle_diff(a.phase0_rad, target_phase);
                double dgd = (a.group_delay_fs - target_gd) * domega_half;  // -> radians
                double damp = 1.0 - a.mean_amplitude;
                double cost = dphi * dphi + gd_weight * dgd * dgd +
                              amplitude_weight * damp * damp;
                if (cost < best_cost) { best_cost = cost; best = q; }
            }

            const DispersiveAtom& a = lib.atoms[best];
            std::size_t off = static_cast<std::size_t>(iy) * n + ix;
            d.atom_index[off] = best;
            d.fill_map[off] = a.fill;
            d.thickness_map[off] = a.thickness_um;

            double ephi = angle_diff(a.phase0_rad, target_phase);
            sq_phase += ephi * ephi;
            double egd = a.group_delay_fs - target_gd;
            sq_gd += egd * egd;
            amp_sum += a.mean_amplitude;
            ++count;
            if (!any) { hmin = hmax = a.thickness_um; any = true; }
            else { hmin = std::min(hmin, a.thickness_um); hmax = std::max(hmax, a.thickness_um); }
        }
    }

    d.rms_phase_error_deg = std::sqrt(sq_phase / count) * 180.0 / pi;
    d.rms_group_delay_error_fs = std::sqrt(sq_gd / count);
    d.mean_amplitude = amp_sum / count;
    d.required_gd_span_fs = required_span;
    d.available_gd_span_fs = lib.gd_max_fs - lib.gd_min_fs;
    d.gd_coverage = required_span > 1e-12 ? d.available_gd_span_fs / required_span : 1.0;
    d.min_height_um = hmin;
    d.max_height_um = hmax;
    d.single_height = (hmax - hmin) < 1e-9;
    return d;
}

std::vector<AchromaticFocalPoint> verify_achromatic_focus(
    const DispersiveLibrary& lib, const AchromaticDesign& design) {
    const double p = design.period_um;
    const double center = (design.n_cells - 1) / 2.0;
    const double R_ap = design.diameter_um / 2.0;
    const int nb = static_cast<int>(lib.wavelengths_um.size());

    // Aperture cell geometry + the chosen atom index at each site.
    struct Cell { double x, y; int atom; };
    std::vector<Cell> cells;
    for (int iy = 0; iy < design.n_cells; ++iy)
        for (int ix = 0; ix < design.n_cells; ++ix) {
            const double x = (ix - center) * p, y = (iy - center) * p;
            if (std::sqrt(x * x + y * y) > R_ap) continue;
            cells.push_back(
                {x, y, design.atom_index[static_cast<std::size_t>(iy) * design.n_cells + ix]});
        }

    auto on_axis = [&](double z, double k, int band_j) {
        cdouble E{0.0, 0.0};
        for (const Cell& c : cells) {
            const DispersiveAtom& a = lib.atoms[c.atom];
            cdouble t = std::polar(a.amplitude[band_j], a.phase_rad[band_j]);
            double Rr = std::sqrt(c.x * c.x + c.y * c.y + z * z);
            E += t * std::polar(1.0 / Rr, k * Rr);
        }
        return std::norm(E);
    };

    // Reference peak at the center wavelength / design focus.
    double k0 = 2.0 * pi / lib.center_wavelength_um;
    double design_peak = on_axis(design.focal_length_um, k0, lib.center_index);

    std::vector<AchromaticFocalPoint> out;
    for (int j = 0; j < nb; ++j) {
        double lam = lib.wavelengths_um[j];
        double k = 2.0 * pi / lam;
        // Search a generous z window around the design focus (achromatic designs
        // intentionally keep f roughly fixed, so center on the target, not on the
        // chromatic estimate f0*l0/l).
        double zlo = 0.4 * design.focal_length_um, zhi = 1.8 * design.focal_length_um;
        const int NZ = 200;
        double best_z = design.focal_length_um, best_I = -1.0;
        for (int q = 0; q < NZ; ++q) {
            double z = zlo + (zhi - zlo) * q / (NZ - 1);
            double I = on_axis(z, k, j);
            if (I > best_I) { best_I = I; best_z = z; }
        }
        out.push_back({lam, best_z, design_peak > 0 ? best_I / design_peak : 0.0});
    }
    return out;
}

MetalensDesign to_metalens_design(const AchromaticDesign& d) {
    MetalensDesign m;
    m.n_cells = d.n_cells;
    m.period_um = d.period_um;
    m.fill_map = d.fill_map;
    m.rms_phase_error_deg = d.rms_phase_error_deg;
    m.mean_amplitude = d.mean_amplitude;
    return m;
}

} // namespace celeris
