#include "celeris/design/pb_achromatic.hpp"

#include <algorithm>
#include <cmath>
#include <future>
#include <thread>

namespace celeris {
namespace {

// Speed of light, microns per femtosecond -- with wavelengths in um, the angular
// frequency omega = 2*pi*c/lambda comes out in rad/fs and a group delay d(phi)/
// d(omega) is in femtoseconds (see achromatic.cpp for the same convention).
constexpr double c_um_per_fs = 0.299792458;

double angle_diff(double a, double b) {
    double d = a - b;
    while (d > pi) d -= 2.0 * pi;
    while (d <= -pi) d += 2.0 * pi;
    return d;
}

void unwrap(std::vector<double>& p) {
    for (std::size_t i = 1; i < p.size(); ++i) {
        double d = p[i] - p[i - 1];
        while (d > pi) { p[i] -= 2.0 * pi; d -= 2.0 * pi; }
        while (d < -pi) { p[i] += 2.0 * pi; d += 2.0 * pi; }
    }
}

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

DispersivePbLibrary build_dispersive_pb_library(
    const Material& pillar, const Material& background, const Material& incident,
    const Material& substrate, double period_um,
    const std::vector<double>& band_wavelengths_um, double center_wavelength_um,
    double fill_min, double fill_max, int n_fills, double thickness_um, int M,
    double min_amplitude) {
    DispersivePbLibrary lib;
    lib.wavelengths_um = band_wavelengths_um;
    lib.center_wavelength_um = center_wavelength_um;
    lib.period_um = period_um;
    lib.thickness_um = thickness_um;
    lib.n_fill = n_fills;
    const int nb = static_cast<int>(band_wavelengths_um.size());

    // Center sample = the band wavelength nearest the requested center.
    lib.center_index = 0;
    double best_d = std::abs(band_wavelengths_um[0] - center_wavelength_um);
    for (int j = 1; j < nb; ++j) {
        double d = std::abs(band_wavelengths_um[j] - center_wavelength_um);
        if (d < best_d) { best_d = d; lib.center_index = j; }
    }

    std::vector<double> omega(nb);
    for (int j = 0; j < nb; ++j)
        omega[j] = 2.0 * pi * c_um_per_fs / band_wavelengths_um[j];

    // (fill_x, fill_y) grid -- birefringence comes from fill_x != fill_y, and the
    // size/aspect ratio sets the group delay. Both halves of the grid (fy<fx and
    // fy>fx) are kept: they give a_cross of opposite sign but distinct dispersion.
    struct Spec { double fx, fy; };
    std::vector<Spec> specs;
    specs.reserve(static_cast<std::size_t>(n_fills) * n_fills);
    auto fill_at = [&](int i) {
        return n_fills <= 1 ? fill_min
                            : fill_min + (fill_max - fill_min) * i / (n_fills - 1);
    };
    for (int iy = 0; iy < n_fills; ++iy)
        for (int ix = 0; ix < n_fills; ++ix)
            specs.push_back({fill_at(ix), fill_at(iy)});

    const int n_atoms = static_cast<int>(specs.size());
    lib.atoms.resize(n_atoms);

    // Each atom: two RCWA solves (x- and y-pol) per band wavelength give the
    // diagonal Jones transmissions t_x, t_y (an axis-aligned rectangle has no
    // cross terms), then a_cross = (t_x - t_y)/2 is the spin-flip amplitude.
    auto solve_one = [&](int idx) {
        const Spec& s = specs[idx];
        DispersivePbAtom a;
        a.fill_x = s.fx;
        a.fill_y = s.fy;
        a.thickness_um = thickness_um;
        a.cross.resize(nb);
        std::vector<double> ph(nb);
        double amp_sum = 0.0, conv_sum = 0.0;
        for (int j = 0; j < nb; ++j) {
            RectCell2D cell{pillar, background, s.fx, s.fy, thickness_um,
                            MetaShape::Rectangle, 0.5};
            Rcwa2DStack stack{period_um, period_um, {cell}};
            auto rx = solve_rcwa_2d(incident, stack, substrate, band_wavelengths_um[j],
                                    0.0, 0.0, /*Ex0=*/1.0, /*Ey0=*/0.0, M, M);
            auto ry = solve_rcwa_2d(incident, stack, substrate, band_wavelengths_um[j],
                                    0.0, 0.0, /*Ex0=*/0.0, /*Ey0=*/1.0, M, M);
            cdouble tx = rx.tx0, ty = ry.ty0;
            cdouble acr = (tx - ty) / 2.0;
            a.cross[j] = acr;
            ph[j] = std::arg(acr);
            amp_sum += std::abs(acr);
            conv_sum += std::norm(acr);
            if (j == lib.center_index)
                a.retardance_center_deg =
                    angle_diff(std::arg(tx), std::arg(ty)) * 180.0 / pi;
        }
        a.phase0_rad = ph[lib.center_index];
        a.mean_amplitude = amp_sum / nb;
        a.mean_conversion = conv_sum / nb;
        std::vector<double> u = ph;
        unwrap(u);
        a.group_delay_fs = lstsq_slope(omega, u);
        lib.atoms[idx] = std::move(a);
    };

    unsigned hw = std::max(1u, std::thread::hardware_concurrency());
    int workers = std::min<int>(static_cast<int>(hw), std::max(1, n_atoms));
    std::vector<std::future<void>> jobs;
    for (int w = 0; w < workers; ++w)
        jobs.push_back(std::async(std::launch::async, [&, w] {
            for (int i = w; i < n_atoms; i += workers) solve_one(i);
        }));
    for (auto& j : jobs) j.get();

    // Keep only genuinely birefringent atoms: a near-square pillar has a_cross ~ 0,
    // so its phase (and group delay) is meaningless noise that could otherwise win
    // the group-delay selection with zero spin-flip amplitude. Drop atoms below the
    // mean-|a_cross| floor (but never empty the library).
    std::vector<DispersivePbAtom> kept;
    kept.reserve(lib.atoms.size());
    for (auto& a : lib.atoms)
        if (a.mean_amplitude >= min_amplitude) kept.push_back(std::move(a));
    if (!kept.empty()) lib.atoms = std::move(kept);

    lib.gd_min_fs = lib.gd_max_fs = lib.atoms.empty() ? 0.0 : lib.atoms[0].group_delay_fs;
    for (const auto& a : lib.atoms) {
        lib.gd_min_fs = std::min(lib.gd_min_fs, a.group_delay_fs);
        lib.gd_max_fs = std::max(lib.gd_max_fs, a.group_delay_fs);
    }
    return lib;
}

PbAchromaticDesign design_pb_achromatic_metalens(
    const DispersivePbLibrary& lib, double focal_length_um, double diameter_um,
    int handedness, double gd_weight, double amplitude_weight) {
    const double p = lib.period_um;
    const double lambda0 = lib.center_wavelength_um;
    const double k0 = 2.0 * pi / lambda0;
    const double f = focal_length_um;
    const double R = diameter_um / 2.0;
    const int n = std::max(1, static_cast<int>(std::round(diameter_um / p)));
    const double center = (n - 1) / 2.0;
    handedness = handedness >= 0 ? 1 : -1;

    PbAchromaticDesign d;
    d.n_cells = n;
    d.period_um = p;
    d.thickness_um = lib.thickness_um;
    d.handedness = handedness;
    d.center_wavelength_um = lambda0;
    d.focal_length_um = f;
    d.diameter_um = diameter_um;
    d.atom_index.resize(static_cast<std::size_t>(n) * n);
    d.rotation_rad.resize(static_cast<std::size_t>(n) * n);
    d.fill_x_map.resize(static_cast<std::size_t>(n) * n);
    d.fill_y_map.resize(static_cast<std::size_t>(n) * n);

    // Required group-delay span over the aperture (the achromatic budget) and the
    // radius-independent reference that centers the target range within the library
    // (the freedom to add a constant reference group delay; see achromatic.cpp).
    const double required_span = (std::sqrt(R * R + f * f) - f) / c_um_per_fs;
    const double lib_mid = 0.5 * (lib.gd_min_fs + lib.gd_max_fs);
    const double gd_ref = lib_mid + 0.5 * required_span;

    // Map a group-delay error to the band-edge phase error it causes, so the GD
    // term is in radians like the amplitude term.
    double domega_half = 0.0;
    if (lib.wavelengths_um.size() >= 2) {
        double w_lo = 2.0 * pi * c_um_per_fs / lib.wavelengths_um.front();
        double w_hi = 2.0 * pi * c_um_per_fs / lib.wavelengths_um.back();
        domega_half = 0.5 * std::abs(w_lo - w_hi);
    }

    double sq_phase = 0.0, sq_gd = 0.0, amp_sum = 0.0, conv_sum = 0.0;
    int count = 0;
    for (int iy = 0; iy < n; ++iy) {
        for (int ix = 0; ix < n; ++ix) {
            const double x = (ix - center) * p, y = (iy - center) * p;
            const double r = std::sqrt(x * x + y * y);
            const double sag = std::sqrt(r * r + f * f) - f;
            const double target_phase = -k0 * sag;                 // base phase (mod 2pi)
            const double target_gd = -sag / c_um_per_fs + gd_ref;  // group delay [fs]

            // Atom pick: the base phase is free (the rotation will hit it exactly
            // for ANY atom), so the only objectives are the group-delay match and
            // the spin-flip conversion. gd_weight=0 -> best-conversion atom only.
            int best = 0;
            double best_cost = 1e300;
            for (int q = 0; q < static_cast<int>(lib.atoms.size()); ++q) {
                const DispersivePbAtom& a = lib.atoms[q];
                double dgd = (a.group_delay_fs - target_gd) * domega_half;  // -> radians
                double damp = 1.0 - a.mean_amplitude;
                double cost = gd_weight * dgd * dgd + amplitude_weight * damp * damp;
                if (cost < best_cost) { best_cost = cost; best = q; }
            }

            const DispersivePbAtom& a = lib.atoms[best];
            // Geometric phase phi_geo must carry (target_phase - atom base phase) so
            // the realized base phase = atom_phase0 + phi_geo = target_phase exactly.
            // phi_geo = -2*handedness*theta  =>  theta = -handedness*phi_geo/2.
            const double phi_geo = target_phase - a.phase0_rad;
            const double theta = -handedness * phi_geo / 2.0;

            std::size_t off = static_cast<std::size_t>(iy) * n + ix;
            d.atom_index[off] = best;
            d.rotation_rad[off] = theta;
            d.fill_x_map[off] = a.fill_x;
            d.fill_y_map[off] = a.fill_y;

            // Realized base phase at center: atom_phase0 + phi_geo == target_phase,
            // so this residual is identically zero (geometric phase is exact).
            double ephi = angle_diff(a.phase0_rad + phi_geo, target_phase);
            sq_phase += ephi * ephi;
            double egd = a.group_delay_fs - target_gd;
            sq_gd += egd * egd;
            amp_sum += a.mean_amplitude;
            conv_sum += a.mean_conversion;
            ++count;
        }
    }

    d.rms_phase_error_deg = std::sqrt(sq_phase / count) * 180.0 / pi;
    d.rms_group_delay_error_fs = std::sqrt(sq_gd / count);
    d.mean_amplitude = amp_sum / count;
    d.mean_conversion = conv_sum / count;
    d.required_gd_span_fs = required_span;
    d.available_gd_span_fs = lib.gd_max_fs - lib.gd_min_fs;
    d.gd_coverage = required_span > 1e-12 ? d.available_gd_span_fs / required_span : 1.0;
    return d;
}

std::vector<AchromaticFocalPoint> verify_pb_achromatic_focus(
    const DispersivePbLibrary& lib, const PbAchromaticDesign& design) {
    const double p = design.period_um;
    const double center = (design.n_cells - 1) / 2.0;
    const double R_ap = design.diameter_um / 2.0;
    const int nb = static_cast<int>(lib.wavelengths_um.size());

    // Aperture cells: chosen atom index + the geometric phase the rotation imprints
    // (wavelength-independent), realized t_cross(omega) = a_cross(omega)*exp(i*phi_geo).
    struct Cell { double x, y; int atom; double phi_geo; };
    std::vector<Cell> cells;
    for (int iy = 0; iy < design.n_cells; ++iy)
        for (int ix = 0; ix < design.n_cells; ++ix) {
            const double x = (ix - center) * p, y = (iy - center) * p;
            if (std::sqrt(x * x + y * y) > R_ap) continue;
            std::size_t off = static_cast<std::size_t>(iy) * design.n_cells + ix;
            double phi_geo = -2.0 * design.handedness * design.rotation_rad[off];
            cells.push_back({x, y, design.atom_index[off], phi_geo});
        }

    auto on_axis = [&](double z, double k, int band_j) {
        cdouble E{0.0, 0.0};
        for (const Cell& c : cells) {
            cdouble t = lib.atoms[c.atom].cross[band_j] * std::polar(1.0, c.phi_geo);
            double Rr = std::sqrt(c.x * c.x + c.y * c.y + z * z);
            E += t * std::polar(1.0 / Rr, k * Rr);
        }
        return std::norm(E);
    };

    double k0 = 2.0 * pi / lib.center_wavelength_um;
    double design_peak = on_axis(design.focal_length_um, k0, lib.center_index);

    std::vector<AchromaticFocalPoint> out;
    for (int j = 0; j < nb; ++j) {
        double lam = lib.wavelengths_um[j];
        double k = 2.0 * pi / lam;
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

} // namespace celeris
