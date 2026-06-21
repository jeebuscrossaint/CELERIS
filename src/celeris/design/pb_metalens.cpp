#include "celeris/design/pb_metalens.hpp"

#include "celeris/design/polar_metalens.hpp"

#include <algorithm>
#include <cmath>
#include <future>
#include <thread>

namespace celeris {
namespace {

double wrap_pi(double a) {
    while (a > pi) a -= 2.0 * pi;
    while (a <= -pi) a += 2.0 * pi;
    return a;
}

// Circular-polarization decomposition of an output field [Ex; Ey] under RCP
// illumination. Basis: e_R = (x - i y)/sqrt2, e_L = (x + i y)/sqrt2, so the
// co (RCP) and cross (LCP) amplitudes are <e_R|E> and <e_L|E>.
//   a_co   = (Ex + i Ey) / sqrt2
//   a_cross= (Ex - i Ey) / sqrt2
struct Circular { cdouble co, cross; };
Circular circular_components(const JonesMatrix& J) {
    constexpr cdouble I{0.0, 1.0};
    const double inv = 1.0 / std::sqrt(2.0);
    // Input RCP: E_in = (1/sqrt2)[1; -i]; E_out = J E_in.
    cdouble ex = (J.xx - I * J.xy) * inv;
    cdouble ey = (J.yx - I * J.yy) * inv;
    return {(ex + I * ey) * inv, (ex - I * ey) * inv};
}

} // namespace

JonesMatrix solve_jones(const Material& incident, const Rcwa2DStack& stack,
                        const Material& substrate, double wavelength_um, int M) {
    auto rx = solve_rcwa_2d(incident, stack, substrate, wavelength_um, 0.0, 0.0,
                            1.0, 0.0, M, M);  // x-polarized input
    auto ry = solve_rcwa_2d(incident, stack, substrate, wavelength_um, 0.0, 0.0,
                            0.0, 1.0, M, M);  // y-polarized input
    // Columns are the response to each input polarization.
    return {rx.tx0, ry.tx0, rx.ty0, ry.ty0};
}

HwpAtom find_hwp_atom(const Material& pillar, const Material& background,
                      const Material& incident, const Material& substrate,
                      double period_um, double wavelength_um, double thickness_um,
                      double fill_min, double fill_max, int n_samples, int M) {
    // The (fill_x, fill_y) sweep already solves both polarizations per cell.
    PolarizationLibrary lib = build_polarization_library(
        pillar, background, incident, substrate, period_um, wavelength_um,
        thickness_um, fill_min, fill_max, n_samples, M);

    HwpAtom best;
    double best_eff = -1.0;
    for (const PolarCell& c : lib.cells) {
        cdouble tx = std::polar(c.amp_x, c.phase_x);
        cdouble ty = std::polar(c.amp_y, c.phase_y);
        // Spin-flip (RCP->LCP) conversion amplitude of an x/y birefringent atom
        // is (t_x - t_y)/2; efficiency is its modulus squared.
        double eff = std::norm(tx - ty) / 4.0;
        if (eff > best_eff) {
            best_eff = eff;
            best.fill_x = c.fill_x;
            best.fill_y = c.fill_y;
            best.thickness_um = thickness_um;
            best.t_x = tx;
            best.t_y = ty;
            best.retardance_deg = wrap_pi(c.phase_x - c.phase_y) * 180.0 / pi;
            best.conversion_efficiency = eff;
        }
    }
    return best;
}

std::vector<PbVerifyPoint> verify_pb_phase(
    const Material& pillar, const Material& background, const Material& incident,
    const Material& substrate, double period_um, double wavelength_um,
    const HwpAtom& atom, const std::vector<double>& rotations_rad, int M) {
    std::vector<PbVerifyPoint> out(rotations_rad.size());

    auto solve_one = [&](std::size_t i) {
        RectCell2D cell{pillar,           background, atom.fill_x,
                        atom.fill_y,      atom.thickness_um,
                        MetaShape::Rectangle, 0.5, rotations_rad[i]};
        Rcwa2DStack stack{period_um, period_um, {cell}};
        JonesMatrix J = solve_jones(incident, stack, substrate, wavelength_um, M);
        Circular cc = circular_components(J);
        out[i] = {rotations_rad[i] * 180.0 / pi, std::arg(cc.cross) * 180.0 / pi,
                  std::norm(cc.cross), std::norm(cc.co)};
    };

    unsigned hw = std::max(1u, std::thread::hardware_concurrency());
    int workers = std::min<int>((int)hw, std::max<int>(1, (int)rotations_rad.size()));
    std::vector<std::future<void>> jobs;
    for (int w = 0; w < workers; ++w)
        jobs.push_back(std::async(std::launch::async, [&, w] {
            for (std::size_t i = w; i < rotations_rad.size(); i += workers)
                solve_one(i);
        }));
    for (auto& j : jobs) j.get();
    return out;
}

PbMetalensDesign design_pb_metalens(const HwpAtom& atom, double period_um,
                                    double wavelength_um, const PhaseProfile& profile,
                                    double diameter_um, int handedness) {
    const int n = std::max(1, (int)std::round(diameter_um / period_um));
    const double center = (n - 1) / 2.0;

    PbMetalensDesign d;
    d.n_cells = n;
    d.period_um = period_um;
    d.atom = atom;
    d.rotation_rad.resize((std::size_t)n * n);
    d.t_cross.resize((std::size_t)n * n);

    // Geometric phase imprinted on the cross output is -handedness * 2*theta, so
    // theta = -handedness * phi / 2 hits any target phi(x,y) exactly.
    // The spin-flip amplitude of the un-rotated atom (uniform across the surface)
    // contributes only a global piston, irrelevant to the optical function -- the
    // only thing that matters is that the per-site phase tracks phi(x,y).
    cdouble a_cross0 = (atom.t_x - atom.t_y) / 2.0;
    const double piston = std::arg(a_cross0);

    double sq = 0;
    int count = 0;
    for (int iy = 0; iy < n; ++iy)
        for (int ix = 0; ix < n; ++ix) {
            double x = (ix - center) * period_um, y = (iy - center) * period_um;
            double phi = phase_profile_value(profile, x, y, wavelength_um);
            double theta = -handedness * phi / 2.0;
            std::size_t off = (std::size_t)iy * n + ix;
            d.rotation_rad[off] = theta;
            // Realized cross transmission = a_cross0 * exp(i * phi) (geometric).
            d.t_cross[off] = a_cross0 * std::polar(1.0, phi);
            // Phase error vs the target (modulo the global piston) -- identically
            // zero by construction: geometric phase is set by an angle, exactly.
            double err = wrap_pi(std::arg(d.t_cross[off]) - phi - piston);
            sq += err * err;
            ++count;
        }
    if (count) d.rms_phase_error_deg = std::sqrt(sq / count) * 180.0 / pi;
    d.conversion_efficiency = std::norm(a_cross0);
    return d;
}

PbMetalensDesign design_pb_metalens(const HwpAtom& atom, double period_um,
                                    double wavelength_um, double focal_length_um,
                                    double diameter_um, int handedness) {
    PhaseProfile p;
    p.kind = PhaseProfileKind::Focusing;
    p.focal_length_um = focal_length_um;
    return design_pb_metalens(atom, period_um, wavelength_um, p, diameter_um,
                              handedness);
}

} // namespace celeris
