#pragma once
/// @file
/// @brief Phase unwrapping across a wavelength band, with adaptive refinement.
///
/// The dispersive libraries (`achromatic.hpp`, `pb_achromatic.hpp`) take an
/// atom's group delay from the slope of its unwrapped phase vs angular frequency.
/// Plain unwrapping forces each step between adjacent band samples into
/// (−π, π]. A real step larger than π (a large group delay on a coarse band, or a
/// resonance between samples) is aliased by 2π, and the fitted group delay is
/// then wrong by 2π/Δω (issue #21). This helper solves extra samples at the
/// angular-frequency midpoint of any interval whose wrapped step is larger than a
/// threshold, so the 2π branch at every band sample is decided on a finer grid.
/// The refined samples are used only to unwrap and to find the smallest
/// amplitude; the caller still fits the group delay on the band samples.

#include "celeris/core.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

namespace celeris::detail {

/// @brief Result of `unwrap_band_refined`.
struct BandPhase {
    std::vector<double> unwrapped;  ///< Unwrapped phase at each band sample (rad).
    double min_amplitude = 0;       ///< Smallest |response| over band AND refined samples.
    int n_refined = 0;              ///< Extra samples solved during refinement.
};

/// @brief Unwrap the phase of a complex band response, refining coarse intervals.
///
/// For each pair of adjacent band samples, if the wrapped phase step exceeds
/// `max_step_rad` the response is sampled at the angular-frequency midpoint and
/// both halves are checked again, up to `max_depth` levels (at most
/// 2^max_depth − 1 extra samples per interval). The unwrapped increments of the
/// halves are summed, so the band samples get the 2π branch the finer grid sees.
/// A narrow feature that falls between samples and never produces a large step
/// (a sharp resonance or a band zero) can still be missed; only denser band
/// sampling finds those.
/// @param band_um      Band wavelengths (µm), ascending.
/// @param at_band      Response at each band sample (already solved).
/// @param sample       Callable `cdouble(double lambda_um)` that solves one extra sample.
/// @param max_step_rad Wrapped step above which an interval is refined (rad).
/// @param max_depth    Maximum bisection depth per interval.
template <class Sample>
BandPhase unwrap_band_refined(const std::vector<double>& band_um,
                              const std::vector<cdouble>& at_band, Sample&& sample,
                              double max_step_rad = pi / 2, int max_depth = 4) {
    constexpr double c_um_per_fs = 0.299792458;
    auto wrap = [](double d) {
        while (d > pi) d -= 2.0 * pi;
        while (d <= -pi) d += 2.0 * pi;
        return d;
    };
    BandPhase out;
    const std::size_t n = at_band.size();
    out.unwrapped.resize(n);
    if (n == 0) return out;
    out.min_amplitude = std::abs(at_band[0]);
    for (const auto& v : at_band) out.min_amplitude = std::min(out.min_amplitude, std::abs(v));

    // Unwrapped phase increment from (wa, va) to (wb, vb), bisecting in omega.
    auto increment = [&](auto&& self, double wa, cdouble va, double wb, cdouble vb,
                         int depth) -> double {
        const double d = wrap(std::arg(vb) - std::arg(va));
        if (std::abs(d) <= max_step_rad || depth <= 0) return d;
        const double wm = 0.5 * (wa + wb);
        const cdouble vm = sample(2.0 * pi * c_um_per_fs / wm);
        ++out.n_refined;
        out.min_amplitude = std::min(out.min_amplitude, std::abs(vm));
        return self(self, wa, va, wm, vm, depth - 1) + self(self, wm, vm, wb, vb, depth - 1);
    };

    out.unwrapped[0] = std::arg(at_band[0]);
    for (std::size_t j = 1; j < n; ++j) {
        const double wa = 2.0 * pi * c_um_per_fs / band_um[j - 1];
        const double wb = 2.0 * pi * c_um_per_fs / band_um[j];
        out.unwrapped[j] = out.unwrapped[j - 1] +
                           increment(increment, wa, at_band[j - 1], wb, at_band[j], max_depth);
    }
    return out;
}

} // namespace celeris::detail
