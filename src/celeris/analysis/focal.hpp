#pragma once
// Focal-plane analysis — does the designed lens actually focus, and how well?
//
// Takes a finished metalens design, builds its complex aperture field (each
// pillar's transmission from the library), propagates it to the focal plane by
// Rayleigh-Sommerfeld summation, and reports the standard optical metrics:
//   - Strehl ratio: peak intensity vs an ideal (perfect-phase) lens (1.0 = perfect)
//   - focal-spot FWHM, compared to the diffraction limit ~ lambda*f/D
//   - encircled energy within the first Airy null (focusing quality)
// This is the end-to-end proof that the whole pipeline produces a working lens.

#include "celeris/design/metalens.hpp"

namespace celeris {

struct FocalAnalysis {
    // Transmission-weighted Strehl: peak(design) / peak(ideal lens, |t|=1
    // everywhere). Folds the aperture-averaged transmission loss into the number,
    // so strehl ~= (mean |t|)^2 * phase_strehl. Conservative; 1.0 = perfect.
    double strehl;
    // Phase Strehl: peak(design) / peak(same-amplitude, perfectly-phased lens).
    // Divides out the transmission loss, isolating wavefront/phase quality — the
    // classic optical Strehl ratio (relative to an aberration-free system with
    // the SAME amplitude apodization). 1.0 = diffraction-limited phase.
    double phase_strehl;
    double fwhm_um;            // focal-spot full width at half maximum
    double diffraction_limit_um;  // lambda * f / D for reference
    double encircled_energy;   // fraction of focal-window power within first null
};

// Analyze the focus of `lens` at the given focal length / wavelength / aperture
// diameter. The aperture is masked to a circle of the given diameter.
FocalAnalysis analyze_focus(const MetalensDesign& lens,
                            const UnitCellLibrary& lib, double focal_length_um,
                            double wavelength_um, double diameter_um);

// Full 2D focal-plane intensity map (the point-spread function), sampled on an
// n x n grid spanning +/- half_window_um about the axis at the focal plane.
struct PsfMap {
    int n;
    double half_window_um;
    std::vector<double> intensity;  // row-major, n*n, |E|^2
};
PsfMap compute_psf(const MetalensDesign& lens, const UnitCellLibrary& lib,
                   double focal_length_um, double wavelength_um,
                   double diameter_um, int n, double half_window_um);

// Low-level focal-plane propagation from an explicit pillar list (position +
// complex transmission), to an n x n window centered at (cx,cy) at distance z.
// Uses the GPU kernel when built with CUDA, else parallel CPU. This is the
// shared core used by polarization designs (independent t per polarization).
PsfMap propagate_pillars(const std::vector<double>& px,
                         const std::vector<double>& py,
                         const std::vector<cdouble>& t, double cx, double cy,
                         double z, double wavelength_um, int n,
                         double half_window_um);

// Off-axis PSF: illuminate the aperture with a plane wave tilted by angle_deg
// (in the x-z plane) and sample the focal plane in a window centered on the
// expected chief-ray landing point (x ~ f*tan(theta)). This is the per-field
// spot of a classic spot-vs-field diagram. The returned map's `cx_um` is that
// window-center x offset so callers can report the lateral spot shift.
struct FieldPsf {
    PsfMap psf;
    double angle_deg;
    double cx_um;        // window center (chief-ray landing) along x
    double rel_strehl;   // peak intensity relative to the on-axis peak
};
FieldPsf compute_psf_field(const MetalensDesign& lens, const UnitCellLibrary& lib,
                           double focal_length_um, double wavelength_um,
                           double diameter_um, double angle_deg, int n,
                           double half_window_um, double on_axis_peak = 0.0);

} // namespace celeris
