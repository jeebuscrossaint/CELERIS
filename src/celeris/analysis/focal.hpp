#pragma once
/// @file
/// @brief Focal-plane analysis: does the designed lens actually focus, and how well?
///
/// Takes a finished metalens design, builds its complex aperture field (each
/// pillar's transmission from the library), propagates it to the focal plane by
/// Rayleigh-Sommerfeld summation, and reports the standard optical metrics:
///   - Strehl ratio: peak intensity vs an ideal (perfect-phase) lens (1.0 = perfect)
///   - focal-spot FWHM, compared to the diffraction limit ~ lambda*f/D
///   - encircled energy within the first Airy null (focusing quality)
///
/// This is the end-to-end proof that the whole pipeline produces a working lens.

#include "celeris/design/metalens.hpp"

namespace celeris {

/// @brief Scalar focus metrics of one lens at one wavelength.
struct FocalAnalysis {
    /// Transmission-weighted Strehl: peak(design) / peak(ideal lens, |t|=1
    /// everywhere). Folds the aperture-averaged transmission loss into the number,
    /// so strehl ~= (mean |t|)^2 * phase_strehl. Conservative; 1.0 = perfect.
    double strehl;
    /// Phase Strehl: peak(design) / peak(same-amplitude, perfectly-phased lens).
    /// Divides out the transmission loss, isolating wavefront/phase quality — the
    /// classic optical Strehl ratio (relative to an aberration-free system with
    /// the SAME amplitude apodization). 1.0 = diffraction-limited phase.
    double phase_strehl;
    double fwhm_um;               ///< Focal-spot full width at half maximum, central row (µm).
    double diffraction_limit_um;  ///< lambda * f / D, for reference (µm).
    double encircled_energy;      ///< Fraction of focal-window power within the first Airy null (1.22 lambda f / D).
};

/// @brief Analyze the focus of `lens`.
///
/// The aperture is masked to a circle of the given diameter. The focal plane is
/// sampled on a 121 x 121 grid of half-width max(6 lambda f / D, 4 µm); the GPU
/// kernel is used when built with CUDA.
/// @param lens            The finished design.
/// @param lib             The library the design was built from.
/// @param focal_length_um Focal length (µm); the analysis plane.
/// @param wavelength_um   Vacuum wavelength (µm).
/// @param diameter_um     Circular aperture diameter (µm).
/// @return The focus metrics.
FocalAnalysis analyze_focus(const MetalensDesign& lens,
                            const UnitCellLibrary& lib, double focal_length_um,
                            double wavelength_um, double diameter_um);

/// @brief A sampled focal-plane intensity map (the point-spread function).
///
/// Samples span ±half_window_um about the window centre on an n x n grid.
struct PsfMap {
    int n;                          ///< Samples per side.
    double half_window_um;          ///< Half-width of the sampled window (µm).
    std::vector<double> intensity;  ///< Row-major n*n map of |E|^2 (unnormalized).
};

/// @brief Full 2D on-axis PSF of `lens` at its focal plane.
/// @param lens            The finished design.
/// @param lib             The library the design was built from.
/// @param focal_length_um Distance to the sampled plane (µm).
/// @param wavelength_um   Vacuum wavelength (µm).
/// @param diameter_um     Circular aperture diameter (µm).
/// @param n               Samples per side (n >= 2).
/// @param half_window_um  Half-width of the window about the axis (µm).
/// @return The n x n intensity map.
PsfMap compute_psf(const MetalensDesign& lens, const UnitCellLibrary& lib,
                   double focal_length_um, double wavelength_um,
                   double diameter_um, int n, double half_window_um);

/// @brief Low-level focal-plane propagation from an explicit pillar list.
///
/// Each pillar is a position plus a complex transmission; the field is summed to
/// an n x n window centered at (cx, cy) at distance z. Uses the GPU kernel when
/// built with CUDA, else propagate_pillars_cpu(). This is the shared core used by
/// polarization designs (independent t per polarization).
/// @param px,py  Pillar centre coordinates (µm), one entry per pillar.
/// @param t      Complex transmission of each pillar.
/// @param cx,cy  Centre of the sampled window in the observation plane (µm).
/// @param z      Propagation distance to the observation plane (µm).
/// @param wavelength_um  Free-space wavelength (µm).
/// @param n      Samples per side of the square window (n >= 2).
/// @param half_window_um Half-width of the window (µm).
/// @return n x n map of |E|^2.
PsfMap propagate_pillars(const std::vector<double>& px,
                         const std::vector<double>& py,
                         const std::vector<cdouble>& t, double cx, double cy,
                         double z, double wavelength_um, int n,
                         double half_window_um);

/// CPU-only focal-plane propagation: the double-precision Rayleigh-Sommerfeld
/// sum behind propagate_pillars, parallel over focal-plane rows, never
/// dispatched to the GPU even in a CUDA build. propagate_pillars and
/// compute_psf fall back to it; `celeris psfbench` times it as the CPU
/// baseline against the float GPU kernel.
/// @param px,py  Pillar centre coordinates (µm), one entry per pillar.
/// @param t      Complex transmission of each pillar.
/// @param cx,cy  Centre of the sampled window in the observation plane (µm).
/// @param z      Propagation distance to the observation plane (µm).
/// @param wavelength_um  Free-space wavelength (µm).
/// @param n      Samples per side of the square window (n >= 2).
/// @param half_window_um Half-width of the window (µm).
/// @return n x n map of |E|^2; all zeros if there are no pillars.
PsfMap propagate_pillars_cpu(const std::vector<double>& px,
                             const std::vector<double>& py,
                             const std::vector<cdouble>& t, double cx, double cy,
                             double z, double wavelength_um, int n,
                             double half_window_um);

/// @brief Off-axis PSF at one field angle.
///
/// The aperture is lit by a plane wave tilted by angle_deg (in the x-z plane),
/// and the focal plane is sampled in a window centered on the expected chief-ray
/// landing point (x ~ f*tan(theta)). This is the per-field spot of a classic
/// spot-vs-field diagram.
struct FieldPsf {
    PsfMap psf;          ///< The sampled off-axis PSF.
    double angle_deg;    ///< Field angle (degrees).
    double cx_um;        ///< Window center (chief-ray landing) along x (µm); gives the lateral spot shift.
    double rel_strehl;   ///< Peak intensity relative to the on-axis peak.
};

/// @brief Compute the off-axis PSF of `lens` at one field angle.
/// @param lens            The finished design.
/// @param lib             The library the design was built from.
/// @param focal_length_um Focal length (µm); the sampled plane.
/// @param wavelength_um   Vacuum wavelength (µm).
/// @param diameter_um     Circular aperture diameter (µm).
/// @param angle_deg       Field angle in the x-z plane (degrees).
/// @param n               Samples per side (n >= 2).
/// @param half_window_um  Half-width of the window about the chief ray (µm).
/// @param on_axis_peak    On-axis peak intensity for normalizing rel_strehl;
///                        0 (the default) sets rel_strehl = 1.
/// @return The off-axis PSF and its metadata.
FieldPsf compute_psf_field(const MetalensDesign& lens, const UnitCellLibrary& lib,
                           double focal_length_um, double wavelength_um,
                           double diameter_um, double angle_deg, int n,
                           double half_window_um, double on_axis_peak = 0.0);

} // namespace celeris
