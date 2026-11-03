#pragma once
/// @file
/// @brief Through-focus analysis: on-axis intensity vs z, depth of focus, and the x-z caustic.
///
/// The "focus" axis of optical analysis. Scans the on-axis intensity along z
/// around the nominal focus to show how sharply the lens focuses and over what
/// axial range (depth of focus). Standard tab in real optical-design tools.

#include "celeris/design/metalens.hpp"

#include <vector>

namespace celeris {

/// @brief On-axis through-focus scan plus the longitudinal caustic.
struct ThroughFocus {
    std::vector<float> z_um;        ///< Axial positions of the scan (µm).
    std::vector<float> intensity;   ///< On-axis intensity at each z, normalized to its peak.
    double z_peak_um;               ///< Where on-axis intensity actually peaks (µm).
    double dof_um;                  ///< Depth of focus: axial FWHM of the peak (µm).

    /// @name Longitudinal caustic
    /// |E|^2 on an (nz x nx) grid in the y=0 plane (row = z, col = x), normalized
    /// to its peak. Empty if it could not be computed. Shows the focusing cone --
    /// the standard "x-z" focus view.
    ///@{
    int caustic_nx = 0;             ///< Caustic samples along x.
    int caustic_nz = 0;             ///< Caustic samples along z.
    double caustic_xmin = 0;        ///< Lateral extent, low edge (µm).
    double caustic_xmax = 0;        ///< Lateral extent, high edge (µm).
    double caustic_zmin = 0;        ///< Axial extent, low edge (µm).
    double caustic_zmax = 0;        ///< Axial extent, high edge (µm).
    std::vector<float> caustic;     ///< Row-major nz x nx intensity, normalized to its peak.
    ///@}
};

/// @brief Scan the focus of `lens` along z.
///
/// The scan covers focal_length_um ± max(2 lambda / NA², f/4) (lower end clipped
/// at f/10) with 161 points, so low-NA lenses with long depth of focus are
/// covered. The caustic is computed on the GPU when available, else on a coarser
/// CPU grid.
/// @param lens            The finished design.
/// @param lib             The library the design was built from.
/// @param focal_length_um Nominal focal length (µm).
/// @param wavelength_um   Vacuum wavelength (µm).
/// @param diameter_um     Circular aperture diameter (µm).
/// @return The scan, its peak, depth of focus, and caustic.
ThroughFocus analyze_through_focus(const MetalensDesign& lens,
                                   const UnitCellLibrary& lib,
                                   double focal_length_um, double wavelength_um,
                                   double diameter_um);

} // namespace celeris
