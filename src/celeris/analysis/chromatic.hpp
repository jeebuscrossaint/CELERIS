#pragma once
/// @file
/// @brief Chromatic analysis: how the focus of a finished lens moves with wavelength.
///
/// A standard phase-profile metalens is strongly chromatic: its focal length scales as
/// f(λ) ≈ f₀·λ₀/λ. This is the #1 limitation of metalenses (and why achromatic
/// designs are a hot research area), so quantifying it is essential.
///
/// Two models are provided. analyze_chromatic() holds the pillars' imparted phase
/// fixed at the design value (geometric-phase approximation) and varies only the
/// propagation wavenumber k(λ); this captures the dominant chromatic focal shift.
/// analyze_chromatic_dispersive() rebuilds the RCWA library at each wavelength, so
/// the meta-atoms' own material and waveguide dispersion is included.

#include "celeris/design/metalens.hpp"

#include <functional>
#include <vector>

namespace celeris {

/// @brief The located focus of a lens at one wavelength.
struct ChromaticPoint {
    double wavelength_um;     ///< Vacuum wavelength of this sample (µm).
    double focal_length_um;   ///< Axial position where on-axis intensity actually peaks (µm).
    double rel_peak;          ///< Peak on-axis intensity relative to the design-wavelength peak at the design focus.
};

/// @brief Fast chromatic sweep: fixed pillar transmission, wavelength-dependent propagation.
///
/// Sweeps wavelength from lambda_min to lambda_max (n points) and locates the
/// focus at each by scanning on-axis intensity (Rayleigh-Sommerfeld sum) over
/// 0.55–1.6 × the analytic estimate f₀·λ₀/λ.
///
/// FAST model: the pillars' imparted transmission is held fixed at the design
/// value and only the propagation wavenumber k(λ) varies. Captures the dominant
/// chromatic focal shift but ignores the meta-atoms' own dispersion.
/// @param lens               The finished design (fill map on the lens grid).
/// @param lib                The library the design was built from (fill → transmission).
/// @param design_focal       Nominal focal length (µm).
/// @param design_wavelength  Nominal design wavelength (µm).
/// @param diameter           Circular aperture diameter (µm).
/// @param lambda_min,lambda_max  Sweep range (µm).
/// @param n                  Number of wavelengths (n >= 2).
/// @return One ChromaticPoint per wavelength, in sweep order.
std::vector<ChromaticPoint> analyze_chromatic(const MetalensDesign& lens,
                                              const UnitCellLibrary& lib,
                                              double design_focal,
                                              double design_wavelength,
                                              double diameter, double lambda_min,
                                              double lambda_max, int n);

/// @brief Rigorous chromatic sweep: the RCWA library is rebuilt at every wavelength.
///
/// RIGOROUS model: rebuild the RCWA unit-cell library at every wavelength via
/// build_lib_at(λ), so each pillar's transmission reflects the true material AND
/// waveguide dispersion at that wavelength (not just the propagation phase). This
/// is heavier (one RCWA library sweep per wavelength) but is the physically
/// correct chromatic response of a real material metalens.
/// @param lens          The finished design.
/// @param build_lib_at  Returns the library at a given wavelength (µm). It must use
///                      the same fill grid as the design (the same builder that
///                      created the lens), so transmission_for_fill maps each pillar
///                      correctly.
/// @param design_focal,design_wavelength  Nominal design point (µm).
/// @param diameter      Circular aperture diameter (µm).
/// @param lambda_min,lambda_max  Sweep range (µm).
/// @param n             Number of wavelengths (n >= 2).
/// @return One ChromaticPoint per wavelength, in sweep order.
std::vector<ChromaticPoint> analyze_chromatic_dispersive(
    const MetalensDesign& lens,
    const std::function<UnitCellLibrary(double)>& build_lib_at,
    double design_focal, double design_wavelength, double diameter,
    double lambda_min, double lambda_max, int n);

} // namespace celeris
