#pragma once
/// @file
/// @brief Modulation Transfer Function (MTF): image contrast vs spatial frequency.
///
/// The MTF is the headline imaging-quality metric in every lens datasheet: how
/// much image contrast survives at each spatial frequency. The optical transfer
/// function is the autocorrelation of the exit pupil P = |t|·exp(i·2π·OPD);
/// MTF = |OTF|. We also return the diffraction-limited MTF (perfect pupil) as the
/// reference ceiling.

#include "celeris/design/metalens.hpp"

#include <vector>

namespace celeris {

/// @brief A sampled MTF curve with its diffraction-limited reference.
struct MtfCurve {
    std::vector<float> freq_cyc_mm;  ///< Spatial frequency (cycles/mm).
    std::vector<float> mtf;          ///< Lens MTF (0..1).
    std::vector<float> mtf_ideal;    ///< Diffraction-limited MTF at the same frequencies.
    double cutoff_cyc_mm;            ///< Diffraction cutoff D/(lambda*f) (cycles/mm).
};

/// @brief Compute the MTF of `lens` along x.
///
/// The pupil autocorrelation is taken along x at shifts of whole lens cells, so
/// frequencies are spaced by period/(lambda f) and run from 0 to the cutoff.
/// @param lens            The finished design.
/// @param lib             The library the design was built from.
/// @param focal_length_um Focal length (µm).
/// @param wavelength_um   Vacuum wavelength (µm).
/// @param diameter_um     Circular aperture diameter (µm).
/// @return The MTF curve; both curves are normalized to 1 at zero frequency.
MtfCurve analyze_mtf(const MetalensDesign& lens, const UnitCellLibrary& lib,
                     double focal_length_um, double wavelength_um,
                     double diameter_um);

} // namespace celeris
