#pragma once
// Modulation Transfer Function (MTF) — the headline imaging-quality metric in
// every lens datasheet: how much image contrast survives at each spatial
// frequency. The optical transfer function is the autocorrelation of the exit
// pupil P = |t|·exp(i·2π·OPD); MTF = |OTF|. We also return the diffraction-
// limited MTF (perfect pupil) as the reference ceiling.

#include "celeris/design/metalens.hpp"

#include <vector>

namespace celeris {

struct MtfCurve {
    std::vector<float> freq_cyc_mm;  // spatial frequency, cycles/mm
    std::vector<float> mtf;          // lens MTF (0..1)
    std::vector<float> mtf_ideal;    // diffraction-limited MTF at same freqs
    double cutoff_cyc_mm;            // diffraction cutoff = D/(lambda*f)
};

MtfCurve analyze_mtf(const MetalensDesign& lens, const UnitCellLibrary& lib,
                     double focal_length_um, double wavelength_um,
                     double diameter_um);

} // namespace celeris
