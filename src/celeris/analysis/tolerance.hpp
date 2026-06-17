#pragma once
// Fabrication tolerance / yield analysis — how does the focus survive real-world
// manufacturing error? Lithography + etch leave each pillar a few nanometers off
// its design size. This Monte-Carlo perturbs every pillar's width by a Gaussian
// of a given sigma and measures the resulting Strehl, so you can answer "what
// CD (critical-dimension) control does this design need to hit spec?" — a hard
// requirement for defense/aerospace procurement.

#include "celeris/design/metalens.hpp"

#include <vector>

namespace celeris {

struct ToleranceResult {
    double sigma_nm;      // 1-sigma fabrication error on pillar size (nm)
    double mean_strehl;   // mean Strehl over the Monte-Carlo trials
    double std_strehl;    // standard deviation of Strehl
    double worst_strehl;  // worst-case trial
};

// For each sigma, run `trials` Monte-Carlo realizations (independent Gaussian CD
// error per pillar) and report the Strehl statistics. sigma_nm = 0 reproduces
// the nominal design. Trials run in parallel; `seed` makes it reproducible.
std::vector<ToleranceResult> analyze_tolerance(
    const MetalensDesign& lens, const UnitCellLibrary& lib, double focal_length_um,
    double wavelength_um, double diameter_um, const std::vector<double>& sigmas_nm,
    int trials, unsigned seed);

} // namespace celeris
