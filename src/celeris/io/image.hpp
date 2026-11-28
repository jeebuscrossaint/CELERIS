#pragma once
/// @file
/// @brief Grayscale PGM image output for PSF and field maps.
///
/// Minimal grayscale image output (PGM, P5 binary) — no dependencies. Used to
/// dump the focal-plane PSF / field maps so results can be viewed in any image
/// tool without a GUI.

#include <string>
#include <vector>

namespace celeris {

/// @brief Write an 8-bit grayscale PGM (P5) image.
///
/// The values are normalized to 0–255 by their own maximum. A `gamma` other
/// than 1 applies value^(1/gamma) before scaling; gamma > 1 brightens faint
/// detail such as PSF sidelobes.
/// @param path   Output file.
/// @param nx     Image width (pixels).
/// @param ny     Image height (pixels).
/// @param values nx · ny values, row-major.
/// @param gamma  Display gamma.
/// @return False on a file error.
bool write_pgm(const std::string& path, int nx, int ny,
               const std::vector<double>& values, double gamma = 1.0);

} // namespace celeris
