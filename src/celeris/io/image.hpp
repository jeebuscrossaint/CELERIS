#pragma once
// Minimal grayscale image output (PGM, P5 binary) — no dependencies. Used to
// dump the focal-plane PSF / field maps so results can be viewed in any image
// tool without a GUI.

#include <string>
#include <vector>

namespace celeris {

// Write an n_x by n_y grayscale PGM. `values` is row-major; it is normalized to
// 0..255 by its own maximum. If `gamma` != 1, applies value^(1/gamma) before
// scaling (gamma > 1 brightens faint detail like PSF sidelobes). Returns false
// on file error.
bool write_pgm(const std::string& path, int nx, int ny,
               const std::vector<double>& values, double gamma = 1.0);

} // namespace celeris
