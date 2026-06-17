#pragma once
// GDSII export — the fabrication deliverable. A metalens design becomes a
// foundry-ready .gds file: each nanopillar is written as a polygon (BOUNDARY)
// on a mask layer. GDSII is the universal stream format every photonics/MEMS
// foundry accepts, so this is what turns a CELERIS design into something you
// can actually have made. No external dependency — we emit the binary records
// directly.

#include "celeris/design/metalens.hpp"

#include <string>

namespace celeris {

// Write a metalens pillar layout to a GDSII file. Each cell with fill >=
// min_fill becomes a square pillar (side = fill * period) on the given layer.
// Returns the number of pillars written, or -1 on file error.
int write_metalens_gds(const MetalensDesign& lens, const std::string& path,
                       int layer = 1, double min_fill = 0.05);

// Re-read a GDSII file and count BOUNDARY (polygon) records. Used to validate
// that an exported file is well-formed and has the expected element count.
// Returns -1 on error.
int gds_count_boundaries(const std::string& path);

} // namespace celeris
