#pragma once
// GDSII export — the fabrication deliverable. A metalens design becomes a
// foundry-ready .gds file: each nanopillar is written as a polygon (BOUNDARY)
// on a mask layer. GDSII is the universal stream format every photonics/MEMS
// foundry accepts, so this is what turns a CELERIS design into something you
// can actually have made. No external dependency — we emit the binary records
// directly.

#include "celeris/design/metalens.hpp"

#include <string>
#include <vector>

namespace celeris {

// A polygon read back from a GDSII file, vertices in microns.
struct GdsPolygon {
    std::vector<std::pair<double, double>> pts;  // (x,y) microns, closed ring
    int layer = 0;
};

// A parsed GDSII layout: all BOUNDARY polygons plus their bounding box (microns).
struct GdsLayout {
    std::vector<GdsPolygon> polygons;
    double min_x = 0, max_x = 0, min_y = 0, max_y = 0;
    bool ok = false;  // false if the file could not be parsed
};

// Parse a GDSII file into polygons (microns), reading the UNITS record so any
// foundry file scales correctly. Handles the subset CELERIS emits (BOUNDARY on
// layers) and ignores cell references / paths. Returns ok=false on read error.
GdsLayout read_gds(const std::string& path);

// Write a metalens pillar layout to a GDSII file. Each cell with fill >=
// min_fill becomes a square pillar (side = fill * period) on the given layer.
// Returns the number of pillars written, or -1 on file error.
int write_metalens_gds(const MetalensDesign& lens, const std::string& path,
                       int layer = 1, double min_fill = 0.05);

// Write a rectangular-pillar layout (independent x/y widths per cell, e.g. a
// polarization-multiplexed design). fill_x/fill_y are n_cells^2 row-major; each
// cell with max(fill) >= min_fill becomes a rectangle of size fill_x*period by
// fill_y*period. Returns the number of pillars written, or -1 on error.
int write_rect_gds(const std::string& path, int n_cells, double period_um,
                   const std::vector<double>& fill_x,
                   const std::vector<double>& fill_y, int layer = 1,
                   double min_fill = 0.05);

// Write a Pancharatnam-Berry layout: every site holds the SAME rectangular
// pillar (fill_x*period by fill_y*period) but ROTATED in-plane by rotation_rad
// (n_cells^2 row-major, CCW about the cell center). Each pillar is emitted as a
// rotated 4-vertex polygon -- the rotation IS the design, so a square footprint
// would be wrong. Returns the number of pillars written, or -1 on error.
int write_pb_gds(const std::string& path, int n_cells, double period_um,
                 double fill_x, double fill_y,
                 const std::vector<double>& rotation_rad, int layer = 1);

// Re-read a GDSII file and count BOUNDARY (polygon) records. Used to validate
// that an exported file is well-formed and has the expected element count.
// Returns -1 on error.
int gds_count_boundaries(const std::string& path);

} // namespace celeris
