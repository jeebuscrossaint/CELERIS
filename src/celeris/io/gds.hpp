#pragma once
/// @file
/// @brief GDSII export of metalens layouts, and a reader for checking the output.
///
/// GDSII export — the fabrication deliverable. A metalens design becomes a
/// foundry-ready .gds file: each nanopillar is written as a polygon (BOUNDARY)
/// on a mask layer. GDSII is the universal stream format every photonics/MEMS
/// foundry accepts, so this is what turns a CELERIS design into something you
/// can actually have made. No external dependency — we emit the binary records
/// directly.

#include "celeris/design/metalens.hpp"

#include <string>
#include <vector>

namespace celeris {

/// @brief A polygon read back from a GDSII file, vertices in µm.
struct GdsPolygon {
    std::vector<std::pair<double, double>> pts;  ///< (x, y) vertices in µm, closed ring.
    int layer = 0;                               ///< GDSII layer number.
};

/// @brief A parsed GDSII layout: all BOUNDARY polygons and their bounding box (µm).
struct GdsLayout {
    std::vector<GdsPolygon> polygons;  ///< Every BOUNDARY element in the file.
    double min_x = 0;                  ///< Bounding box, smallest x (µm).
    double max_x = 0;                  ///< Bounding box, largest x (µm).
    double min_y = 0;                  ///< Bounding box, smallest y (µm).
    double max_y = 0;                  ///< Bounding box, largest y (µm).
    bool ok = false;                   ///< False if the file could not be parsed.
};

/// @brief Parse a GDSII file into polygons in µm.
///
/// Reads the UNITS record, so files with any database unit scale correctly.
/// Handles the subset CELERIS writes (BOUNDARY elements on layers) and ignores
/// cell references and paths.
/// @param path GDSII file to read.
/// @return The layout; `ok` is false on a read or parse error.
GdsLayout read_gds(const std::string& path);

/// @brief Write a propagation-phase metalens as square pillars to a GDSII file.
///
/// Each cell with fill ≥ `min_fill` becomes a square of side fill · period.
/// @param lens     The designed lens (cell grid, period, per-cell fill).
/// @param path     Output file.
/// @param layer    GDSII layer for the pillars.
/// @param min_fill Cells below this fill are left empty.
/// @return Number of pillars written, or −1 on a file error.
int write_metalens_gds(const MetalensDesign& lens, const std::string& path,
                       int layer = 1, double min_fill = 0.05);

/// @brief Write a rectangular-pillar layout with independent x and y widths per cell.
///
/// Used for polarization-multiplexed designs. Each cell with
/// max(fill_x, fill_y) ≥ `min_fill` becomes a fill_x · period by fill_y · period
/// rectangle.
/// @param path      Output file.
/// @param n_cells   Cells per side of the square grid.
/// @param period_um Cell period (µm).
/// @param fill_x    Per-cell x fill, n_cells² values, row-major.
/// @param fill_y    Per-cell y fill, n_cells² values, row-major.
/// @param layer     GDSII layer for the pillars.
/// @param min_fill  Cells below this fill are left empty.
/// @return Number of pillars written, or −1 on error.
int write_rect_gds(const std::string& path, int n_cells, double period_um,
                   const std::vector<double>& fill_x,
                   const std::vector<double>& fill_y, int layer = 1,
                   double min_fill = 0.05);

/// @brief Write a Pancharatnam–Berry layout: one rectangular pillar, rotated per site.
///
/// Every site holds the same fill_x · period by fill_y · period rectangle,
/// rotated in-plane about the cell centre. Each pillar is written as a rotated
/// 4-vertex polygon, because the rotation carries the phase.
/// @param path         Output file.
/// @param n_cells      Cells per side of the square grid.
/// @param period_um    Cell period (µm).
/// @param fill_x       Pillar x fill (fraction of the period).
/// @param fill_y       Pillar y fill (fraction of the period).
/// @param rotation_rad Per-site rotation, counter-clockwise, n_cells² values, row-major.
/// @param layer        GDSII layer for the pillars.
/// @return Number of pillars written, or −1 on error.
int write_pb_gds(const std::string& path, int n_cells, double period_um,
                 double fill_x, double fill_y,
                 const std::vector<double>& rotation_rad, int layer = 1);

/// @brief Write an achromatic Pancharatnam–Berry layout: per-site footprint and rotation.
///
/// Generalizes write_pb_gds() to the per-site dispersive atom an achromatic PB
/// design picks. All pillars share one etch depth, so the result is a
/// single-layer, single-etch mask. Each pillar is a rotated 4-vertex polygon.
/// @param path         Output file.
/// @param n_cells      Cells per side of the square grid.
/// @param period_um    Cell period (µm).
/// @param fill_x       Per-site x fill, n_cells² values, row-major.
/// @param fill_y       Per-site y fill, n_cells² values, row-major.
/// @param rotation_rad Per-site rotation, counter-clockwise, n_cells² values, row-major.
/// @param layer        GDSII layer for the pillars.
/// @param min_fill     Sites whose x and y fill are both below this are left empty.
/// @return Number of pillars written, or −1 on error.
int write_pb_rect_gds(const std::string& path, int n_cells, double period_um,
                      const std::vector<double>& fill_x,
                      const std::vector<double>& fill_y,
                      const std::vector<double>& rotation_rad, int layer = 1,
                      double min_fill = 0.05);

/// @brief Count the BOUNDARY (polygon) records in a GDSII file.
///
/// Used to check that an exported file is well formed and has the expected
/// number of elements.
/// @param path GDSII file to read.
/// @return Number of BOUNDARY records, or −1 on error.
int gds_count_boundaries(const std::string& path);

} // namespace celeris
