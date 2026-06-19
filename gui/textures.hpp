#pragma once
// OpenGL texture uploads for the GUI's rendered maps, plus the cyclic-phase
// colormap helper. All run on the main/GL thread.
#include <cstdint>

#include "celeris/analysis/focal.hpp"         // PsfMap
#include "celeris/analysis/throughfocus.hpp"  // ThroughFocus
#include "celeris/analysis/wavefront.hpp"     // WavefrontAnalysis
#include "celeris/design/metalens.hpp"        // MetalensDesign, UnitCellLibrary

namespace celeris::gui {

// Map hue in [0,1) (S=V=1) to RGB — the cyclic phase colormap.
void hue_to_rgb(double h, std::uint8_t& r, std::uint8_t& g, std::uint8_t& b);

void upload_psf_texture(const PsfMap& psf, unsigned int& tex);
void upload_wavefront_texture(const WavefrontAnalysis& wf, unsigned int& tex);
void upload_layout_texture(const MetalensDesign& d, const UnitCellLibrary& lib,
                           int mode, unsigned int& tex);
void upload_caustic_texture(const ThroughFocus& tf, unsigned int& tex);

} // namespace celeris::gui
