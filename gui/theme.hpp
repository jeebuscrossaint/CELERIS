#pragma once
// CELERIS GUI theme. Utilitarian engineering-tool look (think OpticStudio/MFC):
// square corners, thin borders, dense spacing. Light and dark variants.
#include "imgui.h"

namespace celeris::gui {

// Build an ImVec4 from 0-255 RGB(A).
ImVec4 rgb(int r, int g, int b, float a = 1.0f);

// Apply the theme to the current ImGui context. dark = dark-mode palette.
void apply_theme(bool dark);

} // namespace celeris::gui
