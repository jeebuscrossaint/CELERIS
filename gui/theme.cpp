#include "theme.hpp"

namespace celeris::gui {

ImVec4 rgb(int r, int g, int b, float a) {
    return ImVec4(r / 255.0f, g / 255.0f, b / 255.0f, a);
}

// CELERIS is a photonics design instrument, so the theme reads as one: a calm,
// cool-slate ground with a generous spacing rhythm and crisp (not bubbly) corners,
// and a single SPECTRAL accent -- laser cyan-teal -- spent only on interactive and
// active states (grabs, checks, selection, plots). The accent is the subject
// (spectral light); everything else stays quiet so it carries the identity alone.
void apply_theme(bool dark) {
    if (dark) ImGui::StyleColorsDark();
    else ImGui::StyleColorsLight();

    ImGuiStyle& s = ImGui::GetStyle();

    // Geometry: subtle rounding keeps it precise rather than toy-like; borders are
    // dropped on frames/children (contrast does the separating) so it stops reading
    // as a stack of Win32 dialog boxes. Spacing is a consistent, roomy rhythm.
    s.WindowRounding = 4; s.ChildRounding = 4; s.FrameRounding = 3;
    s.PopupRounding = 4;  s.GrabRounding = 3;  s.TabRounding = 4;
    s.ScrollbarRounding = 4;
    s.WindowBorderSize = 1; s.ChildBorderSize = 0; s.FrameBorderSize = 0;
    s.PopupBorderSize = 1;  s.SeparatorTextBorderSize = 1;
    s.WindowPadding = ImVec2(12, 10); s.FramePadding = ImVec2(10, 6);
    s.ItemSpacing = ImVec2(9, 7);     s.ItemInnerSpacing = ImVec2(8, 6);
    s.CellPadding = ImVec2(8, 5);     s.ScrollbarSize = 12; s.GrabMinSize = 11;
    s.WindowTitleAlign = ImVec2(0.0f, 0.5f);

    ImVec4* c = s.Colors;

    if (dark) {
        // Cool-slate instrument ground; two panel depths for quiet hierarchy.
        const ImVec4 ground = rgb(22, 25, 29),   panel = rgb(28, 32, 37),
                     field  = rgb(33, 38, 44),    field_hi = rgb(42, 48, 55),
                     text   = rgb(228, 231, 235), dim = rgb(138, 146, 155),
                     border = rgb(42, 47, 54),    row_alt = rgb(26, 30, 35),
                     button = rgb(38, 44, 51),    button_hi = rgb(47, 54, 62);
        // Spectral accent (cyan-teal) + a muted resting tint for selection fills.
        const ImVec4 accent = rgb(52, 211, 196), accent_dim = rgb(36, 179, 166),
                     sel_rest = rgb(29, 74, 71);
        c[ImGuiCol_WindowBg] = ground;       c[ImGuiCol_ChildBg] = ImVec4(0, 0, 0, 0);
        c[ImGuiCol_PopupBg] = panel;         c[ImGuiCol_MenuBarBg] = panel;
        c[ImGuiCol_Text] = text;             c[ImGuiCol_TextDisabled] = dim;
        c[ImGuiCol_Border] = border;         c[ImGuiCol_BorderShadow] = ImVec4(0, 0, 0, 0);
        c[ImGuiCol_FrameBg] = field;         c[ImGuiCol_FrameBgHovered] = field_hi;
        c[ImGuiCol_FrameBgActive] = field_hi;
        c[ImGuiCol_Button] = button;         c[ImGuiCol_ButtonHovered] = button_hi;
        c[ImGuiCol_ButtonActive] = sel_rest;
        c[ImGuiCol_Header] = sel_rest;       c[ImGuiCol_HeaderHovered] = button_hi;
        c[ImGuiCol_HeaderActive] = accent_dim;
        c[ImGuiCol_SliderGrab] = accent;     c[ImGuiCol_SliderGrabActive] = accent_dim;
        c[ImGuiCol_CheckMark] = accent;
        c[ImGuiCol_Tab] = panel;             c[ImGuiCol_TabHovered] = button_hi;
        c[ImGuiCol_TabSelected] = field_hi;    c[ImGuiCol_TabDimmed] = ground;
        c[ImGuiCol_TabDimmedSelected] = panel;
        c[ImGuiCol_TitleBg] = panel;         c[ImGuiCol_TitleBgActive] = panel;
        c[ImGuiCol_TitleBgCollapsed] = ground;
        c[ImGuiCol_TableHeaderBg] = panel;
        c[ImGuiCol_TableBorderStrong] = border;
        c[ImGuiCol_TableBorderLight] = rgb(35, 40, 46);
        c[ImGuiCol_TableRowBg] = ImVec4(0, 0, 0, 0);
        c[ImGuiCol_TableRowBgAlt] = row_alt;
        c[ImGuiCol_Separator] = border;      c[ImGuiCol_SeparatorHovered] = accent_dim;
        c[ImGuiCol_SeparatorActive] = accent;
        c[ImGuiCol_ResizeGrip] = border;     c[ImGuiCol_ResizeGripHovered] = accent_dim;
        c[ImGuiCol_ResizeGripActive] = accent;
        c[ImGuiCol_ScrollbarBg] = ImVec4(0, 0, 0, 0);
        c[ImGuiCol_ScrollbarGrab] = rgb(46, 52, 59);
        c[ImGuiCol_ScrollbarGrabHovered] = rgb(58, 66, 75);
        c[ImGuiCol_ScrollbarGrabActive] = accent_dim;
        c[ImGuiCol_PlotLines] = accent;      c[ImGuiCol_PlotLinesHovered] = accent_dim;
        c[ImGuiCol_PlotHistogram] = accent;  c[ImGuiCol_PlotHistogramHovered] = accent_dim;
        c[ImGuiCol_DragDropTarget] = accent;
        c[ImGuiCol_NavCursor] = accent;
    } else {
        // Clean "lab" light: cool off-white ground, white fields, deeper accent so
        // the spectral cyan-teal keeps contrast on a bright background.
        const ImVec4 ground = rgb(245, 246, 247), panel = rgb(255, 255, 255),
                     field  = rgb(255, 255, 255),  field_hi = rgb(236, 240, 242),
                     text   = rgb(26, 30, 35),      dim = rgb(138, 145, 153),
                     border = rgb(214, 218, 222),   row_alt = rgb(240, 242, 244),
                     button = rgb(236, 238, 240),   button_hi = rgb(226, 233, 232);
        const ImVec4 accent = rgb(15, 181, 170), accent_dim = rgb(12, 143, 134),
                     sel_rest = rgb(205, 237, 233);
        c[ImGuiCol_WindowBg] = ground;       c[ImGuiCol_ChildBg] = ImVec4(0, 0, 0, 0);
        c[ImGuiCol_PopupBg] = panel;         c[ImGuiCol_MenuBarBg] = panel;
        c[ImGuiCol_Text] = text;             c[ImGuiCol_TextDisabled] = dim;
        c[ImGuiCol_Border] = border;         c[ImGuiCol_BorderShadow] = ImVec4(0, 0, 0, 0);
        c[ImGuiCol_FrameBg] = field;         c[ImGuiCol_FrameBgHovered] = field_hi;
        c[ImGuiCol_FrameBgActive] = field_hi;
        c[ImGuiCol_Button] = button;         c[ImGuiCol_ButtonHovered] = button_hi;
        c[ImGuiCol_ButtonActive] = sel_rest;
        c[ImGuiCol_Header] = sel_rest;       c[ImGuiCol_HeaderHovered] = button_hi;
        c[ImGuiCol_HeaderActive] = accent;
        c[ImGuiCol_SliderGrab] = accent;     c[ImGuiCol_SliderGrabActive] = accent_dim;
        c[ImGuiCol_CheckMark] = accent_dim;
        c[ImGuiCol_Tab] = ground;            c[ImGuiCol_TabHovered] = button_hi;
        c[ImGuiCol_TabSelected] = panel;       c[ImGuiCol_TabDimmed] = ground;
        c[ImGuiCol_TabDimmedSelected] = panel;
        c[ImGuiCol_TitleBg] = panel;         c[ImGuiCol_TitleBgActive] = panel;
        c[ImGuiCol_TitleBgCollapsed] = ground;
        c[ImGuiCol_TableHeaderBg] = rgb(238, 240, 242);
        c[ImGuiCol_TableBorderStrong] = border;
        c[ImGuiCol_TableBorderLight] = rgb(228, 231, 234);
        c[ImGuiCol_TableRowBg] = ImVec4(0, 0, 0, 0);
        c[ImGuiCol_TableRowBgAlt] = row_alt;
        c[ImGuiCol_Separator] = border;      c[ImGuiCol_SeparatorHovered] = accent_dim;
        c[ImGuiCol_SeparatorActive] = accent_dim;
        c[ImGuiCol_ResizeGrip] = border;     c[ImGuiCol_ResizeGripHovered] = accent_dim;
        c[ImGuiCol_ResizeGripActive] = accent_dim;
        c[ImGuiCol_ScrollbarBg] = ImVec4(0, 0, 0, 0);
        c[ImGuiCol_ScrollbarGrab] = rgb(206, 211, 216);
        c[ImGuiCol_ScrollbarGrabHovered] = rgb(184, 190, 196);
        c[ImGuiCol_ScrollbarGrabActive] = accent_dim;
        c[ImGuiCol_PlotLines] = accent_dim;  c[ImGuiCol_PlotLinesHovered] = accent;
        c[ImGuiCol_PlotHistogram] = accent_dim; c[ImGuiCol_PlotHistogramHovered] = accent;
        c[ImGuiCol_DragDropTarget] = accent_dim;
        c[ImGuiCol_NavCursor] = accent_dim;
    }
}

} // namespace celeris::gui
