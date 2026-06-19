#include "theme.hpp"

namespace celeris::gui {

ImVec4 rgb(int r, int g, int b, float a) {
    return ImVec4(r / 255.0f, g / 255.0f, b / 255.0f, a);
}

void apply_theme(bool dark) {
    if (dark) ImGui::StyleColorsDark();
    else ImGui::StyleColorsLight();

    ImGuiStyle& s = ImGui::GetStyle();
    s.WindowRounding = 0; s.ChildRounding = 0; s.FrameRounding = 0;
    s.PopupRounding = 0;  s.GrabRounding = 0;  s.TabRounding = 0;
    s.ScrollbarRounding = 0;
    s.WindowBorderSize = 1; s.ChildBorderSize = 1; s.FrameBorderSize = 1;
    s.PopupBorderSize = 1;
    s.WindowPadding = ImVec2(6, 6);  s.FramePadding = ImVec2(6, 3);
    s.ItemSpacing = ImVec2(6, 4);    s.ItemInnerSpacing = ImVec2(5, 4);
    s.CellPadding = ImVec2(6, 3);    s.ScrollbarSize = 14; s.GrabMinSize = 10;

    ImVec4* c = s.Colors;
    const ImVec4 sel = rgb(0, 120, 215);  // classic Windows-blue selection (both)
    if (!dark) {
        const ImVec4 face = rgb(238, 238, 238), field = rgb(255, 255, 255),
                     text = rgb(20, 20, 20), border = rgb(158, 158, 158),
                     head = rgb(222, 222, 222),
                     hov = rgb(229, 241, 251), prs = rgb(204, 228, 247);
        c[ImGuiCol_WindowBg] = face;   c[ImGuiCol_ChildBg] = face;
        c[ImGuiCol_MenuBarBg] = head;  c[ImGuiCol_PopupBg] = field;
        c[ImGuiCol_Text] = text;       c[ImGuiCol_Border] = border;
        c[ImGuiCol_FrameBg] = field;   c[ImGuiCol_FrameBgHovered] = hov;
        c[ImGuiCol_FrameBgActive] = prs;
        c[ImGuiCol_Button] = rgb(225, 225, 225); c[ImGuiCol_ButtonHovered] = hov;
        c[ImGuiCol_ButtonActive] = prs;
        c[ImGuiCol_Header] = prs;      c[ImGuiCol_HeaderHovered] = hov;
        c[ImGuiCol_HeaderActive] = sel;
        c[ImGuiCol_SliderGrab] = sel;  c[ImGuiCol_SliderGrabActive] = rgb(0, 102, 184);
        c[ImGuiCol_CheckMark] = sel;
        c[ImGuiCol_TableHeaderBg] = head;
        c[ImGuiCol_TableBorderStrong] = border;
        c[ImGuiCol_TableBorderLight] = rgb(200, 200, 200);
        c[ImGuiCol_TableRowBg] = field; c[ImGuiCol_TableRowBgAlt] = rgb(247, 247, 247);
        c[ImGuiCol_TitleBg] = head;     c[ImGuiCol_TitleBgActive] = head;
        c[ImGuiCol_Separator] = border; c[ImGuiCol_PlotLines] = sel;
        c[ImGuiCol_ScrollbarBg] = face;
    } else {
        const ImVec4 face = rgb(37, 37, 38), field = rgb(24, 24, 25),
                     text = rgb(225, 225, 225), border = rgb(70, 70, 72),
                     head = rgb(45, 45, 47),
                     hov = rgb(58, 70, 90), prs = rgb(38, 79, 120);
        c[ImGuiCol_WindowBg] = face;   c[ImGuiCol_ChildBg] = face;
        c[ImGuiCol_MenuBarBg] = head;  c[ImGuiCol_PopupBg] = field;
        c[ImGuiCol_Text] = text;       c[ImGuiCol_Border] = border;
        c[ImGuiCol_FrameBg] = field;   c[ImGuiCol_FrameBgHovered] = hov;
        c[ImGuiCol_FrameBgActive] = prs;
        c[ImGuiCol_Button] = rgb(55, 55, 57); c[ImGuiCol_ButtonHovered] = hov;
        c[ImGuiCol_ButtonActive] = prs;
        c[ImGuiCol_Header] = prs;      c[ImGuiCol_HeaderHovered] = hov;
        c[ImGuiCol_HeaderActive] = sel;
        c[ImGuiCol_SliderGrab] = sel;  c[ImGuiCol_SliderGrabActive] = rgb(0, 102, 184);
        c[ImGuiCol_CheckMark] = sel;
        c[ImGuiCol_TableHeaderBg] = head;
        c[ImGuiCol_TableBorderStrong] = border;
        c[ImGuiCol_TableBorderLight] = rgb(60, 60, 62);
        c[ImGuiCol_TableRowBg] = field; c[ImGuiCol_TableRowBgAlt] = rgb(32, 32, 33);
        c[ImGuiCol_TitleBg] = head;     c[ImGuiCol_TitleBgActive] = head;
        c[ImGuiCol_Separator] = border; c[ImGuiCol_PlotLines] = rgb(90, 170, 255);
        c[ImGuiCol_ScrollbarBg] = face;
    }
}

} // namespace celeris::gui
