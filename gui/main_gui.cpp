// CELERIS desktop GUI (Dear ImGui + GLFW + OpenGL).
//
// A thin front-end over celeris_core: set lens parameters, hit Design, and the
// full pipeline (RCWA library -> lens -> focal/chromatic analysis -> PSF) runs
// on a worker thread while the UI stays live. Results show as metrics, a
// rendered point-spread-function image, and a chromatic focal-shift plot.

#include <GLFW/glfw3.h>
#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"

#include <algorithm>
#include <atomic>
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <format>
#include <fstream>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#include "celeris/analysis/chromatic.hpp"
#include "celeris/analysis/field.hpp"
#include "celeris/analysis/focal.hpp"
#include "celeris/analysis/tolerance.hpp"
#include "celeris/design/metalens.hpp"
#include "celeris/io/gds.hpp"
#include "celeris/materials/database.hpp"

using namespace celeris;

namespace {

struct Params {
    float focal = 50, diameter = 20, wavelength = 0.532f, period = 0.35f;
    float thickness = 0.6f, pillar_n = 2.4f;
    int harmonics = 6, fill_samples = 18;
};

struct Results {
    double strehl = 0, fwhm = 0, dl = 0, encircled = 0, rms = 0, meanT = 0;
    double coverage_deg = 0, na = 0;
    int n_cells = 0, pillars = 0;
    PsfMap psf;
    std::vector<float> chrom_wl, chrom_focus;
    // Kept so the UI can export GDSII and run further analyses on demand.
    MetalensDesign design;
    UnitCellLibrary lib;
    Params used;
};

ImVec4 rgb(int r, int g, int b, float a = 1.0f) {
    return ImVec4(r / 255.0f, g / 255.0f, b / 255.0f, a);
}

// Utilitarian engineering-tool theme (think OpticStudio/MFC): light gray,
// square corners, thin borders everywhere, dense spacing, classic Windows-blue
// selection. Deliberately plain — it should read "serious instrument," not "app".
void apply_theme() {
    ImGui::StyleColorsLight();
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
    const ImVec4 face = rgb(238, 238, 238), field = rgb(255, 255, 255),
                 text = rgb(20, 20, 20), border = rgb(158, 158, 158),
                 head = rgb(222, 222, 222), sel = rgb(0, 120, 215),
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
}

std::atomic<bool> g_running{false};
std::atomic<bool> g_pending{false};  // worker finished, main thread must ingest
std::atomic<float> g_progress{0.0f};
std::mutex g_mtx;
std::string g_status = "Ready — set parameters and click Design.";
Results g_res;

// On-demand analyses (computed from the stored design when requested).
std::atomic<bool> g_tol_pending{false};
std::atomic<bool> g_fov_pending{false};
std::vector<ToleranceResult> g_tol;
std::vector<FieldPoint> g_fov;

void set_phase(const char* msg, float progress) {
    std::lock_guard<std::mutex> lk(g_mtx);
    g_status = msg;
    g_progress = progress;
}

void run_design(Params p) {
    g_running = true;
    set_phase("Building unit-cell library (RCWA sweep)...", 0.05f);
    const auto pillar = Material::constant(cdouble{p.pillar_n, 0.0}, "pillar");
    auto lib = build_unit_cell_library(pillar, materials::air(), materials::air(),
                                       materials::bk7(), p.period, p.wavelength,
                                       p.thickness, 0.08, 0.92, p.fill_samples,
                                       p.harmonics);
    set_phase("Assembling lens (phase profile -> pillar map)...", 0.45f);
    auto lens = design_metalens(lib, p.focal, p.diameter);
    set_phase("Analyzing focus (Rayleigh-Sommerfeld)...", 0.55f);
    auto foc = analyze_focus(lens, lib, p.focal, p.wavelength, p.diameter);
    double dl = p.wavelength * p.focal / p.diameter;
    set_phase("Rendering point-spread function...", 0.72f);
    auto psf = compute_psf(lens, lib, p.focal, p.wavelength, p.diameter, 161,
                           std::max(5.0 * dl, 4.0));
    set_phase("Chromatic sweep...", 0.9f);
    auto chrom = analyze_chromatic(lens, lib, p.focal, p.wavelength, p.diameter,
                                   p.wavelength * 0.85, p.wavelength * 1.25, 11);

    Results r;
    r.strehl = foc.strehl;
    r.fwhm = foc.fwhm_um;
    r.dl = foc.diffraction_limit_um;
    r.encircled = foc.encircled_energy;
    r.rms = lens.rms_phase_error_deg;
    r.meanT = lens.mean_amplitude;
    r.coverage_deg = lib.phase_span() * 180.0 / 3.14159265358979;
    r.na = std::sin(std::atan((p.diameter / 2.0) / p.focal));  // numerical aperture (in air)
    r.n_cells = lens.n_cells;
    r.pillars = lens.n_cells * lens.n_cells;
    r.psf = std::move(psf);
    for (auto& c : chrom) {
        r.chrom_wl.push_back(static_cast<float>(c.wavelength_um * 1000.0));
        r.chrom_focus.push_back(static_cast<float>(c.focal_length_um));
    }
    r.design = lens;
    r.lib = std::move(lib);
    r.used = p;
    {
        std::lock_guard<std::mutex> lk(g_mtx);
        g_res = std::move(r);
        g_status = "Done.";
        g_progress = 1.0f;
    }
    g_pending = true;
    g_running = false;
}

void run_tolerance() {
    g_running = true;
    set_phase("Monte-Carlo fabrication tolerance...", 0.2f);
    MetalensDesign d; UnitCellLibrary lib; Params p;
    { std::lock_guard<std::mutex> lk(g_mtx); d = g_res.design; lib = g_res.lib; p = g_res.used; }
    auto tol = analyze_tolerance(d, lib, p.focal, p.wavelength, p.diameter,
                                 {0.0, 5.0, 10.0, 20.0}, 12, 12345);
    { std::lock_guard<std::mutex> lk(g_mtx); g_tol = std::move(tol);
      g_status = "Tolerance analysis done."; g_progress = 1.0f; }
    g_tol_pending = true;
    g_running = false;
}

void run_fov() {
    g_running = true;
    set_phase("Field-of-view (off-axis) analysis...", 0.2f);
    MetalensDesign d; UnitCellLibrary lib; Params p;
    { std::lock_guard<std::mutex> lk(g_mtx); d = g_res.design; lib = g_res.lib; p = g_res.used; }
    auto fov = analyze_field_of_view(d, lib, p.focal, p.wavelength, p.diameter,
                                     {0.0, 1.0, 2.0, 5.0, 10.0});
    { std::lock_guard<std::mutex> lk(g_mtx); g_fov = std::move(fov);
      g_status = "Field-of-view analysis done."; g_progress = 1.0f; }
    g_fov_pending = true;
    g_running = false;
}

// Upload a PSF map to an OpenGL texture (grayscale, gamma-boosted). Runs on the
// main/GL thread.
void upload_psf_texture(const PsfMap& psf, unsigned int& tex) {
    if (psf.n <= 0 || psf.intensity.empty()) return;
    double mx = 0.0;
    for (double v : psf.intensity) mx = std::max(mx, v);
    if (mx <= 0.0) mx = 1.0;
    std::vector<std::uint8_t> rgba(static_cast<std::size_t>(psf.n) * psf.n * 4);
    for (std::size_t i = 0; i < psf.intensity.size(); ++i) {
        double v = std::pow(psf.intensity[i] / mx, 1.0 / 2.2);  // gamma
        auto b = static_cast<std::uint8_t>(std::lround(std::clamp(v, 0.0, 1.0) * 255.0));
        rgba[i * 4 + 0] = b;
        rgba[i * 4 + 1] = b;
        rgba[i * 4 + 2] = b;
        rgba[i * 4 + 3] = 255;
    }
    if (tex == 0) glGenTextures(1, &tex);
    glBindTexture(GL_TEXTURE_2D, tex);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, psf.n, psf.n, 0, GL_RGBA,
                 GL_UNSIGNED_BYTE, rgba.data());
}

} // namespace

int main() {
    if (!glfwInit()) return 1;
    GLFWwindow* window = glfwCreateWindow(1100, 720, "CELERIS — Metalens Designer",
                                          nullptr, nullptr);
    if (!window) { glfwTerminate(); return 1; }
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO();
    // Use the system Segoe UI font (present on all Windows) for a real-software
    // look instead of the pixelated default. Falls back to the default if absent.
    {
        const char* font_path = "C:\\Windows\\Fonts\\segoeui.ttf";
        std::ifstream probe(font_path);
        if (probe.good()) io.Fonts->AddFontFromFileTTF(font_path, 19.0f);
    }
    apply_theme();
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 130");

    Params params;
    unsigned int psf_tex = 0;
    bool have_result = false, have_tol = false, have_fov = false;
    char gds_name[256] = "metalens.gds";
    double run_start = 0.0;  // ImGui time when the current run was launched

    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();

        if (g_pending.exchange(false)) {
            std::lock_guard<std::mutex> lk(g_mtx);
            upload_psf_texture(g_res.psf, psf_tex);
            have_result = true;
        }
        if (g_tol_pending.exchange(false)) have_tol = true;
        if (g_fov_pending.exchange(false)) have_fov = true;

        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        bool running = g_running.load();
        auto launch_design = [&] {
            if (!g_running.load()) {
                run_start = ImGui::GetTime();
                std::thread(run_design, params).detach();
            }
        };

        ImGui::SetNextWindowPos({0, 0});
        ImGui::SetNextWindowSize(io.DisplaySize);
        ImGui::Begin("CELERIS", nullptr,
                     ImGuiWindowFlags_MenuBar | ImGuiWindowFlags_NoResize |
                         ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoCollapse |
                         ImGuiWindowFlags_NoTitleBar |
                         ImGuiWindowFlags_NoBringToFrontOnFocus);

        static bool show_about = false;
        if (ImGui::BeginMenuBar()) {
            if (ImGui::BeginMenu("File")) {
                if (ImGui::MenuItem("Exit")) glfwSetWindowShouldClose(window, 1);
                ImGui::EndMenu();
            }
            if (ImGui::BeginMenu("Run")) {
                if (ImGui::MenuItem("Design Lens", "F5", false, !running))
                    launch_design();
                ImGui::EndMenu();
            }
            if (ImGui::BeginMenu("Help")) {
                if (ImGui::MenuItem("About CELERIS")) show_about = true;
                ImGui::EndMenu();
            }
            ImGui::EndMenuBar();
        }
        if (ImGui::IsKeyPressed(ImGuiKey_F5)) launch_design();
        if (show_about) { ImGui::OpenPopup("About CELERIS"); show_about = false; }
        if (ImGui::BeginPopupModal("About CELERIS", nullptr,
                                   ImGuiWindowFlags_AlwaysAutoResize)) {
            ImGui::TextUnformatted(
                "CELERIS  —  metalens design via rigorous coupled-wave analysis.");
            ImGui::TextUnformatted(
                "RCWA engine, inverse design, GDSII export, focal/chromatic analysis.");
            ImGui::Separator();
            if (ImGui::Button("Close", ImVec2(120, 0))) ImGui::CloseCurrentPopup();
            ImGui::EndPopup();
        }

        const float status_h = ImGui::GetFrameHeight() + 10.0f;
        const float body_h = ImGui::GetContentRegionAvail().y - status_h;

        // ---- Lens Data panel ----
        ImGui::BeginChild("lensdata", ImVec2(290, body_h), true);
        ImGui::TextUnformatted("Lens Data");
        ImGui::Separator();
        auto field = [](const char* label) {
            ImGui::TableNextRow();
            ImGui::TableSetColumnIndex(0);
            ImGui::AlignTextToFramePadding();
            ImGui::TextUnformatted(label);
            ImGui::TableSetColumnIndex(1);
            ImGui::SetNextItemWidth(-FLT_MIN);
        };
        if (ImGui::BeginTable("params", 2,
                              ImGuiTableFlags_BordersInner |
                                  ImGuiTableFlags_SizingStretchProp)) {
            field("Focal length (um)"); ImGui::InputFloat("##f", &params.focal, 0, 0, "%.1f");
            field("Aperture (um)");     ImGui::InputFloat("##d", &params.diameter, 0, 0, "%.1f");
            field("Wavelength (um)");   ImGui::InputFloat("##w", &params.wavelength, 0, 0, "%.3f");
            field("Period (um)");       ImGui::InputFloat("##p", &params.period, 0, 0, "%.3f");
            field("Pillar height (um)");ImGui::InputFloat("##h", &params.thickness, 0, 0, "%.2f");
            field("Pillar index n");    ImGui::InputFloat("##n", &params.pillar_n, 0, 0, "%.2f");
            field("RCWA harmonics");    ImGui::InputInt("##m", &params.harmonics);
            field("Library samples");   ImGui::InputInt("##s", &params.fill_samples);
            ImGui::EndTable();
        }
        ImGui::Spacing();
        if (running) ImGui::BeginDisabled();
        if (ImGui::Button("Run Design  (F5)", ImVec2(-FLT_MIN, 30))) launch_design();
        if (running) ImGui::EndDisabled();

        if (running) {
            ImGui::Spacing();
            char ov[48];
            std::snprintf(ov, sizeof(ov), "%.0f%%", g_progress.load() * 100.0f);
            ImGui::ProgressBar(g_progress.load(), ImVec2(-FLT_MIN, 0), ov);
            ImGui::Text("Elapsed: %.1f s", ImGui::GetTime() - run_start);
            std::lock_guard<std::mutex> lk(g_mtx);
            ImGui::TextWrapped("%s", g_status.c_str());
        }

        ImGui::Spacing();
        ImGui::Separator();
        ImGui::TextUnformatted("Output / Analysis");
        const bool can_act = have_result && !running;
        if (!can_act) ImGui::BeginDisabled();
        ImGui::SetNextItemWidth(-FLT_MIN);
        ImGui::InputText("##gds", gds_name, sizeof(gds_name));
        if (ImGui::Button("Save GDSII", ImVec2(-FLT_MIN, 0))) {
            MetalensDesign d;
            { std::lock_guard<std::mutex> lk(g_mtx); d = g_res.design; }
            int np = write_metalens_gds(d, gds_name);
            std::lock_guard<std::mutex> lk(g_mtx);
            g_status = np >= 0 ? std::format("Wrote {} pillars -> {}", np, gds_name)
                               : std::string("ERROR: GDSII write failed");
        }
        if (ImGui::Button("Run Tolerance", ImVec2(-FLT_MIN, 0)))
            std::thread(run_tolerance).detach();
        if (ImGui::Button("Run Field of View", ImVec2(-FLT_MIN, 0)))
            std::thread(run_fov).detach();
        if (!can_act) ImGui::EndDisabled();
        ImGui::EndChild();

        // ---- Results panel ----
        ImGui::SameLine();
        ImGui::BeginChild("results", ImVec2(0, body_h), true);
        if (!have_result) {
            ImGui::TextDisabled("No results. Enter lens data and Run Design (F5).");
        } else {
            std::lock_guard<std::mutex> lk(g_mtx);
            const auto kv_table = [](const char* id) {
                return ImGui::BeginTable(
                    id, 2, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg);
            };
            auto kvrow = [](const char* k, const std::string& v) {
                ImGui::TableNextRow();
                ImGui::TableNextColumn(); ImGui::TextUnformatted(k);
                ImGui::TableNextColumn(); ImGui::TextUnformatted(v.c_str());
            };

            ImGui::TextUnformatted("Design Summary");
            if (kv_table("sum")) {
                kvrow("Array", std::format("{0} x {0} pillars ({1})", g_res.n_cells,
                                           g_res.pillars));
                ImGui::TableNextRow();
                ImGui::TableNextColumn(); ImGui::TextUnformatted("Numerical aperture");
                ImGui::TableNextColumn();
                if (g_res.na >= 0.7)
                    ImGui::TextColored(rgb(190, 40, 40),
                                       "%.2f  (extreme - metrics approximate)", g_res.na);
                else
                    ImGui::Text("%.2f", g_res.na);
                kvrow("Phase coverage", std::format("{:.0f} deg", g_res.coverage_deg));
                kvrow("RMS phase error", std::format("{:.1f} deg", g_res.rms));
                kvrow("Mean transmission", std::format("{:.3f}", g_res.meanT));
                ImGui::EndTable();
            }
            ImGui::Spacing();
            ImGui::TextUnformatted("Focus Performance");
            if (kv_table("foc")) {
                ImGui::TableNextRow();
                ImGui::TableNextColumn(); ImGui::TextUnformatted("Strehl ratio");
                ImGui::TableNextColumn();
                ImVec4 sc = g_res.strehl >= 0.8   ? rgb(20, 130, 40)
                            : g_res.strehl >= 0.5 ? rgb(170, 120, 0)
                                                  : rgb(190, 40, 40);
                ImGui::TextColored(sc, "%.3f", g_res.strehl);
                kvrow("Spot FWHM", std::format("{:.2f} um", g_res.fwhm));
                kvrow("Diffraction limit", std::format("{:.2f} um", g_res.dl));
                kvrow("Encircled energy", std::format("{:.0f} %", g_res.encircled * 100.0));
                ImGui::EndTable();
            }
            ImGui::Spacing();
            ImGui::Columns(2, "viz", false);
            ImGui::TextUnformatted("Focal PSF");
            if (psf_tex)
                ImGui::Image((ImTextureID)(intptr_t)psf_tex, ImVec2(280, 280));
            ImGui::NextColumn();
            ImGui::TextUnformatted("Chromatic focus shift (um vs wavelength)");
            if (!g_res.chrom_focus.empty())
                ImGui::PlotLines("##chrom", g_res.chrom_focus.data(),
                                 static_cast<int>(g_res.chrom_focus.size()), 0,
                                 nullptr, FLT_MAX, FLT_MAX, ImVec2(-FLT_MIN, 200));
            ImGui::Columns(1);

            if (have_tol && !g_tol.empty()) {
                ImGui::Spacing();
                ImGui::TextUnformatted("Fabrication Tolerance (Strehl vs CD error)");
                if (ImGui::BeginTable("tol", 4,
                                      ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg)) {
                    ImGui::TableSetupColumn("sigma (nm)");
                    ImGui::TableSetupColumn("mean");
                    ImGui::TableSetupColumn("std");
                    ImGui::TableSetupColumn("worst");
                    ImGui::TableHeadersRow();
                    for (auto& t : g_tol) {
                        ImGui::TableNextRow();
                        ImGui::TableNextColumn(); ImGui::Text("%.0f", t.sigma_nm);
                        ImGui::TableNextColumn(); ImGui::Text("%.3f", t.mean_strehl);
                        ImGui::TableNextColumn(); ImGui::Text("%.3f", t.std_strehl);
                        ImGui::TableNextColumn(); ImGui::Text("%.3f", t.worst_strehl);
                    }
                    ImGui::EndTable();
                }
            }
            if (have_fov && !g_fov.empty()) {
                ImGui::Spacing();
                ImGui::TextUnformatted("Field of View (off-axis)");
                if (ImGui::BeginTable("fov", 3,
                                      ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg)) {
                    ImGui::TableSetupColumn("angle (deg)");
                    ImGui::TableSetupColumn("rel. Strehl");
                    ImGui::TableSetupColumn("shift (um)");
                    ImGui::TableHeadersRow();
                    for (auto& f : g_fov) {
                        ImGui::TableNextRow();
                        ImGui::TableNextColumn(); ImGui::Text("%.0f", f.angle_deg);
                        ImGui::TableNextColumn(); ImGui::Text("%.3f", f.rel_strehl);
                        ImGui::TableNextColumn(); ImGui::Text("%.2f", f.spot_shift_um);
                    }
                    ImGui::EndTable();
                }
            }
        }
        ImGui::EndChild();

        // ---- Status bar ----
        ImGui::BeginChild("statusbar", ImVec2(0, 0), true);
        ImGui::AlignTextToFramePadding();
        {
            std::lock_guard<std::mutex> lk(g_mtx);
            ImGui::Text("%s%s", running ? "[BUSY]  " : "", g_status.c_str());
        }
        ImGui::EndChild();

        ImGui::End();

        ImGui::Render();
        int w, h;
        glfwGetFramebufferSize(window, &w, &h);
        glViewport(0, 0, w, h);
        glClearColor(0.10f, 0.10f, 0.12f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        glfwSwapBuffers(window);
    }

    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}
