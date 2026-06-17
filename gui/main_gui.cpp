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
#include "imgui_internal.h"  // DockBuilder for the default layout

#include <algorithm>
#include <atomic>
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <format>
#include <fstream>
#include <mutex>
#include <optional>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include "celeris/analysis/chromatic.hpp"
#include "celeris/analysis/field.hpp"
#include "celeris/analysis/focal.hpp"
#include "celeris/analysis/mtf.hpp"
#include "celeris/analysis/throughfocus.hpp"
#include "celeris/analysis/tolerance.hpp"
#include "celeris/analysis/wavefront.hpp"
#include "celeris/design/metalens.hpp"
#include "celeris/design/system_opt.hpp"
#include "celeris/io/gds.hpp"
#include "celeris/io/material_csv.hpp"
#include "celeris/materials/database.hpp"

using namespace celeris;

namespace {

struct LayerRow {  // an extra (unpatterned by default) layer above the pillars
    float n = 1.46f;       // refractive index
    float fill = 1.0f;     // 1.0 = uniform film; <1 = a square patch
    float thickness = 0.1f;
};

struct Params {
    float focal = 50, diameter = 20, wavelength = 0.532f, period = 0.35f;
    float thickness = 0.6f, pillar_n = 2.4f;  // the active (patterned) layer
    int harmonics = 6, fill_samples = 18;
    int pillar_mat = 3;     // index into kPillarMats (default: TiO2 approx)
    int substrate_mat = 0;  // 0 = N-BK7, 1 = air, 2 = fused silica
    std::vector<LayerRow> extra_layers;  // stacked above the pillars (cap/AR…)
};

const char* kPillarMats[] = {"Custom (constant n)", "Silicon nitride (Si3N4)",
                             "Fused silica (SiO2)",  "TiO2 (approx n=2.40)",
                             "a-Si (approx n=3.50)", "GaN (approx n=2.35)",
                             "Loaded CSV"};
const char* kSubstrates[] = {"N-BK7", "Air", "Fused silica (SiO2)"};

// A material loaded from a CSV file at runtime (real n,k data).
std::optional<Material> g_loaded_material;
std::string g_loaded_name = "(none)";

Material make_pillar(const Params& p) {
    switch (p.pillar_mat) {
        case 1: return materials::silicon_nitride();
        case 2: return materials::fused_silica();
        case 3: return Material::constant(cdouble{2.40, 0.0}, "TiO2~");
        case 4: return Material::constant(cdouble{3.50, 0.0}, "a-Si~");
        case 5: return Material::constant(cdouble{2.35, 0.0}, "GaN~");
        case 6:
            if (g_loaded_material) return *g_loaded_material;
            return Material::constant(cdouble{p.pillar_n, 0.0}, "custom");
        default: return Material::constant(cdouble{p.pillar_n, 0.0}, "custom");
    }
}

const Material& make_substrate(const Params& p) {
    switch (p.substrate_mat) {
        case 1: return materials::air();
        case 2: return materials::fused_silica();
        default: return materials::bk7();
    }
}

struct Results {
    double strehl = 0, fwhm = 0, dl = 0, encircled = 0, rms = 0, meanT = 0;
    double coverage_deg = 0, na = 0;
    int n_cells = 0, pillars = 0;
    PsfMap psf;
    std::vector<float> chrom_wl, chrom_focus;
    WavefrontAnalysis wf;
    MtfCurve mtf;
    ThroughFocus tf;
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

// Design optimizer result (best period/height applied to the UI on completion).
std::atomic<bool> g_opt_pending{false};
SystemOptResult g_opt;

void set_phase(const char* msg, float progress) {
    std::lock_guard<std::mutex> lk(g_mtx);
    g_status = msg;
    g_progress = progress;
}

void run_design(Params p) {
    g_running = true;
    set_phase("Building unit-cell library (RCWA sweep)...", 0.05f);
    // Assemble the unit-cell stack: extra layers (caps/AR coatings) above the
    // patterned pillar layer, which is the active (fill-swept) layer.
    Rcwa2DStack stack;
    stack.period_x_um = stack.period_y_um = p.period;
    for (const auto& L : p.extra_layers)
        stack.layers.push_back(RectCell2D{
            Material::constant(cdouble{L.n, 0.0}, "layer"), materials::air(),
            L.fill, L.fill, L.thickness});
    const int active = static_cast<int>(stack.layers.size());
    stack.layers.push_back(RectCell2D{make_pillar(p), materials::air(),
                                      0.5, 0.5, p.thickness});
    auto lib = build_unit_cell_library_stack(stack, active, materials::air(),
                                             make_substrate(p), p.wavelength, 0.08,
                                             0.92, p.fill_samples, p.harmonics);
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
    r.wf = analyze_wavefront(lens, lib, p.focal, p.wavelength, p.diameter);
    r.mtf = analyze_mtf(lens, lib, p.focal, p.wavelength, p.diameter);
    r.tf = analyze_through_focus(lens, lib, p.focal, p.wavelength, p.diameter);
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

void run_optimize(Params p) {
    g_running = true;
    set_phase("Optimizing design (period x height search)...", 0.0f);
    auto res = optimize_system(
        make_pillar(p), materials::air(), materials::air(), make_substrate(p), p.focal,
        p.diameter, p.wavelength, 0.20, 0.45, 0.30, 1.00, /*grid=*/5, /*M=*/5,
        /*fill_samples=*/12, /*efficiency_weight=*/0.3,
        [](float fr) { set_phase("Optimizing design (period x height search)...", fr); });
    {
        std::lock_guard<std::mutex> lk(g_mtx);
        g_opt = res;
        g_status = std::format("Optimized: period {:.3f} um, height {:.3f} um "
                               "(Strehl~{:.3f}). Applied -- press F5 to run.",
                               res.period_um, res.thickness_um, res.strehl);
        g_progress = 1.0f;
    }
    g_opt_pending = true;
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

// Upload the wavefront OPD map with a diverging blue-white-red colormap
// (negative -> blue, 0 -> white, positive -> red), outside-aperture -> gray.
void upload_wavefront_texture(const WavefrontAnalysis& wf, unsigned int& tex) {
    if (wf.n <= 0 || wf.opd.empty()) return;
    double amp = 1e-9;
    for (std::size_t i = 0; i < wf.opd.size(); ++i)
        if (wf.mask[i]) amp = std::max(amp, std::abs(wf.opd[i]));
    std::vector<std::uint8_t> rgba(static_cast<std::size_t>(wf.n) * wf.n * 4);
    for (std::size_t i = 0; i < wf.opd.size(); ++i) {
        std::uint8_t r, g, b;
        if (!wf.mask[i]) { r = g = b = 70; }  // outside aperture
        else {
            double v = std::clamp(wf.opd[i] / amp, -1.0, 1.0);  // -1..1
            if (v >= 0) { r = 255; g = b = static_cast<std::uint8_t>(255 * (1 - v)); }
            else        { b = 255; r = g = static_cast<std::uint8_t>(255 * (1 + v)); }
        }
        rgba[i * 4 + 0] = r; rgba[i * 4 + 1] = g; rgba[i * 4 + 2] = b; rgba[i * 4 + 3] = 255;
    }
    if (tex == 0) glGenTextures(1, &tex);
    glBindTexture(GL_TEXTURE_2D, tex);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, wf.n, wf.n, 0, GL_RGBA,
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
    io.ConfigFlags |= ImGuiConfigFlags_DockingEnable;  // dockable/tabbable workspace
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
    unsigned int psf_tex = 0, wf_tex = 0;
    bool have_result = false, have_tol = false, have_fov = false;
    char gds_name[256] = "metalens.gds";
    double run_start = 0.0;  // ImGui time when the current run was launched

    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();

        if (g_pending.exchange(false)) {
            std::lock_guard<std::mutex> lk(g_mtx);
            upload_psf_texture(g_res.psf, psf_tex);
            upload_wavefront_texture(g_res.wf, wf_tex);
            have_result = true;
        }
        if (g_tol_pending.exchange(false)) have_tol = true;
        if (g_fov_pending.exchange(false)) have_fov = true;
        if (g_opt_pending.exchange(false)) {
            std::lock_guard<std::mutex> lk(g_mtx);
            params.period = static_cast<float>(g_opt.period_um);
            params.thickness = static_cast<float>(g_opt.thickness_um);
        }

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

        static bool show_about = false;
        static bool win_lens = true, win_sum = true, win_foc = true, win_psf = true,
                    win_chr = true, win_tol = true, win_fov = true, win_log = true,
                    win_wf = true, win_mtf = true, win_tf = true, win_stack = true,
                    win_mats = true;
        const bool can_act = have_result && !running;

        if (ImGui::BeginMainMenuBar()) {
            if (ImGui::BeginMenu("File")) {
                if (ImGui::MenuItem("Exit")) glfwSetWindowShouldClose(window, 1);
                ImGui::EndMenu();
            }
            if (ImGui::BeginMenu("Run")) {
                if (ImGui::MenuItem("Design Lens", "F5", false, !running)) launch_design();
                ImGui::Separator();
                if (ImGui::MenuItem("Tolerance Analysis", nullptr, false, can_act))
                    std::thread(run_tolerance).detach();
                if (ImGui::MenuItem("Field of View", nullptr, false, can_act))
                    std::thread(run_fov).detach();
                ImGui::EndMenu();
            }
            if (ImGui::BeginMenu("View")) {
                ImGui::MenuItem("Lens Data", nullptr, &win_lens);
                ImGui::MenuItem("Layer Stack", nullptr, &win_stack);
                ImGui::MenuItem("Materials", nullptr, &win_mats);
                ImGui::MenuItem("Design Summary", nullptr, &win_sum);
                ImGui::MenuItem("Focus Performance", nullptr, &win_foc);
                ImGui::MenuItem("Focal PSF", nullptr, &win_psf);
                ImGui::MenuItem("Wavefront", nullptr, &win_wf);
                ImGui::MenuItem("MTF", nullptr, &win_mtf);
                ImGui::MenuItem("Through Focus", nullptr, &win_tf);
                ImGui::MenuItem("Chromatic", nullptr, &win_chr);
                ImGui::MenuItem("Tolerance", nullptr, &win_tol);
                ImGui::MenuItem("Field of View", nullptr, &win_fov);
                ImGui::MenuItem("Log", nullptr, &win_log);
                ImGui::EndMenu();
            }
            if (ImGui::BeginMenu("Help")) {
                if (ImGui::MenuItem("About CELERIS")) show_about = true;
                ImGui::EndMenu();
            }
            ImGui::EndMainMenuBar();
        }
        if (ImGui::IsKeyPressed(ImGuiKey_F5)) launch_design();
        if (show_about) { ImGui::OpenPopup("About CELERIS"); show_about = false; }
        if (ImGui::BeginPopupModal("About CELERIS", nullptr,
                                   ImGuiWindowFlags_AlwaysAutoResize)) {
            ImGui::TextUnformatted(
                "CELERIS  -  metalens design via rigorous coupled-wave analysis.");
            ImGui::TextUnformatted(
                "RCWA engine, inverse design, GDSII export, focal/chromatic analysis.");
            ImGui::Separator();
            if (ImGui::Button("Close", ImVec2(120, 0))) ImGui::CloseCurrentPopup();
            ImGui::EndPopup();
        }

        // Dockable workspace, with a sensible default layout on first launch.
        ImGuiID dockid = ImGui::DockSpaceOverViewport(0, ImGui::GetMainViewport());
        static bool first_layout = true;
        if (first_layout) {
            first_layout = false;
            ImGui::DockBuilderRemoveNode(dockid);
            ImGui::DockBuilderAddNode(dockid, ImGuiDockNodeFlags_DockSpace);
            ImGui::DockBuilderSetNodeSize(dockid, ImGui::GetMainViewport()->WorkSize);
            ImGuiID center;
            ImGuiID left = ImGui::DockBuilderSplitNode(dockid, ImGuiDir_Left, 0.24f, nullptr, &center);
            ImGuiID cbottom;
            ImGuiID ctop = ImGui::DockBuilderSplitNode(center, ImGuiDir_Up, 0.62f, nullptr, &cbottom);
            ImGuiID ctr;
            ImGuiID ctl = ImGui::DockBuilderSplitNode(ctop, ImGuiDir_Left, 0.5f, nullptr, &ctr);
            ImGui::DockBuilderDockWindow("Lens Data", left);
            ImGui::DockBuilderDockWindow("Layer Stack", left);
            ImGui::DockBuilderDockWindow("Materials", left);
            ImGui::DockBuilderDockWindow("Design Summary", ctl);
            ImGui::DockBuilderDockWindow("Focus Performance", ctl);
            ImGui::DockBuilderDockWindow("Focal PSF", ctr);
            ImGui::DockBuilderDockWindow("Wavefront", ctr);
            ImGui::DockBuilderDockWindow("Chromatic", ctr);
            ImGui::DockBuilderDockWindow("MTF", cbottom);
            ImGui::DockBuilderDockWindow("Through Focus", cbottom);
            ImGui::DockBuilderDockWindow("Tolerance", cbottom);
            ImGui::DockBuilderDockWindow("Field of View", cbottom);
            ImGui::DockBuilderDockWindow("Log", cbottom);
            ImGui::DockBuilderFinish(dockid);
        }

        auto kvrow = [](const char* k, const std::string& v) {
            ImGui::TableNextRow();
            ImGui::TableNextColumn(); ImGui::TextUnformatted(k);
            ImGui::TableNextColumn(); ImGui::TextUnformatted(v.c_str());
        };

        // ---- Lens Data ----
        if (win_lens) {
            if (ImGui::Begin("Lens Data", &win_lens)) {
                auto fld = [](const char* label) {
                    ImGui::TableNextRow(); ImGui::TableSetColumnIndex(0);
                    ImGui::AlignTextToFramePadding(); ImGui::TextUnformatted(label);
                    ImGui::TableSetColumnIndex(1); ImGui::SetNextItemWidth(-FLT_MIN);
                };
                if (ImGui::BeginTable("params", 2,
                                      ImGuiTableFlags_BordersInner | ImGuiTableFlags_SizingStretchProp)) {
                    fld("Focal length (um)"); ImGui::InputFloat("##f", &params.focal, 0, 0, "%.1f");
                    fld("Aperture (um)");     ImGui::InputFloat("##d", &params.diameter, 0, 0, "%.1f");
                    fld("Wavelength (um)");   ImGui::InputFloat("##w", &params.wavelength, 0, 0, "%.3f");
                    fld("Period (um)");       ImGui::InputFloat("##p", &params.period, 0, 0, "%.3f");
                    fld("Pillar height (um)");ImGui::InputFloat("##h", &params.thickness, 0, 0, "%.2f");
                    fld("Pillar material");   ImGui::Combo("##pmat", &params.pillar_mat, kPillarMats, IM_ARRAYSIZE(kPillarMats));
                    fld("  index n (custom)");ImGui::InputFloat("##n", &params.pillar_n, 0, 0, "%.2f");
                    fld("Substrate");         ImGui::Combo("##sub", &params.substrate_mat, kSubstrates, IM_ARRAYSIZE(kSubstrates));
                    fld("RCWA harmonics");    ImGui::InputInt("##m", &params.harmonics);
                    fld("Library samples");   ImGui::InputInt("##s", &params.fill_samples);
                    ImGui::EndTable();
                }
                ImGui::TextDisabled("pillar n(%.3fum) = %.3f + %.3fi", params.wavelength,
                                    make_pillar(params).index(params.wavelength).real(),
                                    make_pillar(params).index(params.wavelength).imag());
                ImGui::Spacing();
                if (running) ImGui::BeginDisabled();
                if (ImGui::Button("Run Design  (F5)", ImVec2(-FLT_MIN, 30))) launch_design();
                if (ImGui::Button("Optimize Design (period x height)", ImVec2(-FLT_MIN, 0)))
                    std::thread(run_optimize, params).detach();
                if (running) ImGui::EndDisabled();
                if (running) {
                    char ov[48];
                    std::snprintf(ov, sizeof(ov), "%.0f%%", g_progress.load() * 100.0f);
                    ImGui::ProgressBar(g_progress.load(), ImVec2(-FLT_MIN, 0), ov);
                    ImGui::Text("Elapsed: %.1f s", ImGui::GetTime() - run_start);
                }
                ImGui::Spacing();
                ImGui::SeparatorText("Output / Analysis");
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
            }
            ImGui::End();
        }

        // ---- Layer Stack (editable system model) ----
        if (win_stack) {
            if (ImGui::Begin("Layer Stack", &win_stack)) {
                ImGui::TextWrapped(
                    "Unit-cell stack, top to bottom. Extra layers (caps / AR "
                    "coatings) sit above the patterned pillar layer.");
                ImGui::Spacing();
                ImGui::TextUnformatted("  incident: air");
                int del = -1;
                for (int i = 0; i < static_cast<int>(params.extra_layers.size()); ++i) {
                    ImGui::PushID(i);
                    LayerRow& L = params.extra_layers[i];
                    ImGui::Text("Layer %d", i + 1); ImGui::SameLine();
                    ImGui::SetNextItemWidth(70); ImGui::InputFloat("n", &L.n, 0, 0, "%.2f"); ImGui::SameLine();
                    ImGui::SetNextItemWidth(70); ImGui::InputFloat("fill", &L.fill, 0, 0, "%.2f"); ImGui::SameLine();
                    ImGui::SetNextItemWidth(80); ImGui::InputFloat("h(um)", &L.thickness, 0, 0, "%.2f"); ImGui::SameLine();
                    if (ImGui::SmallButton("x")) del = i;
                    ImGui::PopID();
                }
                if (del >= 0) params.extra_layers.erase(params.extra_layers.begin() + del);
                if (ImGui::Button("Add layer"))
                    params.extra_layers.push_back(LayerRow{});
                ImGui::Separator();
                ImGui::Text("  >> patterned pillars: n=%.2f  h=%.2f um  (fill swept)",
                            params.pillar_n, params.thickness);
                ImGui::TextUnformatted("  substrate: N-BK7");
                ImGui::Spacing();
                ImGui::TextDisabled("fill = 1.0 -> uniform film; <1 -> square patch. "
                                    "Press F5 to rebuild.");
            }
            ImGui::End();
        }

        // ---- Materials (dispersion + CSV import) ----
        if (win_mats) {
            if (ImGui::Begin("Materials", &win_mats)) {
                ImGui::Text("Pillar material: %s", kPillarMats[params.pillar_mat]);
                Material m = make_pillar(params);
                const float l0 = 0.40f, l1 = 0.75f;
                const int NS = 60;
                std::vector<float> nv(NS), kv(NS);
                for (int i = 0; i < NS; ++i) {
                    double L = l0 + (l1 - l0) * i / (NS - 1);
                    auto idx = m.index(L);
                    nv[i] = static_cast<float>(idx.real());
                    kv[i] = static_cast<float>(idx.imag());
                }
                ImGui::PlotLines("##nl", nv.data(), NS, 0, "n  (0.40-0.75 um)",
                                 FLT_MAX, FLT_MAX, ImVec2(-FLT_MIN, 90));
                ImGui::PlotLines("##kl", kv.data(), NS, 0, "k  (extinction)",
                                 FLT_MAX, FLT_MAX, ImVec2(-FLT_MIN, 60));

                ImGui::SeparatorText("Built-in materials (n at design wavelength)");
                if (ImGui::BeginTable("matlib", 2,
                                      ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg)) {
                    ImGui::TableSetupColumn("material");
                    ImGui::TableSetupColumn("n(lambda)");
                    ImGui::TableHeadersRow();
                    auto row = [&](const char* nm, const Material& mat) {
                        ImGui::TableNextRow();
                        ImGui::TableNextColumn(); ImGui::TextUnformatted(nm);
                        ImGui::TableNextColumn();
                        ImGui::Text("%.3f", mat.index(params.wavelength).real());
                    };
                    row("Air", materials::air());
                    row("N-BK7", materials::bk7());
                    row("Fused silica (SiO2)", materials::fused_silica());
                    row("Silicon nitride (Si3N4)", materials::silicon_nitride());
                    ImGui::EndTable();
                }

                ImGui::SeparatorText("Load real n,k from CSV");
                ImGui::TextDisabled("format: wavelength_um, n [, k]  (refractiveindex.info)");
                static char path[256] = "material.csv";
                ImGui::SetNextItemWidth(-FLT_MIN);
                ImGui::InputText("##mcsv", path, sizeof(path));
                if (ImGui::Button("Load as pillar material", ImVec2(-FLT_MIN, 0))) {
                    try {
                        Material lm = load_material_csv(path, "CSV");
                        g_loaded_material = std::move(lm);
                        g_loaded_name = path;
                        params.pillar_mat = 6;
                        std::lock_guard<std::mutex> lk(g_mtx);
                        g_status = std::format("Loaded material from {}", path);
                    } catch (const std::exception& e) {
                        std::lock_guard<std::mutex> lk(g_mtx);
                        g_status = std::string("CSV load failed: ") + e.what();
                    }
                }
                ImGui::Text("Loaded CSV: %s", g_loaded_name.c_str());
            }
            ImGui::End();
        }

        // ---- Design Summary ----
        if (win_sum) {
            if (ImGui::Begin("Design Summary", &win_sum)) {
                if (!have_result) ImGui::TextDisabled("Run a design (F5).");
                else {
                    std::lock_guard<std::mutex> lk(g_mtx);
                    if (ImGui::BeginTable("sum", 2, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg)) {
                        kvrow("Array", std::format("{0} x {0} pillars ({1})", g_res.n_cells, g_res.pillars));
                        ImGui::TableNextRow();
                        ImGui::TableNextColumn(); ImGui::TextUnformatted("Numerical aperture");
                        ImGui::TableNextColumn();
                        if (g_res.na >= 0.7)
                            ImGui::TextColored(rgb(190, 40, 40), "%.2f  (extreme)", g_res.na);
                        else ImGui::Text("%.2f", g_res.na);
                        kvrow("Phase coverage", std::format("{:.0f} deg", g_res.coverage_deg));
                        kvrow("RMS phase error", std::format("{:.1f} deg", g_res.rms));
                        kvrow("Mean transmission", std::format("{:.3f}", g_res.meanT));
                        ImGui::EndTable();
                    }
                }
            }
            ImGui::End();
        }

        // ---- Focus Performance ----
        if (win_foc) {
            if (ImGui::Begin("Focus Performance", &win_foc)) {
                if (!have_result) ImGui::TextDisabled("Run a design (F5).");
                else {
                    std::lock_guard<std::mutex> lk(g_mtx);
                    if (ImGui::BeginTable("foc", 2, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg)) {
                        ImGui::TableNextRow();
                        ImGui::TableNextColumn(); ImGui::TextUnformatted("Strehl ratio");
                        ImGui::TableNextColumn();
                        ImVec4 sc = g_res.strehl >= 0.8 ? rgb(20, 130, 40)
                                    : g_res.strehl >= 0.5 ? rgb(170, 120, 0) : rgb(190, 40, 40);
                        ImGui::TextColored(sc, "%.3f", g_res.strehl);
                        kvrow("Spot FWHM", std::format("{:.2f} um", g_res.fwhm));
                        kvrow("Diffraction limit", std::format("{:.2f} um", g_res.dl));
                        kvrow("Encircled energy", std::format("{:.0f} %", g_res.encircled * 100.0));
                        ImGui::EndTable();
                    }
                }
            }
            ImGui::End();
        }

        // ---- Focal PSF ----
        if (win_psf) {
            if (ImGui::Begin("Focal PSF", &win_psf)) {
                if (psf_tex) {
                    ImVec2 a = ImGui::GetContentRegionAvail();
                    float s = std::max(64.0f, std::min(a.x, a.y));
                    ImGui::Image((ImTextureID)(intptr_t)psf_tex, ImVec2(s, s));
                } else ImGui::TextDisabled("Run a design (F5).");
            }
            ImGui::End();
        }

        // ---- Wavefront ----
        if (win_wf) {
            if (ImGui::Begin("Wavefront", &win_wf)) {
                if (!have_result) ImGui::TextDisabled("Run a design (F5).");
                else {
                    std::lock_guard<std::mutex> lk(g_mtx);
                    ImGui::Text("RMS %.4f wv   P-V %.3f wv   Strehl(Marechal) %.3f",
                                g_res.wf.rms_waves, g_res.wf.pv_waves,
                                g_res.wf.strehl_marechal);
                    if (wf_tex) {
                        ImGui::Image((ImTextureID)(intptr_t)wf_tex, ImVec2(180, 180));
                        ImGui::SameLine();
                    }
                    if (ImGui::BeginTable("zern", 2,
                                          ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg |
                                              ImGuiTableFlags_ScrollY,
                                          ImVec2(0, 180))) {
                        ImGui::TableSetupColumn("Zernike aberration");
                        ImGui::TableSetupColumn("coeff (waves)");
                        ImGui::TableHeadersRow();
                        for (auto& z : g_res.wf.zernike) {
                            ImGui::TableNextRow();
                            ImGui::TableNextColumn(); ImGui::TextUnformatted(z.name.c_str());
                            ImGui::TableNextColumn();
                            // highlight the dominant aberrations
                            if (std::abs(z.coeff_waves) >= 0.02)
                                ImGui::TextColored(rgb(170, 90, 0), "%+.4f", z.coeff_waves);
                            else
                                ImGui::Text("%+.4f", z.coeff_waves);
                        }
                        ImGui::EndTable();
                    }
                    ImGui::TextDisabled("OPD map: blue = negative, red = positive (piston removed)");
                }
            }
            ImGui::End();
        }

        // ---- Chromatic ----
        if (win_chr) {
            if (ImGui::Begin("Chromatic", &win_chr)) {
                std::lock_guard<std::mutex> lk(g_mtx);
                if (g_res.chrom_focus.empty()) ImGui::TextDisabled("Run a design (F5).");
                else {
                    ImGui::TextUnformatted("Focal length (um) vs wavelength");
                    ImGui::PlotLines("##chrom", g_res.chrom_focus.data(),
                                     static_cast<int>(g_res.chrom_focus.size()), 0,
                                     nullptr, FLT_MAX, FLT_MAX, ImGui::GetContentRegionAvail());
                }
            }
            ImGui::End();
        }

        // ---- Through Focus ----
        if (win_tf) {
            if (ImGui::Begin("Through Focus", &win_tf)) {
                std::lock_guard<std::mutex> lk(g_mtx);
                if (g_res.tf.intensity.empty()) ImGui::TextDisabled("Run a design (F5).");
                else {
                    ImGui::Text("Peak at z = %.2f um   (nominal f = %.1f um)",
                                g_res.tf.z_peak_um, g_res.used.focal);
                    ImGui::Text("Depth of focus (axial FWHM): %.2f um", g_res.tf.dof_um);
                    char ov[64];
                    std::snprintf(ov, sizeof(ov), "on-axis intensity vs z (%.0f..%.0f um)",
                                  g_res.tf.z_um.front(), g_res.tf.z_um.back());
                    ImGui::PlotLines("##tf", g_res.tf.intensity.data(),
                                     static_cast<int>(g_res.tf.intensity.size()), 0,
                                     ov, 0.0f, 1.0f, ImVec2(-FLT_MIN, 160));
                }
            }
            ImGui::End();
        }

        // ---- MTF ----
        if (win_mtf) {
            if (ImGui::Begin("MTF", &win_mtf)) {
                std::lock_guard<std::mutex> lk(g_mtx);
                if (g_res.mtf.mtf.empty()) ImGui::TextDisabled("Run a design (F5).");
                else {
                    ImGui::Text("Diffraction cutoff: %.0f cycles/mm",
                                g_res.mtf.cutoff_cyc_mm);
                    ImGui::PlotLines("##mtf", g_res.mtf.mtf.data(),
                                     static_cast<int>(g_res.mtf.mtf.size()), 0,
                                     "lens MTF (vs spatial frequency)", 0.0f, 1.0f,
                                     ImVec2(-FLT_MIN, 120));
                    ImGui::PlotLines("##mtfi", g_res.mtf.mtf_ideal.data(),
                                     static_cast<int>(g_res.mtf.mtf_ideal.size()), 0,
                                     "diffraction-limited MTF", 0.0f, 1.0f,
                                     ImVec2(-FLT_MIN, 80));
                    if (ImGui::BeginTable("mtft", 3,
                                          ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg |
                                              ImGuiTableFlags_ScrollY, ImVec2(0, 140))) {
                        ImGui::TableSetupColumn("freq (cyc/mm)");
                        ImGui::TableSetupColumn("lens MTF");
                        ImGui::TableSetupColumn("ideal MTF");
                        ImGui::TableHeadersRow();
                        const auto& m = g_res.mtf;
                        int step = std::max<int>(1, static_cast<int>(m.freq_cyc_mm.size()) / 12);
                        for (std::size_t i = 0; i < m.freq_cyc_mm.size(); i += step) {
                            ImGui::TableNextRow();
                            ImGui::TableNextColumn(); ImGui::Text("%.0f", m.freq_cyc_mm[i]);
                            ImGui::TableNextColumn(); ImGui::Text("%.3f", m.mtf[i]);
                            ImGui::TableNextColumn(); ImGui::Text("%.3f", m.mtf_ideal[i]);
                        }
                        ImGui::EndTable();
                    }
                }
            }
            ImGui::End();
        }

        // ---- Tolerance ----
        if (win_tol) {
            if (ImGui::Begin("Tolerance", &win_tol)) {
                std::lock_guard<std::mutex> lk(g_mtx);
                if (!have_tol || g_tol.empty())
                    ImGui::TextDisabled("Run Tolerance (Run menu or Lens Data panel).");
                else if (ImGui::BeginTable("tolt", 4, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg)) {
                    ImGui::TableSetupColumn("sigma (nm)"); ImGui::TableSetupColumn("mean");
                    ImGui::TableSetupColumn("std"); ImGui::TableSetupColumn("worst");
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
            ImGui::End();
        }

        // ---- Field of View ----
        if (win_fov) {
            if (ImGui::Begin("Field of View", &win_fov)) {
                std::lock_guard<std::mutex> lk(g_mtx);
                if (!have_fov || g_fov.empty())
                    ImGui::TextDisabled("Run Field of View (Run menu or Lens Data panel).");
                else if (ImGui::BeginTable("fovt", 3, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg)) {
                    ImGui::TableSetupColumn("angle (deg)"); ImGui::TableSetupColumn("rel. Strehl");
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
            ImGui::End();
        }

        // ---- Log ----
        if (win_log) {
            if (ImGui::Begin("Log", &win_log)) {
                std::lock_guard<std::mutex> lk(g_mtx);
                ImGui::TextWrapped("%s%s", running ? "[BUSY]  " : "", g_status.c_str());
            }
            ImGui::End();
        }

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
