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
#include <sstream>
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
#include "celeris/design/polar_metalens.hpp"
#include "celeris/io/gds.hpp"
#include "celeris/io/image.hpp"
#include "celeris/io/material_csv.hpp"
#include "celeris/materials/database.hpp"
#ifdef CELERIS_USE_CUDA_KERNELS
#include "celeris/cuda/eigensolve.hpp"
#endif

using namespace celeris;

#include "app_state.hpp"
#include "theme.hpp"
#include "textures.hpp"

#ifdef _WIN32
#define NOMINMAX
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <psapi.h>  // process CPU/memory for the internal performance monitor
#define GLFW_EXPOSE_NATIVE_WIN32
#include <GLFW/glfw3native.h>
#include <dwmapi.h>  // immersive dark-mode title bar
#endif

using namespace celeris::gui;

#ifdef _WIN32
// Make the native OS title bar follow the app theme (Win10 2004+/Win11).
static void set_titlebar_dark(GLFWwindow* w, bool dark) {
    HWND hwnd = glfwGetWin32Window(w);
    BOOL v = dark ? TRUE : FALSE;
    // 20 = DWMWA_USE_IMMERSIVE_DARK_MODE (current); 19 on older builds.
    if (FAILED(DwmSetWindowAttribute(hwnd, 20, &v, sizeof(v))))
        DwmSetWindowAttribute(hwnd, 19, &v, sizeof(v));
}
#else
static void set_titlebar_dark(GLFWwindow*, bool) {}
#endif

// A procedural 48x48 RGBA app icon: concentric rings (the focusing/Airy motif)
// in the brand blue on transparent — no asset file, stays a single exe.
static void set_app_icon(GLFWwindow* window) {
    const int N = 48;
    static unsigned char px[N * N * 4];
    const double cx = (N - 1) / 2.0, cy = (N - 1) / 2.0, R = (N - 1) / 2.0;
    for (int y = 0; y < N; ++y)
        for (int x = 0; x < N; ++x) {
            double dx = x - cx, dy = y - cy, r = std::sqrt(dx * dx + dy * dy);
            unsigned char* p = &px[(y * N + x) * 4];
            if (r > R) { p[0] = p[1] = p[2] = p[3] = 0; continue; }
            double t = 0.5 + 0.5 * std::cos(r * 0.95);  // rings
            p[0] = static_cast<unsigned char>(10 + 70 * t);
            p[1] = static_cast<unsigned char>(40 + 130 * t);
            p[2] = static_cast<unsigned char>(90 + 165 * t);
            double edge = std::clamp(R - r, 0.0, 1.5) / 1.5;  // feather rim
            p[3] = static_cast<unsigned char>(255 * edge);
        }
    GLFWimage img{N, N, px};
    glfwSetWindowIcon(window, 1, &img);
}

int main() {
    if (!glfwInit()) return 1;
    GLFWwindow* window = glfwCreateWindow(1100, 720, "CELERIS — Metalens Designer",
                                          nullptr, nullptr);
    if (!window) { glfwTerminate(); return 1; }
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

    set_app_icon(window);
#ifdef CELERIS_GIT_COMMIT
    glfwSetWindowTitle(window,
        "CELERIS  —  Metalens Designer     |     build " CELERIS_GIT_COMMIT);
#endif

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
    bool dark_mode = false;
    apply_theme(dark_mode);
    set_titlebar_dark(window, dark_mode);
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 130");

    Params params;
    unsigned int psf_tex = 0, wf_tex = 0, layout_tex = 0, caustic_tex = 0;
    unsigned int polar_psf_x_tex = 0, polar_psf_y_tex = 0;
    int layout_mode = 0, layout_built_mode = -1;
    std::vector<unsigned int> spot_texs;
    bool have_result = false, have_tol = false, have_fov = false, have_spot = false;
    bool have_polar = false;
    bool gds_need_fit = false;  // refit the GDS viewer when a new design lands
    char gds_name[256] = "metalens.gds";
    double run_start = 0.0;  // ImGui time when the current run was launched

    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();

        if (g_pending.exchange(false)) {
            std::lock_guard<std::mutex> lk(g_mtx);
            upload_psf_texture(g_res.psf, psf_tex);
            upload_wavefront_texture(g_res.wf, wf_tex);
            upload_caustic_texture(g_res.tf, caustic_tex);
            layout_built_mode = -1;  // force layout texture rebuild for new design
            gds_need_fit = true;     // refit the GDS viewer to the new aperture
            have_result = true;
        }
        if (g_tol_pending.exchange(false)) have_tol = true;
        if (g_fov_pending.exchange(false)) have_fov = true;
        if (g_polar_pending.exchange(false)) {
            std::lock_guard<std::mutex> lk(g_mtx);
            upload_psf_texture(g_polar_psf_x, polar_psf_x_tex);
            upload_psf_texture(g_polar_psf_y, polar_psf_y_tex);
            have_polar = true;
        }
        if (g_spot_pending.exchange(false)) {
            std::lock_guard<std::mutex> lk(g_mtx);
            spot_texs.resize(g_spot.size(), 0);
            for (std::size_t i = 0; i < g_spot.size(); ++i)
                upload_psf_texture(g_spot[i].psf, spot_texs[i]);
            have_spot = true;
        }
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

        static bool show_about = false, show_save = false, show_open = false;
        static char proj_path[256] = "metalens.celeris";
        static std::string proj_msg;
        static bool win_lens = true, win_sum = true, win_foc = true, win_psf = true,
                    win_chr = true, win_tol = true, win_fov = true, win_log = true,
                    win_wf = true, win_mtf = true, win_tf = true, win_stack = true,
                    win_mats = true, win_layout = true, win_spot = true,
                    win_gds = true, win_polar = true, win_lib = true,
                    win_perf = false;
        const bool can_act = have_result && !running;

        if (ImGui::BeginMainMenuBar()) {
            if (ImGui::BeginMenu("File")) {
                if (ImGui::MenuItem("Save Project...")) show_save = true;
                if (ImGui::MenuItem("Open Project...")) show_open = true;
                ImGui::Separator();
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
                if (ImGui::MenuItem("Spot vs Field Diagram", nullptr, false, can_act))
                    std::thread(run_spotgrid).detach();
                ImGui::EndMenu();
            }
            if (ImGui::BeginMenu("View")) {
                ImGui::MenuItem("Lens Data", nullptr, &win_lens);
                ImGui::MenuItem("Layer Stack", nullptr, &win_stack);
                ImGui::MenuItem("Materials", nullptr, &win_mats);
                ImGui::MenuItem("Design Summary", nullptr, &win_sum);
                ImGui::MenuItem("Unit-Cell Library", nullptr, &win_lib);
                ImGui::MenuItem("Lens Layout", nullptr, &win_layout);
                ImGui::MenuItem("GDS Layout", nullptr, &win_gds);
                ImGui::MenuItem("Focus Performance", nullptr, &win_foc);
                ImGui::MenuItem("Focal PSF", nullptr, &win_psf);
                ImGui::MenuItem("Wavefront", nullptr, &win_wf);
                ImGui::MenuItem("MTF", nullptr, &win_mtf);
                ImGui::MenuItem("Through Focus", nullptr, &win_tf);
                ImGui::MenuItem("Spot vs Field", nullptr, &win_spot);
                ImGui::MenuItem("Polarization", nullptr, &win_polar);
                ImGui::MenuItem("Chromatic", nullptr, &win_chr);
                ImGui::MenuItem("Tolerance", nullptr, &win_tol);
                ImGui::MenuItem("Field of View", nullptr, &win_fov);
                ImGui::MenuItem("Log", nullptr, &win_log);
                ImGui::MenuItem("Performance", nullptr, &win_perf);
                ImGui::Separator();
                if (ImGui::MenuItem("Dark mode", nullptr, &dark_mode)) {
                    apply_theme(dark_mode);
                    set_titlebar_dark(window, dark_mode);
                }
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

        if (show_save) { ImGui::OpenPopup("Save Project"); show_save = false; proj_msg.clear(); }
        if (show_open) { ImGui::OpenPopup("Open Project"); show_open = false; proj_msg.clear(); }
        if (ImGui::BeginPopupModal("Save Project", nullptr,
                                   ImGuiWindowFlags_AlwaysAutoResize)) {
            ImGui::TextUnformatted("Write current lens parameters to a .celeris file:");
            ImGui::SetNextItemWidth(360);
            ImGui::InputText("##savepath", proj_path, sizeof(proj_path));
            if (!proj_msg.empty()) ImGui::TextColored(rgb(170, 90, 0), "%s", proj_msg.c_str());
            if (ImGui::Button("Save", ImVec2(120, 0))) {
                if (save_project(proj_path, params)) ImGui::CloseCurrentPopup();
                else proj_msg = "Could not write file.";
            }
            ImGui::SameLine();
            if (ImGui::Button("Cancel", ImVec2(120, 0))) ImGui::CloseCurrentPopup();
            ImGui::EndPopup();
        }
        if (ImGui::BeginPopupModal("Open Project", nullptr,
                                   ImGuiWindowFlags_AlwaysAutoResize)) {
            ImGui::TextUnformatted("Load lens parameters from a .celeris file:");
            ImGui::SetNextItemWidth(360);
            ImGui::InputText("##openpath", proj_path, sizeof(proj_path));
            if (!proj_msg.empty()) ImGui::TextColored(rgb(170, 90, 0), "%s", proj_msg.c_str());
            if (ImGui::Button("Open", ImVec2(120, 0))) {
                if (load_project(proj_path, params)) ImGui::CloseCurrentPopup();
                else proj_msg = "Not a valid .celeris file.";
            }
            ImGui::SameLine();
            if (ImGui::Button("Cancel", ImVec2(120, 0))) ImGui::CloseCurrentPopup();
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
            ImGui::DockBuilderDockWindow("Unit-Cell Library", ctl);
            ImGui::DockBuilderDockWindow("Focal PSF", ctr);
            ImGui::DockBuilderDockWindow("Lens Layout", ctr);
            ImGui::DockBuilderDockWindow("GDS Layout", ctr);
            ImGui::DockBuilderDockWindow("Polarization", ctr);
            ImGui::DockBuilderDockWindow("Wavefront", ctr);
            ImGui::DockBuilderDockWindow("Chromatic", ctr);
            ImGui::DockBuilderDockWindow("MTF", cbottom);
            ImGui::DockBuilderDockWindow("Through Focus", cbottom);
            ImGui::DockBuilderDockWindow("Tolerance", cbottom);
            ImGui::DockBuilderDockWindow("Field of View", cbottom);
            ImGui::DockBuilderDockWindow("Spot vs Field", cbottom);
            ImGui::DockBuilderDockWindow("Log", cbottom);
            ImGui::DockBuilderDockWindow("Performance", cbottom);
            ImGui::DockBuilderFinish(dockid);
        }

        auto kvrow = [](const char* k, const std::string& v) {
            ImGui::TableNextRow();
            ImGui::TableNextColumn(); ImGui::TextUnformatted(k);
            ImGui::TableNextColumn(); ImGui::TextUnformatted(v.c_str());
        };

        // Theme-aware colors for the custom-drawn canvases (GDS/polar viewers).
        const ImU32 canvas_bg = dark_mode ? IM_COL32(30, 30, 32, 255)
                                          : IM_COL32(248, 248, 248, 255);
        const ImU32 canvas_fg = dark_mode ? IM_COL32(225, 225, 225, 255)
                                          : IM_COL32(20, 20, 20, 255);
        const ImU32 pillar_col = dark_mode ? IM_COL32(90, 150, 235, 255)
                                           : IM_COL32(38, 78, 150, 255);

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
                static char rpt_prefix[256] = "celeris_report";
                ImGui::SetNextItemWidth(-FLT_MIN);
                ImGui::InputText("##rpt", rpt_prefix, sizeof(rpt_prefix));
                if (ImGui::Button("Save Report (txt + images + gds)", ImVec2(-FLT_MIN, 0))) {
                    Results r; Params up;
                    { std::lock_guard<std::mutex> lk(g_mtx); r = g_res; up = g_res.used; }
                    std::string base = rpt_prefix;
                    bool okall = true;
                    // Metrics report.
                    {
                        std::ofstream f(base + "_report.txt");
                        if (f) {
                            f << "CELERIS metalens design report\n";
                            f << "==============================\n\n";
                            f << std::format("focal length      : {} um\n", up.focal);
                            f << std::format("aperture diameter : {} um\n", up.diameter);
                            f << std::format("wavelength        : {} um\n", up.wavelength);
                            f << std::format("period            : {} um\n", up.period);
                            f << std::format("pillar height     : {} um\n", up.thickness);
                            f << std::format("pillar material   : {}\n", kPillarMats[up.pillar_mat]);
                            f << std::format("substrate         : {}\n\n", kSubstrates[up.substrate_mat]);
                            f << std::format("array             : {} x {} ({} pillars)\n", r.n_cells, r.n_cells, r.pillars);
                            f << std::format("numerical aperture: {:.3f}\n", r.na);
                            f << std::format("phase coverage    : {:.0f} deg\n", r.coverage_deg);
                            f << std::format("RMS phase error   : {:.1f} deg\n", r.rms);
                            f << std::format("mean transmission : {:.3f}\n\n", r.meanT);
                            f << std::format("Strehl ratio      : {:.3f}\n", r.strehl);
                            f << std::format("spot FWHM         : {:.3f} um\n", r.fwhm);
                            f << std::format("diffraction limit : {:.3f} um\n", r.dl);
                            f << std::format("encircled energy  : {:.1f} %\n", r.encircled * 100.0);
                            f << std::format("wavefront RMS     : {:.4f} waves\n", r.wf.rms_waves);
                            f << std::format("depth of focus    : {:.2f} um\n", r.tf.dof_um);
                        } else okall = false;
                    }
                    // PSF image (gamma-boosted for sidelobes).
                    if (!r.psf.intensity.empty())
                        okall &= write_pgm(base + "_psf.pgm", r.psf.n, r.psf.n,
                                           r.psf.intensity, 2.2);
                    // Caustic image (convert normalized floats to doubles).
                    if (!r.tf.caustic.empty()) {
                        std::vector<double> cd(r.tf.caustic.begin(), r.tf.caustic.end());
                        okall &= write_pgm(base + "_caustic.pgm", r.tf.caustic_nx,
                                           r.tf.caustic_nz, cd, 1.6);
                    }
                    okall &= (write_metalens_gds(r.design, base + "_layout.gds") >= 0);
                    std::lock_guard<std::mutex> lk(g_mtx);
                    g_status = okall ? std::format("Report written: {}_report.txt (+ psf/caustic/gds)", base)
                                     : std::string("ERROR: report write failed (check path)");
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

        // ---- Unit-Cell Library (meta-atom design space) ----
        if (win_lib) {
            if (ImGui::Begin("Unit-Cell Library", &win_lib)) {
                if (!have_result) ImGui::TextDisabled("Run a design (F5).");
                else {
                    std::lock_guard<std::mutex> lk(g_mtx);
                    const auto& lib = g_res.lib;
                    int ns = static_cast<int>(lib.fill.size());
                    if (ns < 2) ImGui::TextDisabled("No library data.");
                    else {
                        std::vector<float> ph(ns), am(ns);
                        for (int i = 0; i < ns; ++i) {
                            ph[i] = static_cast<float>(lib.phase[i] * 180.0 / 3.14159265358979);
                            am[i] = static_cast<float>(lib.amplitude[i]);
                        }
                        ImGui::Text("Phase coverage: %.0f deg  (need ~360 for full 2pi control)",
                                    g_res.coverage_deg);
                        if (g_res.coverage_deg < 300.0)
                            ImGui::TextColored(rgb(190, 40, 40),
                                "limited coverage -> taller pillars or higher index for full 2pi");
                        ImGui::PlotLines("##phfill", ph.data(), ns, 0,
                                         "imparted phase (deg) vs fill", -180.0f, 180.0f,
                                         ImVec2(-FLT_MIN, 110));
                        ImGui::PlotLines("##amfill", am.data(), ns, 0,
                                         "transmission |t| vs fill", 0.0f, 1.0f,
                                         ImVec2(-FLT_MIN, 90));
                        ImGui::Text("fill range %.2f..%.2f over %d samples, period %.3f um",
                                    lib.fill.front(), lib.fill.back(), ns, lib.period_um);
                        ImGui::TextDisabled("Each point is one RCWA solve. Phase must span "
                                            "2pi for an ideal lens; |t| sets efficiency.");
                    }
                }
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

        // ---- Polarization (polarization-multiplexed design) ----
        if (win_polar) {
            if (ImGui::Begin("Polarization", &win_polar)) {
                ImGui::TextWrapped("Polarization-multiplexed lens: X-pol and Y-pol "
                                   "focus at independent distances (rectangular "
                                   "pillars).");
                ImGui::SetNextItemWidth(120);
                ImGui::InputFloat("Focal X-pol (um)", &params.focal, 0, 0, "%.1f");
                ImGui::SetNextItemWidth(120);
                ImGui::InputFloat("Focal Y-pol (um)", &params.focal_y, 0, 0, "%.1f");
                bool busy = g_running.load();
                if (busy) ImGui::BeginDisabled();
                if (ImGui::Button("Run Polarization Design", ImVec2(-FLT_MIN, 28)))
                    std::thread(run_polardesign, params).detach();
                if (busy) ImGui::EndDisabled();

                if (!have_polar)
                    ImGui::TextDisabled("Set focal lengths and run.");
                else {
                    std::lock_guard<std::mutex> lk(g_mtx);
                    ImGui::Separator();
                    ImGui::Text("X-pol  @ %.0f um:  RMS %.1f deg   mean |t| %.3f",
                                g_polar_fx, g_polar.rms_phase_error_x_deg, g_polar.mean_amp_x);
                    ImGui::Text("Y-pol  @ %.0f um:  RMS %.1f deg   mean |t| %.3f",
                                g_polar_fy, g_polar.rms_phase_error_y_deg, g_polar.mean_amp_y);
                    ImGui::Text("measured focus:  X %.1f um   Y %.1f um", g_polar_zx, g_polar_zy);
                    if (g_polar_iso_x != 0.0 || g_polar_iso_y != 0.0)
                        ImGui::Text("focal isolation: X-plane %.1f dB   Y-plane %.1f dB",
                                    g_polar_iso_x, g_polar_iso_y);
                    // Focal-plane PSF of each polarization at its target plane:
                    // two tight spots = the bifocal split, proven optically.
                    if (polar_psf_x_tex && polar_psf_y_tex) {
                        float s = 150.0f;
                        ImGui::BeginGroup();
                        ImGui::TextUnformatted("X-pol focus");
                        ImGui::Image((ImTextureID)(intptr_t)polar_psf_x_tex, ImVec2(s, s));
                        ImGui::EndGroup();
                        ImGui::SameLine();
                        ImGui::BeginGroup();
                        ImGui::TextUnformatted("Y-pol focus");
                        ImGui::Image((ImTextureID)(intptr_t)polar_psf_y_tex, ImVec2(s, s));
                        ImGui::EndGroup();
                        ImGui::TextDisabled("each polarization imaged at its own focal plane");
                    }
                    static char ppath[256] = "polar_metalens.gds";
                    ImGui::SetNextItemWidth(-110);
                    ImGui::InputText("##ppath", ppath, sizeof(ppath)); ImGui::SameLine();
                    if (ImGui::Button("Save GDS", ImVec2(-FLT_MIN, 0))) {
                        int np = write_rect_gds(ppath, g_polar.n_cells, g_polar.period_um,
                                                g_polar.fill_x, g_polar.fill_y);
                        g_status = np >= 0 ? std::format("Wrote {} rect pillars -> {}", np, ppath)
                                           : std::string("ERROR: GDS write failed");
                    }
                    // Rectangular-pillar layout canvas (pan/zoom, culled).
                    static float ps = 0.0f; static ImVec2 pcam(0, 0); static bool pfit = true;
                    if (ImGui::Button("Fit")) pfit = true;
                    const int n = g_polar.n_cells; const double pp = g_polar.period_um;
                    const double cen = (n - 1) / 2.0, ext = n * pp;
                    ImVec2 av = ImGui::GetContentRegionAvail();
                    av.x = std::max(av.x, 60.0f); av.y = std::max(av.y, 60.0f);
                    ImVec2 q0 = ImGui::GetCursorScreenPos();
                    ImVec2 q1 = ImVec2(q0.x + av.x, q0.y + av.y);
                    ImGui::InvisibleButton("polcanvas", av);
                    bool hov = ImGui::IsItemHovered();
                    ImVec2 cc = ImVec2((q0.x + q1.x) * 0.5f, (q0.y + q1.y) * 0.5f);
                    if (pfit || ps <= 0) { ps = std::min(av.x, av.y) / float(ext * 1.1); pcam = ImVec2(0, 0); pfit = false; }
                    auto w2s = [&](double wx, double wy) {
                        return ImVec2(cc.x + float((wx - pcam.x) * ps), cc.y - float((wy - pcam.y) * ps));
                    };
                    auto s2w = [&](ImVec2 s) { return ImVec2(pcam.x + (s.x - cc.x) / ps, pcam.y - (s.y - cc.y) / ps); };
                    if (hov && io.MouseWheel != 0) {
                        ImVec2 wb = s2w(io.MousePos); ps *= std::pow(1.15f, io.MouseWheel);
                        ImVec2 wa = s2w(io.MousePos); pcam.x += wb.x - wa.x; pcam.y += wb.y - wa.y;
                    }
                    if (ImGui::IsItemActive() && ImGui::IsMouseDragging(ImGuiMouseButton_Left)) {
                        pcam.x -= io.MouseDelta.x / ps; pcam.y += io.MouseDelta.y / ps;
                    }
                    ImDrawList* dl = ImGui::GetWindowDrawList();
                    dl->PushClipRect(q0, q1, true);
                    dl->AddRectFilled(q0, q1, canvas_bg);
                    ImVec2 wtl = s2w(q0), wbr = s2w(q1);
                    double xmn = std::min(wtl.x, wbr.x), xmx = std::max(wtl.x, wbr.x);
                    double ymn = std::min(wtl.y, wbr.y), ymx = std::max(wtl.y, wbr.y);
                    int ixmin = std::max(0, (int)std::floor(xmn / pp + cen));
                    int ixmax = std::min(n - 1, (int)std::ceil(xmx / pp + cen));
                    int iymin = std::max(0, (int)std::floor(ymn / pp + cen));
                    int iymax = std::min(n - 1, (int)std::ceil(ymx / pp + cen));
                    long vis = (long)std::max(0, ixmax - ixmin + 1) * std::max(0, iymax - iymin + 1);
                    int stride = vis > 60000 ? (int)std::ceil(std::sqrt((double)vis / 60000)) : 1;
                    const ImU32 col = pillar_col;
                    for (int iy = iymin; iy <= iymax; iy += stride)
                        for (int ix = ixmin; ix <= ixmax; ix += stride) {
                            std::size_t off = (std::size_t)iy * n + ix;
                            double fxw = g_polar.fill_x[off], fyw = g_polar.fill_y[off];
                            if (fxw < 0.05 && fyw < 0.05) continue;
                            double cx = (ix - cen) * pp, cy = (iy - cen) * pp;
                            ImVec2 a = w2s(cx - 0.5 * fxw * pp, cy + 0.5 * fyw * pp);
                            ImVec2 b = w2s(cx + 0.5 * fxw * pp, cy - 0.5 * fyw * pp);
                            if (b.x - a.x < 1.0f) dl->AddRectFilled(a, ImVec2(a.x + 1, a.y + 1), col);
                            else dl->AddRectFilled(a, b, col);
                        }
                    dl->PopClipRect();
                    ImGui::SetCursorScreenPos(ImVec2(q0.x + 6, q0.y + 6));
                    ImGui::Text("%d x %d rectangular pillars   aperture %.1f um", n, n, ext);
                }
            }
            ImGui::End();
        }

        // ---- GDS Layout (in-app fab-polygon viewer, pan/zoom) ----
        if (win_gds) {
            if (ImGui::Begin("GDS Layout", &win_gds)) {
                static float gscale = 0.0f;      // pixels per micron
                static ImVec2 gcam(0.0f, 0.0f);  // world point at canvas center (um)
                static bool gphase = false;
                static GdsLayout gloaded;        // a .gds opened from disk
                static char gpath[256] = "metalens.gds";
                static bool show_file = false;   // false = live design, true = file
                static std::string gmsg;

                // Toolbar.
                if (ImGui::Button("Fit")) gds_need_fit = true;
                ImGui::SameLine();
                if (!show_file) {
                    ImGui::Checkbox("Color by phase", &gphase); ImGui::SameLine();
                }
                ImGui::SetNextItemWidth(220);
                ImGui::InputText("##gdspath", gpath, sizeof(gpath)); ImGui::SameLine();
                if (ImGui::Button("Open GDS")) {
                    GdsLayout L = read_gds(gpath);
                    if (L.ok) { gloaded = std::move(L); show_file = true;
                                gds_need_fit = true;
                                gmsg = std::format("Loaded {} polygons", gloaded.polygons.size()); }
                    else gmsg = "Could not parse GDS file.";
                }
                if (gloaded.ok) {
                    ImGui::SameLine();
                    if (ImGui::RadioButton("Design", !show_file)) { show_file = false; gds_need_fit = true; }
                    ImGui::SameLine();
                    if (ImGui::RadioButton("File", show_file)) { show_file = true; gds_need_fit = true; }
                }
                if (!gmsg.empty()) { ImGui::SameLine(); ImGui::TextDisabled("%s", gmsg.c_str()); }

                bool drawing_file = show_file && gloaded.ok;
                if (!have_result && !drawing_file) {
                    ImGui::TextDisabled("Run a design (F5), or type a path and Open GDS.");
                } else {
                    std::lock_guard<std::mutex> lk(g_mtx);
                    const MetalensDesign& d = g_res.design;
                    const UnitCellLibrary& glib = g_res.lib;
                    const int n = d.n_cells;
                    const double pp = d.period_um;
                    const double cen = (n - 1) / 2.0;

                    // World bounding box of whatever we're showing.
                    double bxmin, bxmax, bymin, bymax;
                    if (drawing_file) {
                        bxmin = gloaded.min_x; bxmax = gloaded.max_x;
                        bymin = gloaded.min_y; bymax = gloaded.max_y;
                    } else {
                        double ext = n * pp;
                        bxmin = -ext / 2; bxmax = ext / 2;
                        bymin = -ext / 2; bymax = ext / 2;
                    }
                    double bcx = (bxmin + bxmax) / 2, bcy = (bymin + bymax) / 2;
                    double bsize = std::max(std::max(bxmax - bxmin, bymax - bymin), 1e-6);

                    ImVec2 avail = ImGui::GetContentRegionAvail();
                    avail.x = std::max(avail.x, 60.0f);
                    avail.y = std::max(avail.y, 60.0f);
                    ImVec2 c0 = ImGui::GetCursorScreenPos();
                    ImVec2 c1 = ImVec2(c0.x + avail.x, c0.y + avail.y);
                    ImGui::InvisibleButton("gdscanvas", avail);
                    bool hov = ImGui::IsItemHovered();
                    ImVec2 cc = ImVec2((c0.x + c1.x) * 0.5f, (c0.y + c1.y) * 0.5f);

                    if (gds_need_fit || gscale <= 0.0f) {
                        gscale = std::min(avail.x, avail.y) /
                                 static_cast<float>(bsize * 1.1);
                        gcam = ImVec2(static_cast<float>(bcx), static_cast<float>(bcy));
                        gds_need_fit = false;
                    }

                    auto w2s = [&](double wx, double wy) {
                        return ImVec2(cc.x + static_cast<float>((wx - gcam.x) * gscale),
                                      cc.y - static_cast<float>((wy - gcam.y) * gscale));
                    };
                    auto s2w = [&](ImVec2 s) {
                        return ImVec2(gcam.x + (s.x - cc.x) / gscale,
                                      gcam.y - (s.y - cc.y) / gscale);
                    };

                    if (hov && io.MouseWheel != 0.0f) {
                        ImVec2 wb = s2w(io.MousePos);
                        gscale *= std::pow(1.15f, io.MouseWheel);
                        ImVec2 wa = s2w(io.MousePos);
                        gcam.x += wb.x - wa.x;
                        gcam.y += wb.y - wa.y;
                    }
                    if (ImGui::IsItemActive() &&
                        ImGui::IsMouseDragging(ImGuiMouseButton_Left)) {
                        gcam.x -= io.MouseDelta.x / gscale;
                        gcam.y += io.MouseDelta.y / gscale;
                    }

                    ImDrawList* dl = ImGui::GetWindowDrawList();
                    dl->PushClipRect(c0, c1, true);
                    dl->AddRectFilled(c0, c1, canvas_bg);

                    ImVec2 wtl = s2w(c0), wbr = s2w(c1);
                    double wxmin = std::min(wtl.x, wbr.x), wxmax = std::max(wtl.x, wbr.x);
                    double wymin = std::min(wtl.y, wbr.y), wymax = std::max(wtl.y, wbr.y);
                    const ImU32 mask_col = pillar_col;
                    int stride = 1;

                    if (drawing_file) {
                        // Render parsed polygons (viewport-culled, strided if huge).
                        long np = static_cast<long>(gloaded.polygons.size());
                        const long cap = 80000;
                        if (np > cap) stride = (int)std::ceil((double)np / cap);
                        for (long i = 0; i < np; i += stride) {
                            const auto& poly = gloaded.polygons[i];
                            if (poly.pts.size() < 3) continue;
                            // Quick bbox cull.
                            double pminx = 1e300, pmaxx = -1e300, pminy = 1e300, pmaxy = -1e300;
                            for (auto& q : poly.pts) {
                                pminx = std::min(pminx, q.first); pmaxx = std::max(pmaxx, q.first);
                                pminy = std::min(pminy, q.second); pmaxy = std::max(pmaxy, q.second);
                            }
                            if (pmaxx < wxmin || pminx > wxmax || pmaxy < wymin || pminy > wymax)
                                continue;
                            // Drop the duplicate closing vertex if present.
                            std::size_t m = poly.pts.size();
                            if (m > 1 && poly.pts.front() == poly.pts.back()) --m;
                            static std::vector<ImVec2> sp; sp.clear();
                            for (std::size_t k = 0; k < m; ++k)
                                sp.push_back(w2s(poly.pts[k].first, poly.pts[k].second));
                            dl->AddConvexPolyFilled(sp.data(), (int)sp.size(), mask_col);
                        }
                    } else {
                        // Render the live design's pillar squares (matches the .gds).
                        int ixmin = std::max(0, (int)std::floor(wxmin / pp + cen));
                        int ixmax = std::min(n - 1, (int)std::ceil(wxmax / pp + cen));
                        int iymin = std::max(0, (int)std::floor(wymin / pp + cen));
                        int iymax = std::min(n - 1, (int)std::ceil(wymax / pp + cen));
                        long vis = (long)std::max(0, ixmax - ixmin + 1) *
                                   std::max(0, iymax - iymin + 1);
                        const long cap = 60000;
                        if (vis > cap) stride = (int)std::ceil(std::sqrt((double)vis / cap));
                        for (int iy = iymin; iy <= iymax; iy += stride)
                            for (int ix = ixmin; ix <= ixmax; ix += stride) {
                                double fill = d.fill_map[(std::size_t)iy * n + ix];
                                if (fill < 0.05) continue;  // matches GDS min_fill
                                double cx = (ix - cen) * pp, cy = (iy - cen) * pp;
                                double half = 0.5 * fill * pp;
                                ImVec2 a = w2s(cx - half, cy + half);
                                ImVec2 b = w2s(cx + half, cy - half);
                                ImU32 col = mask_col;
                                if (gphase) {
                                    cdouble t = glib.transmission_for_fill(fill);
                                    double h = (std::arg(t) + pi) / (2.0 * pi);
                                    std::uint8_t r, g, bl;
                                    hue_to_rgb(std::clamp(h, 0.0, 1.0), r, g, bl);
                                    col = IM_COL32(r, g, bl, 255);
                                }
                                if (b.x - a.x < 1.0f)
                                    dl->AddRectFilled(a, ImVec2(a.x + 1.0f, a.y + 1.0f), col);
                                else
                                    dl->AddRectFilled(a, b, col);
                            }
                    }
                    dl->PopClipRect();

                    // Overlay readout.
                    ImGui::SetCursorScreenPos(ImVec2(c0.x + 6, c0.y + 6));
                    if (drawing_file)
                        ImGui::Text("FILE: %zu polygons   extent %.1f x %.1f um",
                                    gloaded.polygons.size(), bxmax - bxmin, bymax - bymin);
                    else
                        ImGui::Text("DESIGN: %d x %d cells   period %.3f um   aperture %.1f um",
                                    n, n, pp, n * pp);
                    if (stride > 1) {
                        ImGui::SetCursorScreenPos(ImVec2(c0.x + 6, c0.y + 26));
                        ImGui::TextColored(rgb(170, 90, 0),
                                           "thinned 1:%d (zoom in for full detail)", stride);
                    }
                    // Scale bar: round micron length ~1/5 of the canvas width.
                    double bar_um = std::pow(10.0, std::floor(std::log10(
                                        (avail.x / gscale) / 5.0)));
                    float bar_px = static_cast<float>(bar_um * gscale);
                    ImVec2 bs(c1.x - bar_px - 12, c1.y - 16);
                    dl->AddLine(bs, ImVec2(bs.x + bar_px, bs.y), canvas_fg, 2.0f);
                    char blab[32]; std::snprintf(blab, sizeof(blab), "%g um", bar_um);
                    dl->AddText(ImVec2(bs.x, bs.y - 16), canvas_fg, blab);
                }
            }
            ImGui::End();
        }

        // ---- Lens Layout (pillar map) ----
        if (win_layout) {
            if (ImGui::Begin("Lens Layout", &win_layout)) {
                if (!have_result) ImGui::TextDisabled("Run a design (F5).");
                else {
                    ImGui::TextUnformatted("View:"); ImGui::SameLine();
                    ImGui::RadioButton("Phase", &layout_mode, 0); ImGui::SameLine();
                    ImGui::RadioButton("Fill fraction", &layout_mode, 1);
                    if (layout_mode != layout_built_mode) {
                        std::lock_guard<std::mutex> lk(g_mtx);
                        upload_layout_texture(g_res.design, g_res.lib, layout_mode,
                                              layout_tex);
                        layout_built_mode = layout_mode;
                    }
                    int nc; double per;
                    { std::lock_guard<std::mutex> lk(g_mtx);
                      nc = g_res.design.n_cells; per = g_res.design.period_um; }
                    ImGui::Text("%d x %d pillars   period %.3f um   aperture %.1f um",
                                nc, nc, per, nc * per);
                    if (layout_tex) {
                        ImVec2 a = ImGui::GetContentRegionAvail();
                        float s = std::max(64.0f, std::min(a.x, a.y));
                        ImGui::Image((ImTextureID)(intptr_t)layout_tex, ImVec2(s, s));
                    }
                    if (layout_mode == 0)
                        ImGui::TextDisabled("Cyclic colormap: hue = imparted phase "
                                            "(-pi..pi). Rings are Fresnel zones.");
                    else
                        ImGui::TextDisabled("Grayscale: pillar fill fraction "
                                            "(black = 0, white = 1).");
                }
            }
            ImGui::End();
        }

        // ---- Spot vs Field (off-axis spot diagram) ----
        if (win_spot) {
            if (ImGui::Begin("Spot vs Field", &win_spot)) {
                bool can = have_result && !g_running.load();
                if (ImGui::Button("Compute spot diagram") && can)
                    std::thread(run_spotgrid).detach();
                ImGui::SameLine();
                ImGui::TextDisabled("off-axis focal spots vs field angle");
                if (!have_spot)
                    ImGui::TextDisabled("Run a design (F5), then Compute.");
                else {
                    std::lock_guard<std::mutex> lk(g_mtx);
                    const float sz = 130.0f;
                    for (std::size_t i = 0; i < g_spot.size() &&
                                           i < spot_texs.size(); ++i) {
                        ImGui::BeginGroup();
                        ImGui::Text("%.0f deg", g_spot[i].angle_deg);
                        if (spot_texs[i])
                            ImGui::Image((ImTextureID)(intptr_t)spot_texs[i],
                                         ImVec2(sz, sz));
                        ImGui::Text("Strehl %.2f", g_spot[i].rel_strehl);
                        ImGui::Text("shift %.1f um", g_spot[i].cx_um);
                        ImGui::EndGroup();
                        if (i + 1 < g_spot.size()) ImGui::SameLine();
                    }
                    ImGui::TextDisabled("Each tile is centered on the chief-ray "
                                        "landing point; coma grows with angle.");
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
                                     ov, 0.0f, 1.0f, ImVec2(-FLT_MIN, 110));
                    if (caustic_tex && g_res.tf.caustic_nx > 0) {
                        ImGui::SeparatorText("Caustic (x-z focal slice)");
                        ImVec2 a = ImGui::GetContentRegionAvail();
                        float ih = std::min(std::max(a.y, 80.0f), 240.0f);
                        // Draw with z horizontal (axial) x vertical: image is
                        // nx wide (lateral) x nz tall (axial); show as-is.
                        ImGui::Image((ImTextureID)(intptr_t)caustic_tex,
                                     ImVec2(a.x, ih));
                        ImGui::TextDisabled("horizontal = lateral x (%.1f..%.1f um), "
                                            "vertical = axial z (%.0f..%.0f um)",
                                            g_res.tf.caustic_xmin, g_res.tf.caustic_xmax,
                                            g_res.tf.caustic_zmin, g_res.tf.caustic_zmax);
                    }
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

        // ---- Performance (internal "task manager") ----
        if (win_perf) {
            if (ImGui::Begin("Performance", &win_perf)) {
                // Frame-rate + a rolling frame-time history.
                static float ft_hist[120] = {0};
                static int ft_ofs = 0;
                float frame_ms = io.DeltaTime * 1000.0f;
                ft_hist[ft_ofs] = frame_ms;
                ft_ofs = (ft_ofs + 1) % IM_ARRAYSIZE(ft_hist);
                ImGui::Text("FPS: %.1f   (%.2f ms/frame)", io.Framerate, 1000.0f / io.Framerate);
                char ov[32]; std::snprintf(ov, sizeof(ov), "%.2f ms", frame_ms);
                ImGui::PlotLines("##ft", ft_hist, IM_ARRAYSIZE(ft_hist), ft_ofs, ov,
                                 0.0f, 40.0f, ImVec2(-FLT_MIN, 60));

                ImGui::SeparatorText("Process");
#ifdef _WIN32
                // Process CPU% (kernel+user time delta / wall delta / cores) and
                // working-set memory, refreshed a few times a second.
                static double cpu_pct = 0.0, mem_mb = 0.0;
                static double last_wall = 0.0;
                static unsigned long long last_kt = 0, last_ut = 0;
                static unsigned ncpu = std::max(1u, std::thread::hardware_concurrency());
                double now = ImGui::GetTime();
                FILETIME fc, fe, fk, fu;
                if (GetProcessTimes(GetCurrentProcess(), &fc, &fe, &fk, &fu)) {
                    auto u64 = [](FILETIME f) {
                        return (static_cast<unsigned long long>(f.dwHighDateTime) << 32) |
                               f.dwLowDateTime;
                    };
                    unsigned long long kt = u64(fk), ut = u64(fu);
                    if (last_wall > 0.0) {
                        double dw = now - last_wall;
                        if (dw >= 0.25) {  // 100ns ticks -> seconds: *1e-7
                            double dproc = static_cast<double>((kt - last_kt) + (ut - last_ut)) * 1e-7;
                            cpu_pct = 100.0 * dproc / (dw * ncpu);
                            last_wall = now; last_kt = kt; last_ut = ut;
                        }
                    } else { last_wall = now; last_kt = kt; last_ut = ut; }
                }
                PROCESS_MEMORY_COUNTERS pmc;
                if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc)))
                    mem_mb = pmc.WorkingSetSize / 1048576.0;
                ImGui::Text("CPU: %.0f%%  of %u cores", cpu_pct, ncpu);
                ImGui::Text("Memory (working set): %.0f MB", mem_mb);
#else
                ImGui::TextDisabled("CPU/memory stats: Windows only");
#endif
                ImGui::SeparatorText("Compute");
                ImGui::Text("Workers: %s", g_running.load() ? "BUSY" : "idle");
#ifdef CELERIS_USE_CUDA_KERNELS
                if (cuda::available())
                    ImGui::TextColored(rgb(20, 130, 40), "Device: GPU (%s)", cuda::device_name());
                else
                    ImGui::TextColored(rgb(170, 90, 0), "Device: CPU (no CUDA device)");
#else
                ImGui::TextDisabled("Device: CPU (CPU build)");
#endif
                ImGui::Text("ImGui: %d verts, %d draw cmds",
                            ImGui::GetIO().MetricsRenderVertices,
                            ImGui::GetIO().MetricsRenderIndices / 3);
            }
            ImGui::End();
        }

        // ---- Log ----
        if (win_log) {
            if (ImGui::Begin("Log", &win_log)) {
                std::lock_guard<std::mutex> lk(g_mtx);
                ImGui::TextWrapped("%s%s", running ? "[BUSY]  " : "", g_status.c_str());
                ImGui::Separator();
#ifdef CELERIS_USE_CUDA_KERNELS
                if (cuda::available())
                    ImGui::TextColored(rgb(20, 130, 40), "Propagation: GPU  (%s)",
                                       cuda::device_name());
                else
                    ImGui::TextColored(rgb(170, 90, 0),
                                       "Propagation: CPU  (no CUDA device found)");
#else
                ImGui::TextDisabled("Propagation: CPU (multi-core; build with CUDA "
                                    "for GPU acceleration)");
#endif
            }
            ImGui::End();
        }

        ImGui::Render();
        int w, h;
        glfwGetFramebufferSize(window, &w, &h);
        glViewport(0, 0, w, h);
        if (dark_mode) glClearColor(0.10f, 0.10f, 0.12f, 1.0f);
        else glClearColor(0.55f, 0.55f, 0.56f, 1.0f);
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
