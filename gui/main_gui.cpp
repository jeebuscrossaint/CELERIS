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
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#include "celeris/analysis/chromatic.hpp"
#include "celeris/analysis/focal.hpp"
#include "celeris/design/metalens.hpp"
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
    double coverage_deg = 0;
    int n_cells = 0, pillars = 0;
    PsfMap psf;
    std::vector<float> chrom_wl, chrom_focus;
};

std::atomic<bool> g_running{false};
std::atomic<bool> g_pending{false};  // worker finished, main thread must ingest
std::mutex g_mtx;
std::string g_status = "Ready — set parameters and click Design.";
Results g_res;

void run_design(Params p) {
    g_running = true;
    {
        std::lock_guard<std::mutex> lk(g_mtx);
        g_status = "Building unit-cell library (RCWA sweep)...";
    }
    const auto pillar = Material::constant(cdouble{p.pillar_n, 0.0}, "pillar");
    auto lib = build_unit_cell_library(pillar, materials::air(), materials::air(),
                                       materials::bk7(), p.period, p.wavelength,
                                       p.thickness, 0.08, 0.92, p.fill_samples,
                                       p.harmonics);
    auto lens = design_metalens(lib, p.focal, p.diameter);
    auto foc = analyze_focus(lens, lib, p.focal, p.wavelength, p.diameter);
    double dl = p.wavelength * p.focal / p.diameter;
    auto psf = compute_psf(lens, lib, p.focal, p.wavelength, p.diameter, 161,
                           std::max(5.0 * dl, 4.0));
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
    r.n_cells = lens.n_cells;
    r.pillars = lens.n_cells * lens.n_cells;
    r.psf = std::move(psf);
    for (auto& c : chrom) {
        r.chrom_wl.push_back(static_cast<float>(c.wavelength_um * 1000.0));
        r.chrom_focus.push_back(static_cast<float>(c.focal_length_um));
    }
    {
        std::lock_guard<std::mutex> lk(g_mtx);
        g_res = std::move(r);
        g_status = "Done.";
    }
    g_pending = true;
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
    ImGui::StyleColorsDark();
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 130");

    Params params;
    unsigned int psf_tex = 0;
    bool have_result = false;

    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();

        if (g_pending.exchange(false)) {
            std::lock_guard<std::mutex> lk(g_mtx);
            upload_psf_texture(g_res.psf, psf_tex);
            have_result = true;
        }

        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        ImGui::SetNextWindowPos({0, 0});
        ImGui::SetNextWindowSize(ImGui::GetIO().DisplaySize);
        ImGui::Begin("CELERIS", nullptr,
                     ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove |
                         ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoTitleBar);

        ImGui::Text("CELERIS — GPU-ready metalens design (RCWA)");
        ImGui::Separator();

        ImGui::BeginChild("controls", ImVec2(330, 0), true);
        ImGui::TextUnformatted("Lens specification");
        ImGui::Spacing();
        ImGui::SliderFloat("focal (um)", &params.focal, 5.0f, 200.0f, "%.1f");
        ImGui::SliderFloat("diameter (um)", &params.diameter, 4.0f, 60.0f, "%.1f");
        ImGui::SliderFloat("wavelength (um)", &params.wavelength, 0.4f, 0.7f, "%.3f");
        ImGui::SliderFloat("period (um)", &params.period, 0.15f, 0.5f, "%.3f");
        ImGui::SliderFloat("pillar height (um)", &params.thickness, 0.2f, 1.2f, "%.2f");
        ImGui::SliderFloat("pillar index", &params.pillar_n, 1.5f, 3.5f, "%.2f");
        ImGui::SliderInt("RCWA harmonics", &params.harmonics, 4, 10);
        ImGui::SliderInt("library samples", &params.fill_samples, 8, 30);
        ImGui::Spacing();

        bool running = g_running.load();
        if (running) ImGui::BeginDisabled();
        if (ImGui::Button("Design lens", ImVec2(-1, 36))) {
            std::thread(run_design, params).detach();
        }
        if (running) ImGui::EndDisabled();

        ImGui::Spacing();
        {
            std::lock_guard<std::mutex> lk(g_mtx);
            ImGui::TextWrapped("%s", g_status.c_str());
        }
        ImGui::EndChild();

        ImGui::SameLine();
        ImGui::BeginChild("results", ImVec2(0, 0), true);
        if (!have_result) {
            ImGui::TextDisabled("Results will appear here after a design run.");
        } else {
            std::lock_guard<std::mutex> lk(g_mtx);
            ImGui::Text("Design: %dx%d pillars (%d total)", g_res.n_cells,
                        g_res.n_cells, g_res.pillars);
            ImGui::Text("Library phase coverage: %.0f deg", g_res.coverage_deg);
            ImGui::Text("RMS phase error: %.1f deg    mean |t|: %.3f", g_res.rms,
                        g_res.meanT);
            ImGui::Separator();
            ImGui::Text("Strehl ratio:    %.3f", g_res.strehl);
            ImGui::Text("Spot FWHM:       %.2f um  (diffraction limit %.2f um)",
                        g_res.fwhm, g_res.dl);
            ImGui::Text("Encircled energy: %.0f%%", g_res.encircled * 100.0);
            ImGui::Separator();
            ImGui::TextUnformatted("Focal point-spread function:");
            if (psf_tex)
                ImGui::Image((ImTextureID)(intptr_t)psf_tex, ImVec2(300, 300));
            ImGui::Separator();
            ImGui::TextUnformatted("Chromatic focus shift (focal length vs wavelength):");
            if (!g_res.chrom_focus.empty())
                ImGui::PlotLines("##chrom", g_res.chrom_focus.data(),
                                 static_cast<int>(g_res.chrom_focus.size()), 0,
                                 nullptr, FLT_MAX, FLT_MAX, ImVec2(-1, 120));
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
