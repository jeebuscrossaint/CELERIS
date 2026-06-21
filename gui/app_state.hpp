#pragma once
// Shared GUI state: the parameter/result data model, the cross-thread globals
// the worker threads and the render loop communicate through, and the worker
// entry points. All in namespace celeris::gui.
#include <atomic>
#include <mutex>
#include <optional>
#include <string>
#include <vector>

#include "celeris/analysis/field.hpp"
#include "celeris/analysis/focal.hpp"
#include "celeris/analysis/mtf.hpp"
#include "celeris/analysis/throughfocus.hpp"
#include "celeris/analysis/tolerance.hpp"
#include "celeris/analysis/wavefront.hpp"
#include "celeris/design/metalens.hpp"
#include "celeris/design/polar_metalens.hpp"
#include "celeris/design/system_opt.hpp"
#include "celeris/materials/material.hpp"

namespace celeris::gui {

struct LayerRow {  // an extra (unpatterned by default) layer above the pillars
    float n = 1.46f;
    float fill = 1.0f;     // 1.0 = uniform film; <1 = a square patch
    float thickness = 0.1f;
};

struct Params {
    float focal = 50, diameter = 20, wavelength = 0.532f, period = 0.35f;
    float focal_y = 80;  // Y-polarization focal length (polarization-multiplexed mode)
    float thickness = 0.6f, pillar_n = 2.4f;  // the active (patterned) layer
    int harmonics = 6, fill_samples = 18;
    int pillar_mat = 3;     // index into kPillarMats (default: TiO2 approx)
    int substrate_mat = 0;  // 0 = N-BK7, 1 = air, 2 = fused silica
    std::vector<LayerRow> extra_layers;  // stacked above the pillars (cap/AR…)
};

struct Results {
    double strehl = 0, phase_strehl = 0, fwhm = 0, dl = 0, encircled = 0, rms = 0, meanT = 0;
    double coverage_deg = 0, na = 0;
    int n_cells = 0, pillars = 0;
    PsfMap psf;
    std::vector<float> chrom_wl, chrom_focus;
    WavefrontAnalysis wf;
    MtfCurve mtf;
    ThroughFocus tf;
    MetalensDesign design;
    UnitCellLibrary lib;
    Params used;
};

inline constexpr const char* kPillarMats[] = {
    "Custom (constant n)", "Silicon nitride (Si3N4)", "Fused silica (SiO2)",
    "TiO2 (approx n=2.40)", "a-Si (approx n=3.50)", "GaN (approx n=2.35)",
    "Loaded CSV"};
inline constexpr const char* kSubstrates[] = {"N-BK7", "Air", "Fused silica (SiO2)"};

// --- cross-thread globals (defined in app_state.cpp) ------------------------
extern std::optional<Material> g_loaded_material;
extern std::string g_loaded_name;

extern std::atomic<bool> g_running, g_pending, g_tol_pending, g_fov_pending,
    g_spot_pending, g_polar_pending, g_opt_pending;
extern std::atomic<float> g_progress;
extern std::mutex g_mtx;
extern std::string g_status;
extern Results g_res;

extern std::vector<ToleranceResult> g_tol;
extern std::vector<FieldPoint> g_fov;
extern std::vector<FieldPsf> g_spot;

extern PolarMetalensDesign g_polar;
extern double g_polar_fx, g_polar_fy, g_polar_zx, g_polar_zy, g_polar_iso_x,
    g_polar_iso_y;
extern PsfMap g_polar_psf_x, g_polar_psf_y;

extern SystemOptResult g_opt;

// --- model helpers + worker entry points ------------------------------------
Material make_pillar(const Params& p);
const Material& make_substrate(const Params& p);
bool save_project(const std::string& path, const Params& p);
bool load_project(const std::string& path, Params& p);

void set_phase(const char* msg, float progress);
void run_design(Params p);
void run_tolerance();
void run_fov();
void run_spotgrid();
void run_polardesign(Params p);
void run_optimize(Params p);

} // namespace celeris::gui
