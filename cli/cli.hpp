#pragma once

// Shared front-end header for the CELERIS CLI. Each subcommand is implemented in
// its own cli/*.cpp translation unit; this header carries the common engine
// includes, the shared argument/material/profile helpers, and the declarations of
// every subcommand entry point dispatched from main.cpp. `using namespace celeris`
// is intentional here: these are all internal CLI TUs, never a public API.

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <format>
#include <fstream>
#include <optional>
#include <print>
#include <string>
#include <vector>

// Either GPU path (cuSOLVER eigensolve = CELERIS_USE_CUDA, or the far-field
// kernel = CELERIS_USE_CUDA_KERNELS) needs these; the headers are independent so
// the kernel can build without the (dead-end, opt-in) cuSOLVER eigensolve.
#if defined(CELERIS_USE_CUDA) || defined(CELERIS_USE_CUDA_KERNELS)
#include <Eigen/Dense>
#include <algorithm>
#include <future>
#include <random>
#include <thread>
#endif
#ifdef CELERIS_USE_CUDA
#include "celeris/cuda/eigensolve.hpp"
#endif
#ifdef CELERIS_USE_CUDA_KERNELS
#include "celeris/cuda/propagate.hpp"
#endif

#include "celeris/analysis/chromatic.hpp"
#include "celeris/analysis/efficiency.hpp"
#include "celeris/analysis/field.hpp"
#include "celeris/analysis/focal.hpp"
#include "celeris/analysis/polarization.hpp"
#include "celeris/analysis/throughfocus.hpp"
#include "celeris/analysis/tolerance.hpp"
#include "celeris/analysis/wavefront.hpp"
#include "celeris/design/achromatic.hpp"
#include "celeris/design/metalens.hpp"
#include "celeris/design/optimize.hpp"
#include "celeris/design/pb_achromatic.hpp"
#include "celeris/design/pb_metalens.hpp"
#include "celeris/design/polar_metalens.hpp"
#include "celeris/io/gds.hpp"
#include "celeris/io/image.hpp"
#include "celeris/io/material_csv.hpp"
#include "celeris/materials/database.hpp"
#include "celeris/optics/tmm.hpp"
#include "celeris/rcwa/rcwa1d.hpp"
#include "celeris/rcwa/rcwa2d.hpp"

using namespace celeris;

// Group delay (fs) of a free-space optical path of length L (um): GD = L/c.
inline constexpr double GD_FS_PER_UM = 1e9 / 2.99792458e8;  // ~3.3356 fs per um

// Where to render the reconstruction image, plus a one-line optical result.
struct ProfileProof {
    std::string summary;
    double psf_cx = 0, psf_cy = 0, psf_z = 0, psf_hw = 0;
};

// --- shared CLI helpers (defined in cli/cli_common.cpp) ---------------------

// Value following `key` in argv (searched from index 2), or `fallback`.
const char* arg_value(int argc, char** argv, const std::string& key,
                      const char* fallback);

// Resolve the substrate material from --substrate <name> (registry name).
Material resolve_substrate(int argc, char** argv, const std::string& def = "bk7");

// Resolve the pillar material: --pillar-csv > --pillar <name> > --pillar-n.
Material resolve_pillar(int argc, char** argv, double def_n = 2.4);

// Parse --profile (+ its parameters) into a PhaseProfile.
std::optional<PhaseProfile> parse_phase_profile(int argc, char** argv,
                                                double focal, double diameter);

// Propagate a realized aperture field and run the profile's optical proof.
ProfileProof profile_optical_proof(const std::vector<double>& px,
                                   const std::vector<double>& py,
                                   const std::vector<cdouble>& tc,
                                   const PhaseProfile& profile, double lambda,
                                   double focal, double diameter, double recon_z);

// --- subcommand entry points (one per cli/*.cpp) ----------------------------

int run_selftest(bool quick = false);
int cmd_design(int argc, char** argv);
int cmd_widefov(int argc, char** argv);
int cmd_efficiency(int argc, char** argv);
int cmd_fieldmap(int argc, char** argv);
int cmd_birefringence(int argc, char** argv);
int cmd_polardesign(int argc, char** argv);
int cmd_pbdesign(int argc, char** argv);
int cmd_achromatic(int argc, char** argv);
int cmd_pb_achromatic(int argc, char** argv);
int cmd_reproduce(int argc, char** argv);
int cmd_validate(int argc, char** argv);
int cmd_materials(int argc, char** argv);
int cmd_shapeconv(int argc, char** argv);
void print_help();
int run_gpubench(int argc, char** argv);
int run_psfbench(int argc, char** argv);
