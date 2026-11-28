#pragma once
/// @file
/// @brief GPU far-field propagation: focal-plane PSF and x-z caustic.
///
/// GPU far-field propagation — the part of RCWA-metalens analysis that the GPU
/// genuinely wins. Building the focal-plane field is a Rayleigh-Sommerfeld sum
///   E(fx,fy) = sum_pillars t * exp(i k R) / R,   R = sqrt((fx-x)^2+(fy-y)^2+z^2)
/// i.e. an N-body-style dense reduction over (pixels x pillars) with no
/// eigensolve. For large apertures (hundreds of thousands of pillars) this is
/// billions of sincos/sqrt evaluations — embarrassingly parallel and bandwidth-
/// friendly, the opposite of the serial nonsymmetric eigensolve.
///
/// Computed in single precision on the device (the PSF feeds Strehl/FWHM/imaging,
/// where float is ample) and returned as double. Pure C++ signature so callers
/// need no CUDA headers.

#include <complex>

/// @brief CUDA backends: far-field propagation and the experimental eigensolver.
///
/// Built only with `CELERIS_USE_CUDA`. Every function returns false (or an
/// empty result) when no device is available, so callers fall back to the CPU.
namespace celeris::cuda {

// GPU presence / identity (cudart only, no cuSOLVER). Defined here in the kernel
// TU so the kernel build has them without the opt-in cuSOLVER eigensolve unit.

/// @brief True if a CUDA device is present and usable.
bool available();
/// @brief Name of CUDA device 0, or an empty string if there is none.
const char* device_name();

/// @brief Focal-plane intensity |E|² on the GPU by a Rayleigh–Sommerfeld sum over pillars.
///
/// Computes on an n × n grid spanning [cx − W, cx + W] × [cy − W, cy + W] at
/// axial distance z, in single precision on the device.
/// @param px          Pillar x positions (µm), `npil` values.
/// @param py          Pillar y positions (µm), `npil` values.
/// @param pt          Pillar complex transmissions, `npil` values.
/// @param npil        Number of pillars.
/// @param cx          Grid centre, x (µm).
/// @param cy          Grid centre, y (µm).
/// @param z           Axial distance from the lens (µm).
/// @param k           Wavenumber 2π/λ (rad/µm).
/// @param n           Grid points per side.
/// @param half_window Half-width W of the grid (µm).
/// @param out         Output: n² intensities, row-major.
/// @return False if no CUDA device is available or a CUDA call fails; the caller
///         should then use the CPU path.
bool propagate_psf(const double* px, const double* py,
                   const std::complex<double>* pt, int npil,
                   double cx, double cy, double z, double k,
                   int n, double half_window, double* out);

/// @brief Longitudinal caustic |E|² in the y = 0 plane on the GPU.
///
/// Shows the focusing cone and depth of focus.
/// @param px   Pillar x positions (µm), `npil` values.
/// @param py   Pillar y positions (µm), `npil` values.
/// @param pt   Pillar complex transmissions, `npil` values.
/// @param npil Number of pillars.
/// @param xmin Smallest x (µm).
/// @param xmax Largest x (µm).
/// @param nx   Points along x.
/// @param zmin Smallest axial distance (µm).
/// @param zmax Largest axial distance (µm).
/// @param nz   Points along z.
/// @param k    Wavenumber 2π/λ (rad/µm).
/// @param out  Output: nz · nx intensities, row-major (row = z index).
/// @return False on failure.
bool propagate_zx(const double* px, const double* py,
                  const std::complex<double>* pt, int npil, double xmin,
                  double xmax, int nx, double zmin, double zmax, int nz, double k,
                  double* out);

} // namespace celeris::cuda
