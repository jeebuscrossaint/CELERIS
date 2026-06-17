#pragma once
// GPU far-field propagation — the part of RCWA-metalens analysis that the GPU
// genuinely wins. Building the focal-plane field is a Rayleigh-Sommerfeld sum
//   E(fx,fy) = sum_pillars t * exp(i k R) / R,   R = sqrt((fx-x)^2+(fy-y)^2+z^2)
// i.e. an N-body-style dense reduction over (pixels x pillars) with no
// eigensolve. For large apertures (hundreds of thousands of pillars) this is
// billions of sincos/sqrt evaluations — embarrassingly parallel and bandwidth-
// friendly, the opposite of the serial nonsymmetric eigensolve.
//
// Computed in single precision on the device (the PSF feeds Strehl/FWHM/imaging,
// where float is ample) and returned as double. Pure C++ signature so callers
// need no CUDA headers.

#include <complex>

namespace celeris::cuda {

// Fill `out` (n*n doubles, row-major) with |E|^2 on an n x n focal-plane grid
// spanning [cx-W, cx+W] x [cy-W, cy+W] at axial distance z, summing over `npil`
// pillars (px,py in microns, pt complex transmission). Returns false if no CUDA
// device is available or a CUDA call fails (caller should fall back to CPU).
bool propagate_psf(const double* px, const double* py,
                   const std::complex<double>* pt, int npil,
                   double cx, double cy, double z, double k,
                   int n, double half_window, double* out);

} // namespace celeris::cuda
