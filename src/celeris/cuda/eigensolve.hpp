#pragma once
// GPU eigensolver backend (CUDA / cuSOLVER). Only built when CELERIS_USE_CUDA
// is enabled. The RCWA per-layer eigenproblem is a general (non-Hermitian)
// complex matrix, so we use cuSOLVER's Xgeev. The CPU path (Eigen) remains the
// cross-platform fallback.

#include <complex>

namespace celeris::cuda {

// Right-eigenvector decomposition of a general complex(double) n×n matrix,
// stored column-major. Writes n eigenvalues to `w` and the n×n column-major
// right eigenvectors to `vr`. Returns true on success.
bool geev(const std::complex<double>* A, int n,
          std::complex<double>* w, std::complex<double>* vr);

// BATCHED general complex eigendecomposition — the real GPU play for RCWA.
// A metalens library sweep solves many independent same-size eigenproblems (one
// per pillar geometry / wavelength); batching amortizes the per-call handle and
// allocation overhead that made single GPU solves slow, and runs several at once
// across CUDA streams to fill the device.
//
// `As` points to `batch` column-major n×n matrices laid out contiguously
// (matrix b at As + b*n*n). Writes `batch` eigenvalue sets to `ws` (b at
// ws + b*n) and the right eigenvectors to `vrs` (b at vrs + b*n*n). `streams`
// sets how many solves run concurrently. Returns true if every solve succeeded.
bool geev_batched(const std::complex<double>* As, int n, int batch,
                  std::complex<double>* ws, std::complex<double>* vrs,
                  int streams = 4);

// True if a CUDA device is present and usable.
bool available();

} // namespace celeris::cuda
