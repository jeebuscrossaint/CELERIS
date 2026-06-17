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

// True if a CUDA device is present and usable.
bool available();

} // namespace celeris::cuda
