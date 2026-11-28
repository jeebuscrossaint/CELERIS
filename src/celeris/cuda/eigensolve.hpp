#pragma once
/// @file
/// @brief Experimental GPU eigensolver (cuSOLVER), single and batched.
///
/// GPU eigensolver backend (CUDA / cuSOLVER). Only built when CELERIS_USE_CUDA
/// is enabled. The RCWA per-layer eigenproblem is a general (non-Hermitian)
/// complex matrix, so we use cuSOLVER's Xgeev. The CPU path (Eigen) remains the
/// cross-platform fallback.
///
/// Opt-in and off by default: measured, it was no faster than the CPU path.

#include <complex>

namespace celeris::cuda {

/// @brief Eigenvalues and right eigenvectors of a general complex n × n matrix (cuSOLVER).
/// @param A  Input matrix, column-major, n² values.
/// @param n  Matrix size.
/// @param w  Output: n eigenvalues.
/// @param vr Output: right eigenvectors, column-major, n² values.
/// @return True on success.
bool geev(const std::complex<double>* A, int n,
          std::complex<double>* w, std::complex<double>* vr);

/// @brief Batched eigendecomposition of many same-size general complex matrices.
///
/// A library sweep solves many independent eigenproblems of one size (one per
/// geometry and wavelength). Batching amortizes the per-call handle and
/// allocation cost and runs several solves at once on separate CUDA streams.
/// @param As      `batch` column-major n × n matrices, contiguous (matrix b at As + b·n²).
/// @param n       Matrix size.
/// @param batch   Number of matrices.
/// @param ws      Output: eigenvalues, n per matrix (set b at ws + b·n).
/// @param vrs     Output: right eigenvectors, column-major (set b at vrs + b·n²).
/// @param streams Number of solves run concurrently.
/// @return True if every solve succeeded.
bool geev_batched(const std::complex<double>* As, int n, int batch,
                  std::complex<double>* ws, std::complex<double>* vrs,
                  int streams = 4);

/// @brief True if a CUDA device is present and usable.
bool available();

/// @brief Name of CUDA device 0, or an empty string if there is none.
///
/// Shown in the UI as the active accelerator.
const char* device_name();

} // namespace celeris::cuda
