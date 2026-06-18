#pragma once
// General complex eigendecomposition used by the RCWA per-layer solve.
//
// Default backend is Eigen's ComplexEigenSolver (header-only, zero dependency --
// keeps CELERIS a single self-contained binary, and measured the fastest option
// readily available on Windows here). When built with CELERIS_USE_LAPACK it
// instead calls LAPACK's zgeev (statically linked, so still a single exe). NOTE:
// only a *re-optimized* LAPACK (Intel MKL) actually beats Eigen for these sizes
// -- reference LAPACK (CLAPACK/f2c) is ~1.5x slower, since the cost is the zgeev
// algorithm itself, not the BLAS. The shim exists so MKL can be dropped in.
//
// RCWA physics is invariant to per-eigenvector column scaling and to eigenvalue
// ordering (both are absorbed into the modal coefficients), so the two backends
// produce equivalent diffraction efficiencies / transmission.

#include <Eigen/Dense>

namespace celeris {

// Right-eigenvector decomposition: A * evecs.col(i) = evals(i) * evecs.col(i).
void eig_general(const Eigen::MatrixXcd& A, Eigen::VectorXcd& evals,
                 Eigen::MatrixXcd& evecs);

} // namespace celeris
