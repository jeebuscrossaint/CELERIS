#pragma once
// General complex eigendecomposition used by the RCWA per-layer solve.
//
// Default backend is Eigen's ComplexEigenSolver (header-only, zero dependency --
// keeps CELERIS a single self-contained binary). When built with
// CELERIS_USE_LAPACK (e.g. a statically linked OpenBLAS/MKL), it dispatches to
// LAPACK's zgeev, which is markedly faster for the matrix sizes here -- and the
// exe stays a single file because the LAPACK code is linked in statically.
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
