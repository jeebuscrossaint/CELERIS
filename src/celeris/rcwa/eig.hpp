#pragma once
// General complex eigendecomposition used by the RCWA per-layer solve.
//
// Default backend is Eigen's ComplexEigenSolver (header-only, zero dependency --
// keeps CELERIS a single self-contained binary). With CELERIS_USE_LAPACK it
// instead calls LAPACK's zgeev (the EIGEN_USE_MKL_ALL build routes this to MKL).
//
// MEASURED REALITY (don't be misled by "the eigensolve dominates"): for these
// sizes zgeev is only ~15% of solve_rcwa_2d -- the real cost is the dense complex
// matmuls / inverses in the operator assembly and S-matrix recursion. So the
// speed levers, in order, are: (1) AVX2 (the default build; ~2.7x free over the
// SSE2 baseline), (2) Intel MKL via EIGEN_USE_MKL_ALL, which multithreads ALL the
// dense BLAS *and* zgeev (up to ~11x over SSE2 at high M). Reference LAPACK alone
// is a wash. See CMakeLists CELERIS_USE_MKL.
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
