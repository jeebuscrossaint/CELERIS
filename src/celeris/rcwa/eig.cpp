#include "celeris/rcwa/eig.hpp"

#ifdef CELERIS_USE_LAPACK
// Make LAPACKE's complex type std::complex<double> so Eigen buffers pass through
// with no reinterpret/copy. Must be defined before <lapacke.h>.
#define LAPACK_COMPLEX_CPP
#include <lapacke.h>

#include <Eigen/Eigenvalues>

namespace celeris {

void eig_general(const Eigen::MatrixXcd& A, Eigen::VectorXcd& evals,
                 Eigen::MatrixXcd& evecs) {
    const lapack_int n = static_cast<lapack_int>(A.rows());
    Eigen::MatrixXcd a = A;  // zgeev overwrites its input; Eigen is column-major
    evals.resize(n);
    evecs.resize(n, n);
    lapack_int info = LAPACKE_zgeev(
        LAPACK_COL_MAJOR, 'N', 'V', n, a.data(), n, evals.data(),
        /*vl=*/nullptr, n, evecs.data(), n);
    if (info != 0) {  // numerical failure -> fall back to the robust Eigen path
        Eigen::ComplexEigenSolver<Eigen::MatrixXcd> ces(A);
        evals = ces.eigenvalues();
        evecs = ces.eigenvectors();
    }
}

} // namespace celeris

#else
#include <Eigen/Eigenvalues>

namespace celeris {

void eig_general(const Eigen::MatrixXcd& A, Eigen::VectorXcd& evals,
                 Eigen::MatrixXcd& evecs) {
    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> ces(A);
    evals = ces.eigenvalues();
    evecs = ces.eigenvectors();
}

} // namespace celeris
#endif
