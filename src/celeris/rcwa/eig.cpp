#include "celeris/rcwa/eig.hpp"

#ifdef CELERIS_USE_LAPACK
// CLAPACK (f2c) LAPACK interface. We declare zgeev_ by hand rather than include
// f2c.h/clapack.h, which #define abs/min/max and clash with Eigen. f2c types:
// integer == 32-bit (long) == int on MSVC; doublecomplex is {double re, im},
// layout-compatible with std::complex<double>; column-major like Eigen.
#include <complex>
#include <vector>

#include <Eigen/Eigenvalues>

extern "C" {
int zgeev_(const char* jobvl, const char* jobvr, const int* n,
           std::complex<double>* a, const int* lda, std::complex<double>* w,
           std::complex<double>* vl, const int* ldvl, std::complex<double>* vr,
           const int* ldvr, std::complex<double>* work, const int* lwork,
           double* rwork, int* info);
}

namespace celeris {

void eig_general(const Eigen::MatrixXcd& A, Eigen::VectorXcd& evals,
                 Eigen::MatrixXcd& evecs) {
    const int n = static_cast<int>(A.rows());
    Eigen::MatrixXcd a = A;  // zgeev overwrites its input
    evals.resize(n);
    evecs.resize(n, n);
    const char jobvl = 'N', jobvr = 'V';
    int info = 0;
    std::vector<double> rwork(static_cast<std::size_t>(2) * n);

    // Workspace query (lwork = -1): optimal size returned in work[0].real().
    std::complex<double> wkopt;
    int lwork = -1;
    zgeev_(&jobvl, &jobvr, &n, a.data(), &n, evals.data(), nullptr, &n,
           evecs.data(), &n, &wkopt, &lwork, rwork.data(), &info);
    lwork = std::max(2 * n, static_cast<int>(wkopt.real()));
    std::vector<std::complex<double>> work(lwork);
    zgeev_(&jobvl, &jobvr, &n, a.data(), &n, evals.data(), nullptr, &n,
           evecs.data(), &n, work.data(), &lwork, rwork.data(), &info);

    if (info != 0) {  // numerical failure -> robust Eigen fallback
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
