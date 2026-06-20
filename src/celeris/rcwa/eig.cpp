#include "celeris/rcwa/eig.hpp"

#ifdef CELERIS_USE_LAPACK
// CLAPACK (f2c) LAPACK interface. We declare zgeev_ by hand rather than include
// f2c.h/clapack.h, which #define abs/min/max and clash with Eigen. f2c types:
// integer == 32-bit (long) == int on MSVC; doublecomplex is {double re, im},
// layout-compatible with std::complex<double>; column-major like Eigen.
#include <complex>
#include <cstdio>
#include <vector>

#include <Eigen/Eigenvalues>

#ifdef EIGEN_USE_MKL_ALL
// MKL build: <mkl.h> (pulled in by Eigen's MKL support) already declares zgeev_
// with MKL_Complex16, so we must NOT redeclare it -- just use it. MKL_Complex16
// is layout-compatible with std::complex<double>, so the pointers reinterpret.
#include <mkl.h>
using lpk_cd = MKL_Complex16;
#else
// Standalone LAPACK (e.g. reference/CLAPACK): declare zgeev_ by hand rather than
// include clapack.h, which #defines abs/min/max and clashes with Eigen. f2c
// integer == int on MSVC; doublecomplex {re,im} is layout-compatible with
// std::complex<double>; column-major like Eigen.
using lpk_cd = std::complex<double>;
extern "C" {
int zgeev_(const char* jobvl, const char* jobvr, const int* n, lpk_cd* a,
           const int* lda, lpk_cd* w, lpk_cd* vl, const int* ldvl, lpk_cd* vr,
           const int* ldvr, lpk_cd* work, const int* lwork, double* rwork, int* info);
}
#endif

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
    auto cx = [](std::complex<double>* p) { return reinterpret_cast<lpk_cd*>(p); };

    // Workspace query (lwork = -1): optimal size returned in work[0].real().
    std::complex<double> wkopt{};
    int lwork = -1;
    zgeev_(&jobvl, &jobvr, &n, cx(a.data()), &n, cx(evals.data()), nullptr, &n,
           cx(evecs.data()), &n, cx(&wkopt), &lwork, rwork.data(), &info);
    lwork = std::max(2 * n, static_cast<int>(wkopt.real()));
    std::vector<std::complex<double>> work(lwork);
    zgeev_(&jobvl, &jobvr, &n, cx(a.data()), &n, cx(evals.data()), nullptr, &n,
           cx(evecs.data()), &n, cx(work.data()), &lwork, rwork.data(), &info);

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
