#pragma once
// Scattering-matrix (Redheffer) machinery shared by the 1D and 2D RCWA
// solvers. A block S-matrix maps incoming -> outgoing modal amplitudes at the
// two faces of a region. Because layer propagation uses only decaying
// exponentials, chaining S-matrices with the Redheffer star product stays
// numerically stable for arbitrarily thick or evanescent layers.
//
// The interface/propagation builders are dimension-agnostic: they work for any
// field modes W and companion modes V (1D scalar fields or 2D vectorial
// [Ex;Ey] / [Hx;Hy] blocks), so the 2D solver reuses them unchanged.

#include <Eigen/Dense>

namespace celeris::detail {

using Eigen::MatrixXcd;

struct SMatrix {
    MatrixXcd S11, S12, S21, S22;
};

inline SMatrix identity_smatrix(int n) {
    return {MatrixXcd::Zero(n, n), MatrixXcd::Identity(n, n),
            MatrixXcd::Identity(n, n), MatrixXcd::Zero(n, n)};
}

// Redheffer star product  A ⋆ B  (A above, B below).
inline SMatrix star(const SMatrix& A, const SMatrix& B) {
    const int n = static_cast<int>(A.S11.rows());
    const MatrixXcd Id = MatrixXcd::Identity(n, n);
    const MatrixXcd D = (Id - B.S11 * A.S22).inverse();
    const MatrixXcd F = (Id - A.S22 * B.S11).inverse();
    SMatrix S;
    S.S11 = A.S11 + A.S12 * D * B.S11 * A.S21;
    S.S12 = A.S12 * D * B.S12;
    S.S21 = B.S21 * F * A.S21;
    S.S22 = B.S22 + B.S21 * F * A.S22 * B.S12;
    return S;
}

// Interface between region "a" (above, modes Wa,Va) and "b" (below, Wb,Vb),
// amplitudes referenced at the interface. From continuity of field (W) and
// companion (V):  Wa(u_a+ + u_a-) = Wb(u_b+ + u_b-),
//                 Va(u_a+ - u_a-) = Vb(u_b+ - u_b-).
inline SMatrix interface_smatrix(const MatrixXcd& Wa, const MatrixXcd& Va,
                                 const MatrixXcd& Wb, const MatrixXcd& Vb) {
    const int n = static_cast<int>(Wa.rows());
    MatrixXcd Mo(2 * n, 2 * n), Mi(2 * n, 2 * n);
    Mo.block(0, 0, n, n) = Wa;
    Mo.block(0, n, n, n) = -Wb;
    Mo.block(n, 0, n, n) = Va;
    Mo.block(n, n, n, n) = Vb;
    Mi.block(0, 0, n, n) = -Wa;
    Mi.block(0, n, n, n) = Wb;
    Mi.block(n, 0, n, n) = Va;
    Mi.block(n, n, n, n) = Vb;

    MatrixXcd Sfull = Mo.partialPivLu().solve(Mi);
    return {Sfull.block(0, 0, n, n), Sfull.block(0, n, n, n),
            Sfull.block(n, 0, n, n), Sfull.block(n, n, n, n)};
}

// Propagation across a layer with modal decay factor X (applied both ways).
inline SMatrix propagation_smatrix(const MatrixXcd& X) {
    const int n = static_cast<int>(X.rows());
    return {MatrixXcd::Zero(n, n), X, X, MatrixXcd::Zero(n, n)};
}

} // namespace celeris::detail
