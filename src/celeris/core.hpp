#pragma once
/// @file
/// @brief Core types shared across the whole CELERIS engine.
///
/// Keep this header tiny and dependency-free: everything in the engine includes
/// it, so anything heavy here slows the whole build.

#include <complex>
#include <numbers>

/// @brief The CELERIS engine: RCWA solvers, metalens design, analysis, and I/O.
namespace celeris {

/// @brief Complex double used for every field, amplitude, and index.
///
/// All fields in CELERIS are complex (amplitude + phase). A real refractive
/// index n with absorption k is stored as the single complex number n + i*k.
using cdouble = std::complex<double>;

/// @brief π as a double.
inline constexpr double pi = std::numbers::pi;

/// @brief Polarization relative to the plane of incidence.
///
/// At normal incidence the two are identical.
enum class Pol {
    TE,  ///< s-polarization: electric field perpendicular to the plane of incidence.
    TM   ///< p-polarization: electric field lies in the plane of incidence.
};

} // namespace celeris
