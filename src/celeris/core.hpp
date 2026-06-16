#pragma once
// Core types shared across the whole CELERIS engine.
//
// Keep this header tiny and dependency-free: everything in the engine includes
// it, so anything heavy here slows the whole build.

#include <complex>
#include <numbers>

namespace celeris {

// All fields in CELERIS are complex (amplitude + phase). A real refractive
// index n with absorption k is stored as the single complex number n + i*k.
using cdouble = std::complex<double>;

inline constexpr double pi = std::numbers::pi;

// Polarization relative to the plane of incidence:
//   TE (s): electric field perpendicular to the plane of incidence.
//   TM (p): electric field lies in the plane of incidence.
// At normal incidence the two are identical.
enum class Pol { TE, TM };

} // namespace celeris
