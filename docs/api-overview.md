# CELERIS C++ API {#mainpage}

CELERIS is a validated RCWA (rigorous coupled-wave analysis) pipeline for
metalens and metasurface design. The engine library, `celeris_core`, is plain
C++23. Every routine takes parameters and returns plain data structs, so the
CLI, the Python module, and the GUI all call the same functions.

## Conventions

- Lengths are in micrometres (`_um` suffix), angles in radians (`_rad`).
- Wavelengths are vacuum wavelengths.
- Refractive indices and fields are complex (`celeris::cdouble`); an index
  n + ik carries absorption in k.
- Fourier truncation is set by half-counts: `M` keeps orders −M..M.

## Modules

| Header | What it provides |
|---|---|
| `celeris/core.hpp` | Shared types: `cdouble`, `pi`, `Pol` (TE/TM). |
| `celeris/optics/fresnel.hpp` | Single-interface Fresnel coefficients. |
| `celeris/optics/tmm.hpp` | Transfer-matrix solver for planar multilayers (the closed-form reference). |
| `celeris/rcwa/grating1d.hpp` | 1D binary grating layers and their Fourier coefficients. |
| `celeris/rcwa/rcwa1d.hpp` | 1D RCWA, single layer and multilayer (S-matrix), TE and TM. |
| `celeris/rcwa/rcwa2d.hpp` | 2D vectorial RCWA for biperiodic meta-atom cells. |
| `celeris/rcwa/eig.hpp` | The complex eigensolver backend (Eigen or LAPACK/MKL). |

The design, analysis, materials, and I/O modules are being added to this
reference (issue #1). Until then, their headers under `src/celeris/` carry
plain-comment documentation.

## A minimal 2D solve

```cpp
#include "celeris/rcwa/rcwa2d.hpp"
#include "celeris/materials/database.hpp"

using namespace celeris;

Material tio2 = Material::constant({2.40, 0.0}, "TiO2");
Material air = materials::air();

Rcwa2DStack stack{0.35, 0.35, {RectCell2D{tio2, air, 0.5, 0.5, 0.6}}};
Rcwa2DResult r = solve_rcwa_2d(air, stack, materials::fused_silica(),
                               /*wavelength_um=*/0.532, /*theta*/ 0.0, /*phi*/ 0.0,
                               /*Ex0*/ 1.0, /*Ey0*/ 0.0, /*Mx*/ 8, /*My*/ 8);
double phase = std::arg(r.tx0);   // pillar phase delay
```

## Validation

The physics behind these functions is locked by `celeris selftest`
(energy conservation, 2D = 1D, RCWA = TMM, agreement with grcwa and S⁴).
