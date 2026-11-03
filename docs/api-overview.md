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
| `celeris/design/phase_profile.hpp` | Target phase profiles: focusing, quadratic, vortex, deflector, axicon, freeform. |
| `celeris/design/metalens.hpp` | Meta-atom phase libraries, full-2π height search, propagation-phase lens design. |
| `celeris/design/optimize.hpp` | Single-pillar inverse design (Adam, finite-difference gradients). |
| `celeris/design/system_opt.hpp` | Period × height search on a focal-Strehl merit. |
| `celeris/design/polar_metalens.hpp` | Polarization-multiplexed lenses from rectangular pillars. |
| `celeris/design/pb_metalens.hpp` | Pancharatnam-Berry (geometric-phase) lenses: Jones matrices, half-wave-plate atom search, rotation maps. |
| `celeris/design/achromatic.hpp` | Achromatic lenses by dispersion engineering: band-characterized libraries (fill × height or single-etch), two-objective selection, focal-shift check. |
| `celeris/design/pb_achromatic.hpp` | Achromatic geometric-phase lenses: single-etch birefringent libraries chosen for group delay alone. |
| `celeris/analysis/focal.hpp` | Focal-plane metrics (Strehl, phase Strehl, FWHM, encircled energy), PSF maps, Rayleigh-Sommerfeld propagation (CPU and GPU), off-axis PSF. |
| `celeris/analysis/throughfocus.hpp` | On-axis through-focus scan, depth of focus, x-z caustic. |
| `celeris/analysis/wavefront.hpp` | OPD map, RMS/PV wavefront error, Maréchal Strehl, Zernike fit. |
| `celeris/analysis/mtf.hpp` | MTF from the pupil autocorrelation, with the diffraction-limited reference. |
| `celeris/analysis/field.hpp` | Field-of-view sweeps (stop at the lens or offset stop) and the full 2D field grid. |
| `celeris/analysis/chromatic.hpp` | Focal shift vs wavelength, fast (fixed transmission) or rigorous (library rebuilt per wavelength). |
| `celeris/analysis/tolerance.hpp` | Monte-Carlo Strehl under Gaussian pillar-size (CD) error. |
| `celeris/analysis/efficiency.hpp` | Per-order power budget of one meta-atom: zeroth order, stray light, absorption. |
| `celeris/analysis/polarization.hpp` | Form birefringence of rectangular pillars: co-polarized phases and retardance. |

The materials and I/O modules are being added to this
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
