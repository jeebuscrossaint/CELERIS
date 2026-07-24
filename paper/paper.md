---
title: "CELERIS: an integrated, validated RCWA pipeline for metalens design, analysis, and fabrication layout"
# Target: Computer Physics Communications (CPC) + arXiv (physics.comp-ph / physics.optics).
#   Fallback: SoftwareX. (SciPost dropped: it bars undergrad solo submission.)
#   Convert to Elsevier `elsarticle` LaTeX for submission (NOT SciPost.cls).
# STATUS: DRAFT. [TODO] markers flag numbers/figures to generate and claims to verify
# before submission. Do not submit with any [TODO] remaining.
author:
  - name: Amarnath Patel
    orcid: 0009-0008-9460-082X
    affiliation: "University of Central Florida, Orlando, FL, USA"
date: DRAFT
---

## Program summary

*(CPC-mandatory block — fill exact fields on the CPC submission form.)*

- **Program title:** CELERIS
- **Licensing provisions:** Apache License 2.0
- **Programming languages:** C++23 (engine), CUDA (optional GPU far-field kernel), Python
  (pybind11 bindings)
- **Nature of problem:** Designing a metalens/metasurface requires solving Maxwell's
  equations for periodic subwavelength unit cells, assembling a meta-atom library, mapping
  a target optical phase profile onto it, analyzing the resulting device's optical
  performance, and exporting a fabrication mask — a multi-stage workflow that existing RCWA
  solvers do not provide end to end.
- **Solution method:** Rigorous coupled-wave analysis / Fourier modal method (1D and
  2D-vectorial) with Li/Liu–Fan Fourier factorization and a Redheffer S-matrix layer
  recursion, wrapped in a library builder, phase-mapping/inverse-design layer (including
  achromatic and Pancharatnam–Berry recipes), an optical analysis battery, GPU-accelerated
  far-field propagation, and GDSII layout export.
- **Dependencies:** Eigen (header-only, fetched at configure time); optional Intel MKL and
  CUDA for acceleration; pybind11 for the Python module. Validated against `grcwa` and
  Stanford S$^4$.
- **Restrictions:** Scope is metalenses/metasurfaces via RCWA (not FDTD or ray tracing). The
  modal eigensolve runs on the CPU; the GPU accelerates far-field propagation only.
- **Operating systems:** Linux and Windows (CPU build has no external runtime dependencies).

## Summary

Metasurfaces and metalenses — flat optical elements built from subwavelength
scatterers ("meta-atoms") — are designed by (i) solving Maxwell's equations for a
periodic unit cell to build a library of meta-atom responses, (ii) mapping a target
optical phase profile onto that library, (iii) analyzing the resulting device's optical
performance, and (iv) exporting a mask layout for fabrication. The standard rigorous
electromagnetic method for step (i) is **rigorous coupled-wave analysis (RCWA)**, also
known as the Fourier modal method (FMM).

**CELERIS** is a from-scratch, validated implementation of this *entire* workflow in a
single native package. It couples a C++23 RCWA engine (1D and 2D-vectorial, with
Li/Liu–Fan Fourier factorization and a Redheffer S-matrix stack) to a meta-atom library
builder, a phase-mapping and inverse-design layer (including achromatic and
Pancharatnam–Berry designs), a full optical analysis battery
(Strehl ratio, PSF, MTF, wavefront/Zernike, through-focus, chromatic, field-resolved,
and tolerance analyses), GPU-accelerated far-field propagation, and fabrication-ready
GDSII export. It is scriptable from C++, a command-line interface, and Python
(via pybind11), and ships with a desktop GUI. Correctness is enforced by a locked
self-test suite (20+ cases) cross-checked against closed-form physics (Fresnel/TMM,
effective-medium theory), energy conservation, and an independent RCWA solver (`grcwa`),
and the full pipeline is demonstrated by reproducing two canonical published TiO$_2$
metalenses end to end.

## Statement of need

Existing open-source RCWA codes — S$^4$, RETICOLO, `grcwa`, TORCWA, `fmmax`, `meent` —
are **solver kernels**: they compute the electromagnetic response of a given periodic
structure. Turning a solver into a *fabricated, analyzed metalens* still requires the
user to hand-assemble the surrounding workflow: sweeping a meta-atom library, mapping a
phase profile onto it, propagating to a focal plane, computing imaging metrics, running
fabrication-tolerance studies, and emitting a GDSII layout. Each researcher re-implements
this scaffolding, and the glue code is rarely validated or shared.

CELERIS provides that pipeline as a single validated tool. The target audience is
metasurface and metalens researchers who need to go from an optical specification
(focal length, aperture, wavelength, material) to a fabrication-ready, performance-
characterized design without stitching together a solver, a plotting stack, and a
layout tool. The contribution is not a new solver kernel but the **integration** —
solver → library → design → analysis → layout — validated end to end against published
devices.

## State of the field

Table 1 positions CELERIS against widely used RCWA/FMM codes. CELERIS deliberately
*does not* compete on automatic differentiation (see Limitations); it competes on the
integrated design–analysis–fabrication workflow, which the solver kernels do not
provide.

**Table 1.** Feature comparison (`✓` = provided, `~` = partial/possible with user code,
`✗` = not provided). Verified against each tool's documentation/publication; released
versions evolve, so cells should be re-checked at submission.

| Feature | S$^4$ | RETICOLO | grcwa | TORCWA | fmmax | meent | **CELERIS** |
|---|---|---|---|---|---|---|---|
| Language | C++/Lua/Py | MATLAB | Python | PyTorch | JAX | np/JAX/torch | **C++/Py** |
| License | GPL | (MATLAB) | GPL | — | MIT | MIT | **Apache-2.0** |
| 1D + 2D-vectorial RCWA | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | **✓** |
| Native GPU (full solve) | ✗ | ✗ | ✗ | ✓ | ✓ | ✓ | **✗ (GPU far-field only)** |
| Automatic differentiation | ✗ | ✗ | ✓ | ✓ | ✓ | ✓ | **✗ (FD; adjoint planned)** |
| Arbitrary meta-atom shapes | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | **✓** |
| Integrated metalens **design** | ✗ | ✗ | ~ | ~ | ~ | ~ | **✓** |
| Optical **analysis battery** | ✗ | ~ | ✗ | ✗ | ✗ | ✗ | **✓** |
| Achromatic + PB design | ✗ | ✗ | ✗ | ✗ | ✗ | ✗ | **✓** |
| Fabrication-ready **GDSII** export | ✗ | ✗ | ✗ | ✗ | ✗ | ✗ | **✓** |
| Published-device reproductions | — | — | — | — | — | — | **✓** |
| Desktop GUI | ✗ | ✗ | ✗ | ✗ | ✗ | ✗ | **✓** |

We state the trade-off plainly: the GPU/autodiff-native solvers (TORCWA, fmmax, meent)
are *stronger solver kernels* — they run the full modal solve on the GPU with automatic
differentiation, whereas CELERIS GPU-accelerates only far-field propagation and
optimizes by finite differences (an analytic adjoint is planned). CELERIS's contribution
is orthogonal: it is the only tool that wraps the solver in a complete, validated
**design → analysis → fabrication** workflow — integrated design, the optical analysis
battery, achromatic/Pancharatnam–Berry recipes, GDSII export, published-device
reproductions, and a GUI — none of which the solver kernels provide.

## Software design

- **Engine (`celeris_core`, C++23).** 1D RCWA (TE/TM, multilayer, S-matrix) and
  2D-vectorial RCWA. The 2D solver uses the Li/Liu–Fan inverse-rule Fourier
  factorization for axis-aligned rectangles (correct convergence, energy conserved to
  machine precision) and a stable Laurent factorization for arbitrary sampled shapes.
  Layer coupling uses a Redheffer S-matrix recursion. The dense complex linear algebra
  runs on Eigen (header-only), optionally routed through Intel MKL.
- **Design.** A library builder sweeps meta-atom geometry; phase profiles
  (focusing, vortex/OAM, deflector, axicon, quadratic wide-FOV, freeform) drive both a
  propagation-phase path and an exact geometric (Pancharatnam–Berry) path. Achromatic
  design engineers group delay across a band via dispersive libraries (fill×height,
  single-etch shape variety, and a PB variant).
- **Analysis.** Strehl/FWHM/encircled energy, PSF, wavefront (OPD/Zernike/Maréchal),
  MTF, through-focus/caustic, chromatic, spot-vs-field, full-field grid, tolerance,
  birefringence, and per-order efficiency/absorption budgets.
- **Acceleration.** The far-field Rayleigh–Sommerfeld propagation — a dense
  (pixels × pillars) reduction — runs as a CUDA kernel (static runtime; CPU auto-
  fallback when no device is present). **The RCWA eigensolve itself remains on the CPU**
  (see Limitations). The CUDA build is cross-platform (Linux and Windows).
- **Interfaces.** A CLI (`design`, `pbdesign`, `achromatic`, `reproduce`, `selftest`,
  …), Python bindings (pybind11; `pip install .`), and a Dear ImGui desktop GUI.

## Validation and research impact

Correctness is locked by a self-test suite of 20+ cases, cross-checked against:
closed-form Fresnel/TMM and Bragg-mirror results; effective-medium theory in the deep-
subwavelength limit; strict energy conservation ($R+T+A = 1$ to $10^{-6}$); and two
independent RCWA solvers, `grcwa` and **Stanford S$^4$** (the field-standard FMM
reference). The suite runs green on Linux (GCC) and is enforced by continuous
integration.

Running CELERIS and S$^4$ on an *identical* freestanding $n=1.5$ binary grating
($\Lambda = 0.3\,\mu$m, fill 0.5, thickness $0.5\,\mu$m, $\lambda = 0.5\,\mu$m, normal
incidence), the zeroth-order transmittance agrees to **$7\times10^{-8}$ (TE)** at $M=40$
(converging $10^{-5}\to10^{-7}$ with Fourier order), and to $3\times10^{-4}$ (TM), the
latter limited by the well-known slower TM Fourier convergence exhibited by *both*
solvers. This case also equals grcwa's 0.93333 to $\sim10^{-6}$ — a three-way agreement.
The comparison is reproducible via `paper/validation/s4_crosscheck.py` and shown in
Fig. 1. [TODO: extend the table with a metal grating and a 2D meta-atom case.]

![CELERIS vs Stanford S$^4$ on an identical freestanding grating. Left: 0th-order
transmittance vs Fourier order $M$ — CELERIS (TE) and S$^4$ (TE) sit on grcwa's 0.93333.
Right: the CELERIS$-$S$^4$ difference falls to $7\times10^{-8}$ (TE); TM converges more
slowly for both, as expected.](figures/fig_validation.pdf)

The end-to-end pipeline is exercised in Figs. 2–3: sweeping a TiO$_2$ square-pillar's
fill builds the phase library (Fig. 2), and designing an NA$\approx$0.2 focusing lens
from it yields a diffraction-limited focal spot — FWHM equal to $\lambda f/D$ to three
figures (Fig. 3).

![Meta-atom phase library: transmission phase and amplitude $|t|$ vs pillar fill for a
TiO$_2$ square pillar ($\lambda$=532 nm, $\Lambda$=350 nm, $H$=600 nm).](figures/fig_library.pdf)

![Focal-plane PSF of the designed lens (log scale, left) and its central line cut
(right): FWHM = 1.330 $\mu$m, exactly the diffraction limit $\lambda f/D$.](figures/fig_psf.pdf)

The **full pipeline** is demonstrated by reproducing two canonical published devices
end to end:

- **Khorasaninejad et al., *Science* 352, 1190 (2016)** — an NA = 0.80 Pancharatnam–
  Berry TiO$_2$ nanofin metalens. CELERIS solves the published nanofin geometry with the
  deposited ALD-TiO$_2$ optical constants and recovers near-ideal half-wave-plate
  behavior (transmitted-normalized conversion 96–99 %), an exact geometric phase, and a
  diffraction-limited focal spot. [TODO: figure — reproduced efficiency + PSF.]
- **Chen et al., *Nat. Nanotechnol.* 13, 220 (2018)** — a broadband achromatic visible
  metalens. CELERIS reproduces the central physical limit (the required group-delay span
  sits at the ~5 fs single-nanofin ceiling, so the aperture is group-delay-limited) and,
  via single-etch group-delay engineering on a birefringent Pancharatnam–Berry library,
  flattens the chromatic focal shift relative to a dispersion-blind lens built from the
  same library (Fig. 4).

![Broadband achromatic metalens (Chen 2018 recipe): rigorous focal length vs wavelength
for a standard Pancharatnam–Berry lens and a group-delay-engineered achromatic lens
built from the *same* single-etch birefringent library. Group-delay matching flattens
the chromatic focal shift (drift 4.6 → 1.9 $\mu$m). At this small demonstration aperture
($D=10\,\mu$m, low Fresnel number) the *absolute* focus of both designs sits below the
geometric target $f$; the effect shown is the *relative* achromatic flattening, which is
what the group-delay objective controls.](figures/fig_achromatic.pdf)

[TODO: research-impact — add any external adopters / labmate use before submission;
SciPost requires specific, non-aspirational evidence of use.]

## Performance

The GPU accelerates far-field propagation only; the RCWA modal solve runs on the CPU.
We benchmark the propagation kernel (single precision) against CELERIS's own optimized
16-core CPU propagation (`compute_psf`, double precision) using the built-in `psfbench`
command on an RTX 4070 laptop GPU. Results agree to machine precision
(max$|\Delta$ normalized PSF$| = 0$). Measured speedups (GPU vs 16-core CPU):

| Aperture (pillars) | Focal grid | CPU | GPU | Speedup |
|---|---|---|---|---|
| 92 k | 161$^2$ | 178–202 ms | 35–37 ms | **4.8–5.8×** |
| 92 k | 321$^2$ | 274 ms | 140 ms | 2.0× |
| 92 k | 641$^2$ | 666 ms | 519 ms | 1.3× |
| 92 k | 1024$^2$ | 1457 ms | 1312 ms | 1.1× |
| 256 k | 161$^2$ | 241 ms | 111 ms | 2.2× |
| 577 k | 161$^2$ | 434 ms | 252 ms | 1.7× |

So on this hardware the GPU gives a **modest ~1–6× over a well-parallelized 16-core CPU**,
not the order-of-magnitude win a naive single-threaded baseline would suggest. The
speedup *decreasing* with grid size indicates the current kernel is memory-bound (each
thread re-reads the full pillar list from global memory); shared-memory tiling is a clear
future optimization. We report these numbers rather than a headline figure precisely
because an apples-to-apples multicore comparison is the honest one.
[TODO: optionally add single-core CPU and larger-GPU data points; regenerate via
`reproduce_all`.]

The dominant cost, however, is the RCWA modal solve, not propagation. The solve is
bounded by dense complex linear algebra (operator assembly, matrix inverses, and the
S-matrix recursion), which an optional Intel MKL backend accelerates. Because CELERIS
parallelizes its unit-cell sweeps across cores itself, MKL is pinned to a single thread
to avoid nested oversubscription; composed this way it is **2.3× faster than the AVX2
build** on an M = 8 metalens library design (5.9 s vs 13.5 s, 16-core), with larger gains
at higher Fourier order. The default build has no MKL dependency (header-only Eigen with
AVX2); MKL is opt-in.

## Limitations

State plainly (maturity, not weakness):

- **The RCWA eigensolve is CPU-only.** The GPU accelerates far-field propagation, not
  the modal solve. A GPU eigensolve was prototyped (cuSOLVER) and found no faster than
  the CPU path; the promising direction is batching the many independent unit-cell solves
  of a library, which is future work.
- **Inverse design uses finite-difference optimization** over low-dimensional parametric
  meta-atoms; an analytic adjoint (true gradient-based / differentiable RCWA) is planned.
- **Curved shapes (circle/ellipse) converge slowly** (staircasing of a curved interface),
  a known RCWA property confirmed against `grcwa`.
- **Achromatic multi-height designs imply a grayscale/multi-level etch**; the single-etch
  Pancharatnam–Berry variant is the directly fabricable one.
- **Scope is metalenses/metasurfaces via RCWA** — explicitly not FDTD or ray tracing.

## Use of AI assistance

CELERIS was developed with AI coding assistance (Anthropic Claude, via the Claude Code
CLI), used for code generation, refactoring, test scaffolding, build/portability work,
and drafting documentation. All AI-assisted output was reviewed, tested, and validated
by the author, who made all core design and physics decisions; correctness is enforced
by the validation suite described above. [TODO: state tool versions and date range at
submission.]

## References

[TODO: populate paper.bib — Moharam & Gaylord (RCWA); Li (Fourier factorization /
inverse rule); Liu & Fan (S$^4$/FMM); Redheffer (S-matrix); Rumpf (implementation);
Khorasaninejad 2016; Chen 2018; grcwa; TORCWA; fmmax; meent; RETICOLO. Cite the software
archive DOI (Zenodo) once minted.]
