---
title: "CELERIS: an integrated, validated RCWA pipeline for metalens design, analysis, and fabrication layout"
tags:
  - C++
  - Python
  - optics
  - metasurfaces
  - metalenses
  - rigorous coupled-wave analysis
  - Fourier modal method
  - GDSII
authors:
  - name: Amarnath Patel
    orcid: 0009-0008-9460-082X
    corresponding: true
    affiliation: 1
affiliations:
  - name: University of Central Florida, Orlando, FL, United States
    index: 1
date: 19 October 2026
bibliography: paper.bib
---

<!-- JOSS paper (target: Journal of Open Source Software). Body must stay 750-1750 words.
     The long-form methods/validation write-up (former CPC manuscript) lives in
     paper/longform/ and is intended as an arXiv preprint, not a second journal paper. -->

# Summary

Metalenses are flat lenses made of millions of sub-wavelength pillars
("meta-atoms") patterned on a chip. Each pillar delays light by an amount set by
its shape, so a lens is designed by choosing the right pillar at every point of
the aperture. In practice this is a four-step workflow: (i) simulate one periodic
unit cell with a rigorous electromagnetic solver to build a library of pillar
responses, (ii) map a target phase profile onto that library, (iii) predict how
well the resulting device focuses, and (iv) export a mask layout that a
nanofabrication facility can write.

CELERIS implements this whole workflow in one package. Its core is a
from-scratch C++23 solver based on rigorous coupled-wave analysis (RCWA), also
called the Fourier modal method [@moharam1981; @moharam1995], in 1D and
2D-vectorial form. Around the solver CELERIS provides a meta-atom library
builder, propagation-phase, geometric-phase (Pancharatnam–Berry), and broadband
achromatic design recipes, a battery of optical analyses (point-spread function,
Strehl ratio, MTF, wavefront/Zernike, through-focus, chromatic and tolerance
studies), optional GPU-accelerated far-field propagation, and GDSII layout
export. It is driven from a command-line interface, a Python module, or a
desktop GUI, and its physics is locked by a self-test suite that runs in
continuous integration.

# Statement of need

Open RCWA codes such as S$^4$ [@liu2012s4], RETICOLO [@reticolo], `grcwa`
[@grcwa], TORCWA [@torcwa], `fmmax` [@fmmax], and `meent` [@meent] are solver
kernels: given a periodic structure, they return its diffraction efficiencies
and fields. A metalens researcher who starts from one of them still has to write
the surrounding workflow — library sweeps, phase mapping, focal-plane
propagation, imaging metrics, tolerance studies, and layout export — and that
glue code is rarely validated or shared.

CELERIS is aimed at metasurface researchers and students who want to go from an
optical specification (focal length, aperture, wavelength band, material) to an
analyzed, fabrication-ready design with a single, tested tool. Because all stages
share one data model, CELERIS can also answer device-level questions that a
solver called in a loop does not expose directly. The main example is the
*group-delay design limit* of single-etch broadband achromats [@chen2018;
@wang2018]: the aperture–bandwidth trade-off set jointly by the dispersive
library, the design map, and the focal analysis.

# State of the field

The solver kernels above are mature, and several of them (TORCWA, `fmmax`,
`meent`) are stronger than CELERIS as kernels: they run the full modal solve on a
GPU and provide automatic differentiation, whereas CELERIS solves on the CPU,
accelerates only far-field propagation on the GPU, and optimizes by finite
differences. Two open projects also go beyond a kernel: `metabox` [@metabox]
provides RCWA-characterized libraries, differentiable end-to-end optimization
(including broadband objectives), and GDSII export, and `MetaCraft`
[@metacraft] provides database-driven metalens design and optimization. Both
are inverse-design frameworks built around optimizing a merit function.

CELERIS was built as a new package, rather than as a contribution to one of
these, for two reasons. First, its focus is *forward* design and diagnosis — a
broad optical-analysis battery, closed-form geometric-phase and group-delay
recipes, and reproductions of published devices on one validated model — which
does not fit the architecture of a differentiable-optimization framework.
Second, it is written as a native C++ library with thin Python and GUI layers, so
a complete design runs in seconds on a laptop CPU without a GPU or an autodiff
stack. The two approaches are complementary: a lens found by an inverse-design
tool can be checked with CELERIS's analysis battery, and adding an analytic
adjoint to CELERIS is planned work.

# Software design

- **Engine (`celeris_core`, C++23).** 1D RCWA (TE/TM, multilayer) and
  2D-vectorial RCWA. The 2D solver uses Li's inverse-rule Fourier
  factorization [@li1996; @li1997] for axis-aligned rectangles and a Laurent
  factorization for arbitrary sampled shapes; layers are combined with a
  Redheffer scattering-matrix recursion [@redheffer1962; @rumpf2011]. Dense
  linear algebra uses Eigen [@eigen], with an optional Intel MKL backend.
- **Design.** Library sweeps over pillar geometry and height; phase profiles
  (focusing, vortex, deflector, axicon, quadratic wide-field, freeform) drive
  both a propagation-phase path and a Pancharatnam–Berry path [@berry1984;
  @bomzon2002]; achromatic design matches phase and group delay across a band.
- **Analysis.** Rayleigh–Sommerfeld far-field propagation [@goodman], with an
  optional CUDA kernel and automatic CPU fallback, feeding PSF, Strehl,
  encircled energy, MTF, Zernike wavefront, through-focus, field-of-view,
  tolerance, and per-order efficiency analyses.
- **Output and interfaces.** GDSII export (square, rectangular, and rotated
  pillars); a CLI (`design`, `pbdesign`, `achromatic`, `reproduce`, `selftest`,
  …); Python bindings via pybind11 [@pybind11] (`pip install .`); and a Dear
  ImGui desktop GUI. The default build needs only a C++23 compiler and CMake;
  CUDA, MKL, and the GUI are optional.

**Testing.** `celeris selftest` runs 55 locked-tolerance assertions and exits
non-zero if any fails. The checks include energy conservation, agreement of the
2D solver with the 1D solver (~$5\times10^{-12}$), RCWA against the transfer-matrix
method and effective-medium theory, and agreement with external solvers. On an
identical binary grating, CELERIS matches Stanford S$^4$ to $7\times10^{-8}$ in
TE transmittance at Fourier order $M=40$ (\autoref{fig:validation}), and it matches
`grcwa`'s reference value of 0.93333. Continuous integration builds CELERIS on
Linux with CPU only and runs a fast subset on every push, plus the full suite
nightly; `paper/reproduce_all.sh` regenerates every figure.

![CELERIS versus Stanford S$^4$ on an identical freestanding $n=1.5$ grating
($\Lambda=0.3\,\mu$m, fill 0.5, thickness $0.5\,\mu$m, $\lambda=0.5\,\mu$m).
Left: zeroth-order transmittance against Fourier order. Right: the CELERIS–S$^4$
difference falls to $7\times10^{-8}$ (TE); TM converges more slowly in both
solvers.\label{fig:validation}](figures/fig_validation.png)

# Research impact statement

<!-- AUTHOR: JOSS requires this to be "compelling and specific, not aspirational".
     Before submitting, add concrete external use (another group, a course, a lab
     project, a poster/talk, a preprint that uses CELERIS) with links. What is below is
     the honest current state. -->

CELERIS reproduces two published devices end to end, and these serve as its
worked research examples. For the NA = 0.80 Pancharatnam–Berry TiO$_2$ metalens
of @khorasaninejad2016, it solves the published nanofin geometry with deposited
TiO$_2$ optical constants and recovers near-ideal half-wave-plate behavior
(96–99 % transmission-normalized conversion) and a diffraction-limited focus.
For the broadband achromat of @chen2018 (NA = 0.20, $D=26.4\,\mu$m, one
600 nm etch, 470–670 nm), single-etch group-delay engineering cuts the chromatic
focal drift from 22.7 µm to 6.4 µm compared with a dispersion-blind lens built
from the same library (\autoref{fig:achromatic}). The remaining drift is the
group-delay limit: the aperture needs a span of about 4.4 fs, while a single
600 nm nanofin library spans about 3 fs. Version 1.0.0 is
archived on Zenodo (doi:10.5281/zenodo.21891346).

![Chen-2018 achromat at full device scale: rigorous focal length against
wavelength for a standard Pancharatnam–Berry lens and a group-delay-engineered
lens built from the same single-etch library.\label{fig:achromatic}](figures/fig_achromatic.png)

# AI usage disclosure

<!-- AUTHOR: confirm the exact model versions used over the project's lifetime. -->

CELERIS was developed with substantial assistance from generative AI. The author
used Anthropic's Claude Code (Claude Opus-family models, 2026) to help write and
refactor source code, tests, build scripts, documentation, and drafts of this
paper. The author chose the scope, architecture, numerical methods, and
validation targets, and reviewed, edited, and tested all AI-assisted output.
Every numerical result reported here is checked by the locked self-test suite or
regenerated by the scripts in `paper/`. No AI tools were used for correspondence
with editors or reviewers.

# Acknowledgements

This work received no dedicated financial support.

# References
