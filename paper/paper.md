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
    ror: 036nfer12
date: 3 November 2026
bibliography: paper.bib
---

<!-- JOSS paper (target: Journal of Open Source Software). Body must stay 750-1750 words.
     The long-form methods/validation write-up (former CPC manuscript) lives in
     paper/longform/ and is intended as an arXiv preprint, not a second journal paper.
     API and usage documentation belong in README / docs/, not here (JOSS rule). -->

# Summary

A metalens is a flat lens made of millions of tiny pillars patterned on a chip.
Each pillar is smaller than the wavelength of light and delays the light that
passes through it by an amount set by its shape. Choosing the right pillar at
every point of the chip makes the light leaving it converge to a focus, as it
would behind a curved glass lens.

Designing such a lens takes four steps: (i) simulate one pillar with an exact
electromagnetic solver, for many pillar shapes, to build a library of responses;
(ii) pick the pillar from that library that gives the needed delay at each
point; (iii) predict how well the finished lens focuses; and (iv) write a mask
layout that a nanofabrication facility can use to make it.

CELERIS does all four steps in one open-source package. Its solver, written from
scratch in C++, uses rigorous coupled-wave analysis (RCWA), also called the
Fourier modal method [@moharam1981; @moharam1995]. Around it, CELERIS builds
pillar libraries, designs focusing, geometric-phase, and broadband (achromatic)
lenses, computes the standard measures of lens quality, and exports the layout
in the GDSII format used by fabrication tools. It can be used from a
command-line program, from Python, or from a desktop application. A test suite
checks its physics against known results on every change.

# Statement of need

Open RCWA codes such as S$^4$ [@liu2012s4], RETICOLO [@reticolo], `grcwa`
[@grcwa], TORCWA [@torcwa], `fmmax` [@fmmax], and `meent` [@meent] are solver
kernels: given one periodic structure, they return how it diffracts light. A
metalens researcher who starts from one of them must still write the rest of the
workflow: library sweeps, phase mapping, focal-plane propagation, imaging
metrics, tolerance studies, and layout export. That glue code is rarely tested
or shared, and an error in it is easy to miss because the solver itself is
correct.

CELERIS is for metasurface researchers and students who want to go from an
optical specification (focal length, aperture, wavelength band, material) to an
analyzed, fabrication-ready design with one tested tool. Because every stage
uses the same data model, CELERIS can also answer device-level questions that a
solver called in a loop does not expose directly. One example is the limit that
group delay sets on single-etch broadband achromats [@chen2018; @wang2018]: how
much aperture and bandwidth a given pillar library can support, which depends
jointly on the library, the design map, and the focal analysis.

# State of the field

The solver kernels above are mature, and several of them (TORCWA, `fmmax`,
`meent`) are stronger than CELERIS as kernels: they run the full modal solve on a
GPU and provide automatic differentiation. CELERIS solves on the CPU,
accelerates only far-field propagation on the GPU, and optimizes by finite
differences. Two open projects also go beyond a kernel. `metabox` [@metabox]
provides RCWA-characterized libraries, differentiable end-to-end optimization
(including broadband objectives), and GDSII export. `MetaCraft` [@metacraft]
provides database-driven metalens design and optimization. Both are
inverse-design frameworks built around optimizing a merit function.

CELERIS was built as a new package, not as a contribution to one of these, for
two reasons. First, its focus is *forward* design and diagnosis: a broad set of
optical analyses, closed-form geometric-phase and group-delay design recipes, and
reproductions of published devices on one validated model. That does not fit
the architecture of a differentiable-optimization framework, where every stage
must be written in an autodiff library. Second, it is a native C++ library with
thin Python and GUI layers, so a complete design runs in seconds on a laptop CPU
with no GPU or autodiff stack. The approaches are complementary: a lens found by
an inverse-design tool can be checked with CELERIS's analyses, and an analytic
adjoint for CELERIS is planned work.

# Software design

**One library, three front ends.** All physics lives in one C++23 library,
`celeris_core`. Every routine takes plain parameters and returns plain data
structures. The command-line program, the pybind11
[@pybind11] Python module, and the Dear ImGui desktop application are thin
layers that call the same functions. A design made in the GUI and one scripted
in Python therefore cannot disagree, and testing the library tests the physics
behind all three.

**CPU first, accelerators optional.** The default build needs only a C++23
compiler and CMake; Eigen [@eigen] is the linear-algebra backend. CUDA, Intel
MKL, and the GUI are detected at configure time and fall back cleanly when
absent, so a reviewer on a Linux laptop without a GPU builds and tests the same
physics. The GPU is used only where profiling showed a large gain. The
far-field propagation sum (focal-plane pixels × pillars) runs as a CUDA kernel,
about 120–190× faster than the multithreaded CPU path on the author's laptop. A
GPU eigensolve was tried and was not competitive, because the general complex
eigensolver in cuSOLVER runs mostly on the host. That result is kept in the
repository as a benchmark rather than hidden. For the RCWA solve itself, the
costly part is dense matrix algebra (operator assembly and the
scattering-matrix recursion), which MKL speeds up on the CPU.

**Accuracy where the cost is.** The 2D solver uses Li's Fourier factorization
rules [@li1996; @li1997] for axis-aligned rectangles, the dominant metalens
pillar, because they converge much faster in Fourier order than the naive
(Laurent) rule. Arbitrary sampled shapes use the Laurent rule, which is simpler
and general but converges more slowly; the documentation states this. Layers are
combined with a Redheffer scattering-matrix recursion [@redheffer1962;
@rumpf2011], which stays stable for thick and lossy layers. Far fields use a
direct Rayleigh–Sommerfeld sum [@goodman] rather than an FFT, so the focal plane
can be sampled at any resolution and offset (for example, centered on an
off-axis chief ray) without padding.

**Design recipes, not only optimizers.** Propagation-phase lenses choose a
pillar per site from a library. Geometric-phase (Pancharatnam–Berry) lenses
[@berry1984; @bomzon2002] set the phase exactly by rotating one birefringent
pillar. Achromatic lenses match both phase and group delay across a band. These
closed-form recipes give a fast, explainable baseline. A finite-difference
optimizer covers single-pillar inverse design; an analytic adjoint is planned.

**Testing as part of the design.** `celeris selftest` runs 56 assertions with
locked tolerances and exits non-zero if any fails. They check physical
invariants and external references: energy conservation; agreement of the 2D
solver with the 1D solver for TE and TM on a grating that does not vary along
one axis (to about $10^{-14}$); RCWA against the transfer-matrix method and
effective-medium theory; and agreement with `grcwa` (TE transmittance 0.93333,
matched to the five digits quoted). On an identical binary grating, CELERIS
matches Stanford S$^4$ to $7\times10^{-8}$ in TE transmittance at high Fourier
order (\autoref{fig:validation}). Continuous integration builds CELERIS on Linux
with CPU only, runs a fast subset of the suite on every push and pull request,
and runs the full suite nightly. `paper/reproduce_all.sh` regenerates every
number and figure in this paper.

![CELERIS versus Stanford S$^4$ on an identical freestanding $n=1.5$ grating
($\Lambda=0.3\,\mu$m, fill 0.5, thickness $0.5\,\mu$m, $\lambda=0.5\,\mu$m).
Left: zeroth-order transmittance against Fourier order. Right: the CELERIS–S$^4$
TE difference falls to $7\times10^{-8}$. In TM, CELERIS (Li's factorization) is within
$3\times10^{-5}$ of its $M=40$ value by $M=6$, while S$^4$ with default settings
approaches the same value slowly ($6\times10^{-4}$ apart at $M=40$).\label{fig:validation}](figures/fig_validation.png)

# Research impact statement

<!-- AUTHOR: JOSS requires this to be "compelling and specific, not aspirational".
     Before submitting, add concrete use with links (the Zenodo group-delay study,
     ROADMAP R1; any poster/talk/preprint). Log uses in memory/research-use.md only
     from what the author reports. What is below is the honest current state. -->

CELERIS reproduces two published devices end to end; these are its worked
research examples. For the NA = 0.80 Pancharatnam–Berry TiO$_2$ metalenses of
@khorasaninejad2016 (405, 532, and 660 nm), it solves the published nanofin
geometries with measured optical constants of deposited TiO$_2$. The nanofins act
as near-ideal half-wave plates: 96–99.8 % of the transmitted light is converted
to the opposite circular polarization. Lenses built from them, 20 µm across, reach a
diffraction-limited focal spot. For the broadband achromat of @chen2018
(NA = 0.20, $D=26.4\,\mu$m, one 600 nm etch, 470–670 nm), group-delay
engineering with a single-etch library cuts the chromatic focal drift from
22.7 µm to 7.3 µm compared with a dispersion-blind lens from the same library
(\autoref{fig:achromatic}). The cost is mean conversion efficiency, which falls
from 0.52 to 0.31. The remaining drift is the group-delay limit: the aperture
needs a span of about 4.4 fs, while a single 600 nm nanofin library spans about
2.6 fs. Version 1.0.0 is archived on Zenodo [@celeris_v1].

![Chen-2018 achromat at full device scale: rigorous focal length against
wavelength for a standard Pancharatnam–Berry lens and a group-delay-engineered
lens built from the same single-etch library.\label{fig:achromatic}](figures/fig_achromatic.png)

# AI usage disclosure

<!-- AUTHOR: confirm the exact model versions used over the project's lifetime. -->

CELERIS was developed with substantial help from generative AI. The author used
Anthropic's Claude Code (Claude Opus-family models, 2026) to help write and
refactor source code, tests, build scripts, documentation, and drafts of this
paper. The author chose the scope, architecture, numerical methods, and
validation targets, and reviewed, edited, and tested all AI-assisted output.
Correctness is checked in three ways. The locked self-test suite runs in
continuous integration. Every reported number is regenerated by the scripts in
`paper/`. Before submission, an audit re-derived each number claimed in the
documentation and this paper. That audit found, and public issues and pull
requests fixed, a TM Fourier-factorization error and a GPU benchmark that had
timed the GPU against itself. A new regression test guards the TM fix. No AI tools
were used for correspondence with editors or reviewers.

# Acknowledgements

This work received no dedicated financial support.

# References
