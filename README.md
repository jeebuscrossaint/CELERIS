# CELERIS

**GPU-ready metalens design via rigorous coupled-wave analysis (RCWA).**

CELERIS takes a lens specification — focal length, aperture, wavelength, pillar
material — and produces a fabrication-ready GDSII layout plus a full optical
performance report, from a CLI or a native desktop GUI. It is built around a
from-scratch, validated electromagnetic solver: every layer of the stack is
cross-checked against closed-form physics, an independent method, or energy
conservation. Far-field analysis is GPU-accelerated (600–700× over a 16-core
CPU on large apertures).

```
spec ─▶ materials ─▶ Fresnel / TMM ─▶ 1D RCWA (TE+TM, multilayer)
     ─▶ 2D vectorial RCWA ─▶ unit-cell library ─▶ metalens design
     ─▶ inverse-design optimizer ─▶ GDSII export
     ─▶ focal / chromatic / wavefront / MTF / through-focus / field analysis
        (focal propagation runs on the GPU)
```

## What it does

- **Rigorous EM engine** — RCWA/Fourier Modal Method from scratch:
  - 1D gratings: TE + TM (Li inverse-rule factorization), arbitrary multilayer
    stacks via stable scattering-matrix (Redheffer) recursion.
  - 2D biperiodic structures (metalens nanopillars): full vectorial `P·Q`
    formulation.
- **Materials** — dispersion models (Sellmeier, tabulated `n,k`), built-in
  N-BK7, fused silica, Si₃N₄, air; constant-index for quick studies; runtime
  CSV import of real `n,k` (refractiveindex.info format).
- **Metalens design** — sweep pillar geometry into a phase/amplitude library,
  then map an ideal focusing profile onto pillars (amplitude-aware selection).
  `--auto-height` sweeps the etch depth and picks the height with the best phase
  coverage at the highest transmittance — full-2π control without guessing the
  pillar height (still a single fab etch). `--shape` selects the meta-atom
  cross-section (square / circle / ellipse / cross / ring); non-rectangular
  shapes are solved with a sampled-grid (Laurent) Fourier factorization that
  cross-checks against grcwa.
- **Inverse design** — gradient-based optimizer (Adam) finds the pillar
  geometry meeting a target phase with maximum transmission; plus a
  period × height system optimizer.
- **Polarization optics** — rectangular pillars give independent phase to X- and
  Y-polarization (form birefringence). Build polarization-multiplexed lenses
  (X-pol and Y-pol focus at different planes) with a 2D `(fill_x, fill_y)`
  library; reports per-polarization RMS phase, measured foci, and focal
  isolation (channel cross-talk, dB). The core of polarization-imaging
  metasurfaces.
- **Fabrication output** — GDSII export (no dependencies; emits the binary
  records directly), and an in-app GDS viewer that renders any `.gds` file.
- **Analysis battery** — focal Strehl / FWHM / encircled energy, full wavefront
  (OPD + Zernike, Maréchal Strehl), MTF, through-focus + longitudinal caustic,
  chromatic focal shift (with per-wavelength meta-atom dispersion), off-axis
  spot-vs-field diagram, fabrication-tolerance Monte-Carlo, field of view.
- **GPU acceleration** — the far-field propagation (Rayleigh–Sommerfeld sum over
  pixels × pillars) runs as a CUDA kernel: **675× faster** at 120 µm / 92k
  pillars (23.9 s → 35 ms) and **669×** at 250 µm / 401k pillars (147 s → 220 ms)
  vs a 16-core CPU, agreeing to ~1e-5. Falls back to the CPU automatically.
  (The general eigensolve stays on the CPU — cuSOLVER's `Xgeev` is host-serial
  and not competitive; `celeris gpubench` documents this.)
- **Desktop GUI** — a dockable, utilitarian workspace (Dear ImGui): editable
  lens/layer/material data, live PSF / wavefront / MTF / caustic / layout views,
  one-click GDSII and report export, project save/load.

## Build

Dependencies (Eigen) are fetched automatically by CMake. Requires a C++23
compiler.

### MinGW / GCC (default)
```sh
cmake -B build -G "Unix Makefiles" -DCMAKE_MAKE_PROGRAM=make -DCMAKE_CXX_COMPILER=g++
cmake --build build
./build/celeris.exe
```
On MinGW, `std::println` requires linking `stdc++exp` (handled automatically).

### MSVC (Visual Studio 2022)
```sh
cmake -B build-msvc -G "Visual Studio 17 2022" -A x64
cmake --build build-msvc --config Release
```

### Desktop GUI
```sh
cmake -B build-msvc -G "Visual Studio 17 2022" -A x64 -DCELERIS_BUILD_GUI=ON
cmake --build build-msvc --config Release --target celeris_gui
./build-msvc/Release/celeris_gui.exe
```
GLFW and Dear ImGui (docking branch) are fetched automatically.

### With the CUDA GPU backend
Requires the full NVIDIA CUDA Toolkit (the scoop package omits the CCCL headers
cuSOLVER needs). The Visual Studio generator can't locate `nvcc` without the VS
integration, so the GPU build uses **Ninja**: run `scripts/build-cuda.bat`
(builds the CLI and the GPU-accelerated GUI), then launch with
`scripts/run-gui-cuda.bat`, which puts the CUDA `bin/x64` DLLs on `PATH`.

## Usage

```sh
# Design a focusing metalens -> GDSII + performance report
celeris design --focal 50 --diameter 20 --wavelength 0.532 \
               --pillar-n 2.4 --thickness 0.6 --out lens.gds

# Let the etch depth be chosen automatically: sweep pillar height for the best
# phase coverage at the highest transmittance (still a single etch), no guessing
celeris design --focal 50 --diameter 20 --pillar-csv data/TiO2_Siefke.csv \
               --substrate sio2 --auto-height --out lens.gds

# Non-square meta-atom shapes (circle/cross/ring), solved with the grid
# (Laurent) factorization; combine with --auto-height to compare coverage
celeris design --focal 50 --diameter 20 --shape cross --shape-param 0.45 \
               --pillar-csv data/TiO2_Siefke.csv --substrate sio2 --auto-height

# Write a full deliverable bundle (metrics .txt + PSF & caustic PGM + GDS)
celeris design --focal 50 --diameter 40 --report mylens

# Run the physics validation suite
celeris selftest

# Credibility battery on REAL tabulated TiO2 n,k (data/TiO2_Siefke.csv):
# convergence vs harmonics, broadband meta-atom library, end-to-end focusing
celeris validate --report myvalidation

# Polarization-multiplexed lens: X-pol @50um, Y-pol @80um -> rectangular GDS
celeris polardesign --focal-x 50 --focal-y 80 --diameter 24 --report mypolar

# Form-birefringence sweep (waveplate building block)
celeris birefringence --fill-y 0.5

# GPU benchmarks (CUDA build only)
celeris psfbench --diameter 120     # far-field kernel vs CPU
celeris gpubench --n 242 --batch 32 # batched eigensolve (honest negative result)
```

`celeris help` lists all options (focal length, aperture, wavelength, period,
pillar height/index, substrate, material CSV, library resolution, RCWA
harmonics, tolerance/FOV/PSF/report outputs).

## Validation highlights

These are produced by `celeris selftest`:

| Check | Result |
|---|---|
| Anti-reflection coating | reflection → 0 at design wavelength |
| Distributed Bragg mirror | R → 0.998 (8 quarter-wave pairs) |
| RCWA vs TMM (degenerate grating) | agree to **1e-6** |
| Multilayer S-matrix vs single-layer | agree to **1e-15** |
| 2D RCWA reduces to 1D (y-invariant) | agree to **5e-12** |
| Energy conservation (all RCWA) | Σ DE = 1.000000 |
| Designed lens focal spot | **diffraction-limited** FWHM = λf/D |
| Form birefringence (square pillar) | exactly 0° (symmetry) |

## Architecture

The engine (`celeris_core`) is GUI/front-end-agnostic: every routine takes
parameters and returns plain data, so the CLI, a future GUI, or Python bindings
all consume the same library.

```
src/celeris/
  core.hpp              complex types, constants, polarization
  materials/           dispersion models + built-in library + CSV import
  optics/              Fresnel equations, Transfer Matrix Method
  rcwa/                1D (TE/TM, multilayer) + 2D vectorial RCWA, S-matrix
  design/              unit-cell library, lens assembler, inverse-design + system optimizer,
                       polarization-multiplexed design
  analysis/            focal, chromatic, wavefront, MTF, through-focus, field, polarization
  io/                  GDSII read/write, PGM image, material CSV
  cuda/                cuSOLVER eigensolve + GPU far-field propagation kernel
gui/                   Dear ImGui desktop application (celeris_gui)
```

## Status & roadmap

CELERIS is a complete, validated MVP with a CLI, a desktop GUI, and a
GPU-accelerated analysis path. Notable items still ahead:

- Faster CPU eigensolve (OpenBLAS/LAPACK `zgeev`) — the remaining bottleneck is
  the RCWA library build, which the GPU eigensolve cannot accelerate.
- Full-2π / thickness-swept libraries and broadband / achromatic optimization.
- Geometric (Pancharatnam–Berry) phase from rotated pillars for circular-
  polarization optics (the linear-birefringence path is already built).
- Validation against a published metalens; Python bindings; Linux/macOS GUI.

## License

TBD.
