# CELERIS

**GPU-ready metalens design via rigorous coupled-wave analysis (RCWA).**

CELERIS takes a lens specification — focal length, aperture, wavelength, pillar
material — and produces a fabrication-ready GDSII layout plus a full optical
performance report, from a CLI or a native desktop GUI. It is built around a
from-scratch, validated electromagnetic solver: every layer of the stack is
cross-checked against closed-form physics, an independent method, or energy
conservation. Far-field propagation can optionally run on a CUDA GPU (measured
~1–6× over a 16-core CPU on an RTX 4070, depending on problem size; the RCWA
solve itself is CPU).

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

Dependencies (Eigen, and for the GUI: GLFW + Dear ImGui) are fetched
automatically by CMake. Requires a C++23 compiler.

```sh
cmake -B build           # auto-detects this machine and builds the best variant
cmake --build build --config Release
./build/Release/celeris.exe selftest
```

**The configure step auto-detects the machine** (`CELERIS_AUTO=ON` by default)
and enables the fastest backend it can, degrading gracefully — so one build dir
covers every machine instead of one per backend. It prints a summary like:

```
================ CELERIS build configuration ================
  Vectorization : AVX2/FMA
  GPU backend   : CUDA cuSOLVER + far-field kernel (nvcc OK)
  Eigensolve+BLAS: Intel MKL (multithreaded)
  Desktop GUI   : ON (celeris_gui)
  Runtime DLLs  : bundled next to the exe (runs with no PATH setup)
=============================================================
```

What it detects, and the graceful fallback when a piece is missing:

| Capability | Detected from | If absent |
|---|---|---|
| **AVX2/FMA** | on by default (`-DCELERIS_AVX2=OFF` for pre-2013 CPUs) | SSE2 baseline |
| **CUDA** GPU (far-field kernel, ~1–6× vs 16-core CPU; static runtime, no DLLs) | `CUDA_PATH`/`CUDA_HOME`/`/opt/cuda` toolkit + `nvcc` | CPU propagation |
| **Intel MKL** (multithreaded BLAS + eigensolve) | `$MKLROOT`, conda, oneAPI, or a `pip install mkl-devel` prefix | header-only Eigen (AVX2) |

Override any of these explicitly: `-DCELERIS_USE_CUDA=ON|OFF|AUTO`,
`-DCELERIS_USE_MKL=ON|OFF|AUTO`. To have the build *fetch* MKL when it's missing,
add `-DCELERIS_FETCH_MKL=ON` (runs `pip install mkl-devel mkl-include`). The
desktop GUI (`celeris_gui`) builds by default (it fetches GLFW + Dear ImGui);
turn it off for CLI-only / headless / CI builds with `-DCELERIS_BUILD_GUI=OFF`.

**Speed**, in order: AVX2 (free, default) is ~2.7× over the SSE2 baseline; adding
MKL reaches up to ~11× at high harmonic counts. The hot path is dense complex
linear algebra (operator assembly + S-matrix recursion), *not* the eigensolve.

**Deployment**: when MKL is used, the build bundles its runtime DLLs (~470 MB)
next to the executables, so the build folder is self-contained — `celeris.exe`
and `celeris_gui.exe` run with no PATH setup and no launcher (Windows resolves
DLLs from the exe's own directory). Ship the `Release/` folder as one unit. The
GPU far-field kernel links the **static** CUDA runtime, so it needs no CUDA DLLs.
(The cuSOLVER GPU eigensolve is a separate opt-in, `-DCELERIS_USE_CUDA_EIGENSOLVE=ON`
— it's a measured dead-end and would drag in ~630 MB of cuSOLVER/cuBLAS DLLs.)

### Notes per toolchain
- **MinGW/GCC**: `cmake -B build -G "Unix Makefiles" -DCMAKE_CXX_COMPILER=g++`;
  `std::println` links `stdc++exp` automatically.
- **MSVC**: the default generator works; `cmake -B build` then
  `cmake --build build --config Release`.
- **CUDA**: `cmake -B build` builds the GPU kernel automatically — it compiles
  `propagate.cu` by invoking the detected `nvcc` directly (a custom command), so
  it does **not** depend on CMake's CUDA language, the VS generator's MSBuild CUDA
  integration, the compiler-id probe, or a `vcpkg integrate` toolchain — all of
  which historically broke CUDA under the Visual Studio generator. The runtime is
  linked **statically** (`cudart_static`), so the GPU build adds **zero** CUDA
  DLLs. Works the same under the VS generator, Ninja, or Makefiles. (nvcc still
  needs MSVC `cl` as its host compiler — that's an NVIDIA requirement, not ours;
  the VS generator finds `cl` automatically, which is why it's the easy default.)
  `scripts/build-cuda.bat` (Ninja) remains as an alternative.

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

# Pancharatnam-Berry (geometric-phase) lens: one half-wave-plate atom rotated
# per site imprints 2*theta on circularly polarized light. RCWA-verifies the
# 2*theta relation; writes a rotated-pillar GDS.
celeris pbdesign --focal 50 --diameter 20 --report mypb

# Form-birefringence sweep (waveplate building block)
celeris birefringence --fill-y 0.5

# GPU benchmarks (CUDA build only)
celeris psfbench --diameter 120     # far-field kernel vs CPU
celeris gpubench --n 242 --batch 32 # batched eigensolve (honest negative result)
```

`celeris help` lists all options (focal length, aperture, wavelength, period,
pillar height/index, substrate, material CSV, library resolution, RCWA
harmonics, tolerance/FOV/PSF/report outputs).

## Python

The same validated engine is scriptable from Python (numpy/matplotlib), the way
the metasurface community actually works. Build the extension once, then `import
celeris`:

```sh
cmake -B build -DCELERIS_BUILD_PYTHON=ON
cmake --build build --config Release --target _celeris   # -> build/python/celeris
# or:  pip install .        (scikit-build-core drives the same CMake build)
```

```python
import celeris as cel

tio2 = cel.Material.constant(2.40 + 0j, "TiO2")           # or Material.tabulated(...)
lib  = cel.build_unit_cell_library(                       # phase library (RCWA sweep)
    pillar=tio2, background=cel.materials.air(), incident=cel.materials.air(),
    substrate=cel.materials.fused_silica(), period_um=0.35, wavelength_um=0.532,
    thickness_um=0.6, fill_min=0.1, fill_max=0.9, n_samples=32, M=8)
lens = cel.design_metalens(lib, focal_length_um=50.0, diameter_um=20.0)
foc  = cel.analyze_focus(lens, lib, 50.0, 0.532, 20.0)
print(foc.strehl, foc.fwhm_um)                            # diffraction-limited
psf  = cel.compute_psf(lens, lib, 50.0, 0.532, 20.0, 128, 4.0)  # numpy maps
```

Runnable examples (focusing lens, RCWA convergence, real-material dispersion,
full-2π height sweep) are in [`examples/python/`](examples/python/).

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
                       polarization-multiplexed + Pancharatnam-Berry (geometric-phase) design
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
- Arbitrary geometric-phase profiles (vortex/OAM, deflectors) on the
  Pancharatnam–Berry path — focusing is already built (`pbdesign`).
- Validation against a published metalens; Linux/macOS GUI.

## License

TBD.
