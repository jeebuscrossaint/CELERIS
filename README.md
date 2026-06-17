# CELERIS

**GPU-ready metalens design via rigorous coupled-wave analysis (RCWA).**

CELERIS takes a lens specification — focal length, aperture, wavelength, pillar
material — and produces a fabrication-ready GDSII layout plus a full optical
performance report. It is built around a from-scratch, validated electromagnetic
solver: every layer of the stack is cross-checked against closed-form physics,
an independent method, or energy conservation.

```
spec ─▶ materials ─▶ Fresnel / TMM ─▶ 1D RCWA (TE+TM, multilayer)
     ─▶ 2D vectorial RCWA ─▶ unit-cell library ─▶ metalens design
     ─▶ inverse-design optimizer ─▶ GDSII export ─▶ focal & chromatic analysis
```

## What it does

- **Rigorous EM engine** — RCWA/Fourier Modal Method from scratch:
  - 1D gratings: TE + TM (Li inverse-rule factorization), arbitrary multilayer
    stacks via stable scattering-matrix (Redheffer) recursion.
  - 2D biperiodic structures (metalens nanopillars): full vectorial `P·Q`
    formulation.
- **Materials** — dispersion models (Sellmeier, tabulated `n,k`), built-in
  N-BK7, fused silica, air; constant-index for quick studies.
- **Metalens design** — sweep pillar geometry into a phase/amplitude library,
  then map an ideal focusing profile onto pillars (amplitude-aware selection).
- **Inverse design** — gradient-based optimizer (Adam) finds the pillar
  geometry meeting a target phase with maximum transmission.
- **Fabrication output** — GDSII export (no dependencies; emits the binary
  records directly).
- **Analysis** — focal-plane metrics (Strehl ratio, spot FWHM vs the
  diffraction limit, encircled energy) and chromatic focal shift across a band.
- **Optional GPU backend** — cuSOLVER eigensolve (opt-in; see notes).

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

### With the CUDA GPU backend
Requires the full NVIDIA CUDA Toolkit (the scoop package omits the CCCL headers
cuSOLVER needs). See `scripts/build-cuda.bat`, then run with the CUDA
`bin/x64` directory on `PATH`.

## Usage

```sh
# Design a focusing metalens -> GDSII + performance report
celeris design --focal 50 --diameter 20 --wavelength 0.532 \
               --pillar-n 2.4 --thickness 0.6 --out lens.gds

# Run the physics validation suite
celeris selftest
```

`celeris help` lists all options (focal length, aperture, wavelength, period,
pillar height/index, substrate, library resolution, RCWA harmonics, output).

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
  materials/           dispersion models + built-in library
  optics/              Fresnel equations, Transfer Matrix Method
  rcwa/                1D (TE/TM, multilayer) + 2D vectorial RCWA, S-matrix
  design/              unit-cell library, lens assembler, inverse-design optimizer
  analysis/            focal-plane and chromatic analysis
  io/                  GDSII export
  cuda/                optional cuSOLVER eigensolve backend
```

## Status & roadmap

CELERIS is a complete, validated MVP. Notable items still ahead:

- Faster CPU eigensolve (OpenBLAS/LAPACK `zgeev`) and batched GPU sweeps.
- Full-2π / thickness-swept libraries and broadband / achromatic optimization.
- Polarization-multiplexed design (the vectorial solver already supports it).
- Python bindings and a native cross-platform GUI.

## License

TBD.
