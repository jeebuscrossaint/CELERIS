# Changelog

All notable changes to CELERIS are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project aims to
follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed
- Retargeted the code paper to the **Journal of Open Source Software (JOSS)**.
  `paper/paper.md` is rewritten in JOSS format (Summary, Statement of need, State of
  the field, Software design, Research impact statement, AI usage disclosure; ~1,000
  words, PNG figures). The former CPC manuscript moved to `paper/longform/` as an arXiv
  companion report; the CPC cover letter and highlights were removed.
- `paper.bib`: journal names written in full, as JOSS requires.

### Added
- `.github/workflows/draft-pdf.yml` — builds the JOSS paper PDF with the Open Journals
  toolchain on every paper change.
- Issue templates (bug, feature, question) and a pull-request template with the
  selftest checklist; a "Getting help" section in `CONTRIBUTING.md`.
- README sections on testing, getting help, and citing; refreshed status/roadmap.
- Selftest [6](c): a y-invariant 2D cell must reproduce the 1D solver for both TE and
  TM to 1e-10 (high-contrast n=3.5 grating). Restores the 2D = 1D invariant, now
  asserted and covering TM.

### Fixed
- **TM Fourier factorization now follows Li's rules (#8).** The 1D solvers
  (`rcwa1d.cpp`, `rcwa_stack1d.cpp`) had ⟦ε⟧ and ⟦1/ε⟧ swapped in the TM operator and
  its companion field, and the 2D plain-rectangle path used ⟦1/ε⟧ instead of ⟦ε⟧⁻¹ for
  the E_z term (the general-shape path was already correct). TM results converged
  slowly or not at all for high-contrast gratings, and rectangular-pillar phases did
  not converge in M (the selftest [8](b) TiO₂ pillar was off by 15–26° at M=6–12). Results change for every TM / x-polarized
  rectangular-pillar solve; TE is unchanged.

## [1.0.0] - 2026-08-11

First tagged release — the version described by the CELERIS code paper (targeting
Computer Physics Communications) and archived for its Zenodo DOI.

### Added
- Apache-2.0 `LICENSE`, `CITATION.cff`, `CONTRIBUTING.md`, `CODE_OF_CONDUCT.md`, and
  this changelog — publication/open-source readiness (Computer Physics Communications, arXiv).
- `paper/paper.tex` — the manuscript in Elsevier `elsarticle` format (targets CPC), with
  the CPC Program Summary block; compiles clean (`pdflatex` + `bibtex`, 10 pp).
- `selftest --quick` — a fast subset (TMM/1D/Bragg/registry, ~1 s) for CI on every push;
  the full 22-case suite runs on a nightly schedule.
- GitHub Actions CI: Linux CPU-only build + `celeris selftest` + Python examples.

### Changed
- Retargeted the code paper from SciPost Physics Codebases to **Computer Physics
  Communications** (SciPost bars solo-undergraduate submission).
- Split the 3325-line `main.cpp` CLI monolith into per-command translation units under
  `cli/` (dispatcher `main.cpp` is now 38 lines); behavior identical.

### Performance
- Parallelized the physics self-test: **276 s → 88 s** (3.1×) on 16 cores, with output
  byte-identical to the serial suite. Case-level thread pool + `std::async` on the
  independent 2D-RCWA eigensolves of the dominant case (which was 85 % of the runtime).
- **Cross-platform CUDA build.** The GPU far-field kernel now builds on Linux, not just
  Windows: CMake detects `/opt/cuda` (or `CUDA_PATH`) via `lib64/libcudart_static.a`,
  drives `nvcc` with the correct host-compiler flags per platform (`-Xcompiler=-fPIC,…`
  on Linux vs `/MD,/EHsc,…` on Windows), and links the Linux static runtime
  (`libcudart_static.a` + `libculibos.a` + `dl`/`rt`/`pthread`).

### Fixed
- **CI build.** The Linux CI used `ubuntu-latest`'s default `g++-13`, which lacks the
  C++23 `<print>` header the CLI relies on; CI now installs and selects `g++-14`.
- **MKL thread oversubscription on Linux.** CELERIS already parallelizes its workloads
  across cores; MKL's own per-call threading nested on top and ran ~5× *slower* than the
  AVX2 build by default (M=8 design: 33 s nested vs 5.9 s pinned). `eig.cpp` now pins MKL
  to one thread by default (respecting an explicit `MKL_NUM_THREADS`), so MKL composes
  with the outer parallelism: **2.3× faster than the AVX2 build at M=8**, out of the box.
- Python bindings failed to build on Linux: `solve_rcwa_2d` had 10 `py::arg`
  annotations for 11 parameters (the `OrderEfficiency` out-param) — a hard error under
  the fetched pybind11. Bound via a lambda that exposes the clean 10-argument API.
- `celeris_core` was not position-independent, so it could not link into the `_celeris`
  shared module on Linux. Scoped `POSITION_INDEPENDENT_CODE` to the Python build.

### Verified
- Clean **Linux CPU-only build** (GCC 16, CMake 4.4, Ninja) with no CUDA/MKL/MSVC —
  auto-detected CPU-only backend, Eigen fetched at configure time.
- **Full `celeris selftest` passes on Linux** (exit 0, all 21 locked cases): energy
  conservation, RCWA=TMM/EMT, grcwa cross-check, and the Khorasaninejad 2016 / Chen 2018
  device reproductions.
- **Cross-checked against Stanford S$^4$** (the field-standard FMM reference), built on
  Linux and run head-to-head on an identical grating: 0th-order transmittance agrees to
  **7e-8 (TE)** / 3e-4 (TM, slower convergence in both), a three-way match with grcwa
  (0.93333). Reproducible via `paper/validation/s4_crosscheck.py`.
- **CUDA path runs on Linux hardware** — `cuda::available()` true,
  `device_name()` = "NVIDIA GeForce RTX 4070 Laptop GPU". (GPU accelerates far-field
  propagation only; the RCWA eigensolve remains CPU — stated precisely for the paper.)
- Python wheel builds via `pip install .`; `import celeris` + `examples/python/01`
  reproduce the diffraction-limited focus (Strehl 0.64, FWHM at λf/D).

## Prior work (pre-changelog, see git history)

The engine, design paths, analysis battery, GUI, Python bindings, and published-device
reproductions (Khorasaninejad 2016, Chen 2018) were developed before this changelog was
started. See `git log` and `ROADMAP.md` (§P0, "What's built") for the full record.
