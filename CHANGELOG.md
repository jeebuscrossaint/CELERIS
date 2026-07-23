# Changelog

All notable changes to CELERIS are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project aims to
follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Apache-2.0 `LICENSE`, `CITATION.cff`, `CONTRIBUTING.md`, `CODE_OF_CONDUCT.md`, and
  this changelog — publication/open-source readiness (SciPost Physics Codebases, arXiv).
- GitHub Actions CI: Linux CPU-only build + `celeris selftest` + Python examples.
- **Cross-platform CUDA build.** The GPU far-field kernel now builds on Linux, not just
  Windows: CMake detects `/opt/cuda` (or `CUDA_PATH`) via `lib64/libcudart_static.a`,
  drives `nvcc` with the correct host-compiler flags per platform (`-Xcompiler=-fPIC,…`
  on Linux vs `/MD,/EHsc,…` on Windows), and links the Linux static runtime
  (`libcudart_static.a` + `libculibos.a` + `dl`/`rt`/`pthread`).

### Fixed
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
- **CUDA path runs on Linux hardware** — `cuda::available()` true,
  `device_name()` = "NVIDIA GeForce RTX 4070 Laptop GPU". (GPU accelerates far-field
  propagation only; the RCWA eigensolve remains CPU — stated precisely for the paper.)
- Python wheel builds via `pip install .`; `import celeris` + `examples/python/01`
  reproduce the diffraction-limited focus (Strehl 0.64, FWHM at λf/D).

## Prior work (pre-changelog, see git history)

The engine, design paths, analysis battery, GUI, Python bindings, and published-device
reproductions (Khorasaninejad 2016, Chen 2018) were developed before this changelog was
started. See `git log` and `ROADMAP.md` (§P0, "What's built") for the full record.
