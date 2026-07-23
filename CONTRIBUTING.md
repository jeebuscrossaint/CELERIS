# Contributing to CELERIS

Thanks for your interest in CELERIS — a from-scratch RCWA metalens/metasurface design
and analysis pipeline. This guide covers how to build, test, and propose changes.

## Ways to contribute

- **Report a bug or ask a question** — open a [GitHub issue](https://github.com/jeebuscrossaint/CELERIS/issues).
  Include your OS, compiler version, how you built (CPU-only vs GPU), and a minimal
  command or script that reproduces the behavior.
- **Propose a change** — open a pull request (see below).
- **Suggest a feature** — open an issue describing the design problem it solves; note
  that CELERIS's scope is deliberately **metalenses/metasurfaces via RCWA** (not FDTD,
  not general ray tracing).

## Building from source (Linux, CPU-only)

CELERIS needs a C++23 compiler, CMake ≥ 3.20, and network access at configure time
(Eigen is fetched automatically). CUDA and Intel MKL are **optional** accelerators —
the default build is CPU-only and fully functional.

```bash
cmake -B build -G Ninja -DCELERIS_BUILD_GUI=OFF -DCELERIS_USE_CUDA=OFF -DCELERIS_USE_MKL=OFF
cmake --build build --target celeris -j
./build/celeris selftest      # runs the locked physics regression suite
```

To build the Python bindings:

```bash
pip install .                 # scikit-build-core; exposes the `celeris` module
python examples/python/01_focusing_lens.py
```

The desktop GUI (`-DCELERIS_BUILD_GUI=ON`) additionally requires GLFW/OpenGL/X11 dev
libraries and is not needed for the engine, the CLI, or the tests.

## Running the tests

The physics is locked by the built-in self-test, which every change must keep green:

```bash
./build/celeris selftest
```

`selftest` verifies the invariants that define correctness — energy conservation,
2D-solver == 1D-solver agreement (~5e-12), RCWA == TMM, cross-checks against the
external `grcwa` solver, and the reproduced published devices. **A PR that changes
engine behavior must not break any locked case**; if a number legitimately changes,
update the locked value in the same PR and explain why.

## Pull request guidelines

1. Keep changes focused; one logical change per PR.
2. Match the surrounding code's style, naming, and comment density.
3. Ensure a clean CPU-only Linux build and a green `selftest` before pushing.
4. Describe *what* changed and *why*, and note any physics numbers that moved.
5. New engine capabilities should come with a locked `selftest` case where practical.

## Use of AI assistance

CELERIS is developed with AI coding assistance. If your contribution was produced with
AI help, please note it in the PR description (tool and where it was used), and confirm
you reviewed and validated the output. See the project's paper for the full disclosure.

## License

By contributing, you agree that your contributions are licensed under the project's
[Apache License 2.0](LICENSE).
