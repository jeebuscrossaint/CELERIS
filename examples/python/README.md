# CELERIS Python examples

These scripts drive the validated CELERIS C++ engine from Python (numpy +
matplotlib). They cover the core metalens workflow a new user wants first.

| Script | What it shows |
|--------|---------------|
| `01_focusing_lens.py`     | Design a TiO₂-on-glass focusing metalens; print Strehl/FWHM/efficiency; plot the PSF + pillar fill map. |
| `02_rcwa_convergence.py`  | 2D RCWA zeroth-order transmittance + energy vs harmonic count `M` (the solver-credibility check). |
| `03_material_dispersion.py` | Load real tabulated TiO₂ n,k (Siefke 2016) and see how phase coverage shifts across the visible band. |
| `04_height_sweep.py`      | Full-2π search: sweep pillar height to maximize phase coverage at high transmittance. |

## Setup

Build the extension once (from the repo root):

```bash
cmake -B build -DCELERIS_BUILD_PYTHON=ON
cmake --build build --config Release --target _celeris
```

This stages an importable package at `build/python/celeris` (the `.pyd` plus any
MKL/CUDA runtime DLLs). Then either add it to `PYTHONPATH`:

```bash
# Windows PowerShell
$env:PYTHONPATH = "$PWD\build\python"
python examples\python\01_focusing_lens.py

# bash
PYTHONPATH=build/python python examples/python/01_focusing_lens.py
```

…or install the package into your environment with pip (uses scikit-build-core):

```bash
pip install .
python examples/python/01_focusing_lens.py
```

Each plotting script saves a PNG next to itself; set `MPLBACKEND=Agg` for a
headless run.
