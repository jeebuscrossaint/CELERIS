#!/usr/bin/env bash
# reproduce_all.sh — regenerate every number/table in paper.md from scratch on a
# clean Linux CPU-only toolchain (the environment a reviewer uses). Outputs land in
# paper/artifacts/. The optional GPU propagation benchmark runs only if a CUDA build
# is available (reviewers without a GPU can skip it; nothing else needs one).
#
# Usage:  bash paper/reproduce_all.sh            # CPU-only: physics + reproductions
#         CELERIS_REPRO_CUDA=1 bash paper/reproduce_all.sh   # also GPU psfbench
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"                       # data/ must be reachable from the cwd
BUILD="${BUILD:-build-repro}"
ART="paper/artifacts"
mkdir -p "$ART"

echo "==> [1/4] configure + build (CPU-only: no CUDA/MKL/GUI)"
cmake -B "$BUILD" -G Ninja -DCMAKE_BUILD_TYPE=Release \
      -DCELERIS_BUILD_GUI=OFF -DCELERIS_USE_CUDA=OFF -DCELERIS_USE_MKL=OFF
cmake --build "$BUILD" --target celeris -j
CEL="$ROOT/$BUILD/celeris"

echo "==> [2/4] physics self-test (locked regression suite)"
"$CEL" selftest 2>&1 | tee "$ART/selftest.txt"

echo "==> [3/4] credibility battery on real tabulated TiO2 n,k"
"$CEL" validate 2>&1 | tee "$ART/validate.txt"

echo "==> [4/4] published-device reproductions (Khorasaninejad 2016, Chen 2018)"
"$CEL" reproduce --device all      2>&1 | tee "$ART/reproduce_khorasaninejad2016.txt"
"$CEL" reproduce --device chen2018 2>&1 | tee "$ART/reproduce_chen2018.txt"

# Optional: regenerate the paper figures and run the LIVE external-solver
# cross-check. These need the Python module (`pip install .`) + matplotlib/numpy,
# and the S4 panel/script needs a built S4; each is skipped gracefully when its
# dependency is absent, so a reviewer without them still reproduces all the physics
# above (the selftest already locks the grcwa/S4 reference values as literals).
echo "==> [opt] figures + external-solver (S4) cross-check (Python; auto-skips if deps absent)"
if command -v python3 >/dev/null 2>&1 && python3 -c "import celeris" >/dev/null 2>&1; then
  if python3 -c "import matplotlib, numpy" >/dev/null 2>&1; then
    python3 paper/figures/make_figures.py && \
      echo "    figures regenerated in paper/figures/ (S4 panel auto-falls back to CELERIS-only if S4 absent)"
  else
    echo "    SKIP figures: matplotlib/numpy not installed"
  fi
  if python3 -c "import S4" >/dev/null 2>&1; then
    python3 paper/validation/s4_crosscheck.py 2>&1 | tee "$ART/s4_crosscheck.txt"
  else
    echo "    SKIP live S4 cross-check: S4 not importable (selftest [5] locks the same reference values)"
  fi
else
  echo "    SKIP figures + S4 cross-check: the 'celeris' Python module is not installed (run: pip install .)"
fi

# Optional: GPU far-field propagation benchmark (needs a CUDA build + device).
if [[ "${CELERIS_REPRO_CUDA:-0}" == "1" ]]; then
  echo "==> [opt] GPU propagation benchmark (psfbench)"
  CBUILD="${CBUILD:-build-repro-cuda}"
  cmake -B "$CBUILD" -G Ninja -DCMAKE_BUILD_TYPE=Release \
        -DCELERIS_BUILD_GUI=OFF -DCELERIS_USE_CUDA=ON -DCELERIS_USE_MKL=OFF
  cmake --build "$CBUILD" --target celeris -j
  {
    for g in 161 321 641 1024; do "$ROOT/$CBUILD/celeris" psfbench --grid "$g"; done
    for d in 200 300;          do "$ROOT/$CBUILD/celeris" psfbench --diameter "$d" --grid 161; done
  } 2>&1 | tee "$ART/psfbench.txt"
fi

echo "==> done. Artifacts written to $ART/"
