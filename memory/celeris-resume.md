# CELERIS — session resume (handoff)

**Last updated:** 2026-08-07. **HEAD at handoff:** `34b130c` (run `git log --oneline -25`).
**Resuming in ~2 weeks (~late Aug 2026):** everything below is current; the next action is
the v1.0 release (see "NEXT STEP #4" — starts with YOU enabling the repo on Zenodo).

## Update 2026-08-07
- **CPC eligibility audited — no JOSS/SciPost-style trap.** Confirmed against Elsevier's
  guidelines: (1) NO author-level gate — a solo undergrad can be corresponding author;
  (2) NO timing clock (unlike JOSS's 6-mo repo rule); (3) **Apache-2.0 is on CPC's approved
  open-source license list** — the old "CPC Non-Profit Use License Agreement" is the LEGACY
  (pre-2016) program-library model, NOT layered on new Mendeley-Data deposits, so commercial
  optionality is preserved. We are a **CPiP** (Computer Programs in Physics) submission.
  Two real (non-blocking) caveats: review is SLOW (~12 mo to first decision → the arXiv
  preprint does the near-term recruitment work); the real desk-reject risk is novelty/scope
  ("another RCWA solver"), mitigated by the pipeline framing. TODO before submit: a 5-min
  primary read of the Guide for Authors (ScienceDirect 403s bots) to confirm the Program
  Summary field list verbatim + the Mendeley deposit step.
- **CLI refactor (commit `06f547e`).** `main.cpp` was a 3325-LOC monolith (25% of the
  codebase). Split into a 33-line dispatcher + one TU per command family under `cli/`
  (`cli_common`, `selftest`, `cmd_design`, `cmd_polar`, `cmd_achromatic`, `cmd_reproduce`,
  `cmd_validate`, `cmd_materials`, `help`, `bench`); shared decls/includes in `cli/cli.hpp`.
  Verbatim move — behavior identical; builds clean (Ninja/GCC, GUI off), selftest passes all
  22 cases. Build dir used on this box: `build-linux` (CPU-only, no MKL → selftest ~3 min).
- **Paper `[TODO]`s RESOLVED + LaTeX conversion DONE** (commit `03297d8`). `paper/paper.md`
  is now TODO-free (stale SciPost research-impact bullet dropped; AI-disclosure slimmed to a
  one-sentence Elsevier-style "Declaration of generative AI"; optional-figure markers removed).
  **`paper/paper.tex` written in Elsevier `elsarticle`** — full conversion (frontmatter, CPC
  Program Summary block, comparison table, 4 figures, all `\cite` wired) and it **compiles
  clean**: `pdflatex paper && bibtex paper && pdflatex×2` → 10-page PDF, 0 undefined refs.
  Build needs `elsarticle.cls` (NOT in this minimal TeX Live; `texlive-publishers` provides it,
  or regenerate: `latex elsarticle.ins`). I fetched+gitignored the class + `elsarticle-num.bst`
  in `paper/`. **Bib fixes:** TORCWA title was WRONG (fixed); fmmax filled (Schubert & Hammond,
  arXiv:2308.08573); **stripped inline `%` comments after entry keys — BibTeX has no `%`
  comments so they were silently skipping 5 cited entries.** Remaining bib nit: empty `year`
  on `grcwa`/`meent` @misc (cosmetic warnings; verify in the final citation pass).
- **CI fixed + selftest parallelized** (commits `12d058e`, `1dedc9b`). CI had been RED on
  every push (pre-existing, not the refactor): `ubuntu-latest` default g++-13 lacks the C++23
  `<print>` header → now installs/uses **g++-14**. Also parallelized the self-test: it was
  ~serial (16 cores gave only 1.7x). Case [8] (2D RCWA convergence) alone was 85% of runtime;
  its independent solves now run via `std::async`, and the 22 cases run on a thread pool with
  per-case output buffers. **276s → 88s (3.1x) on 16 cores, stdout byte-identical** to the
  serial suite (diff-verified vs a golden run — physics gate intact). Added
  **`selftest --quick`** (fast subset [1]-[7]+[19], <1s); CI runs `--quick` on push/PR
  (build+test now **1m20s**, was failing ~9m43s) and the **full 22-case suite nightly**
  (schedule trigger, 07:00 UTC). NOTE the golden reference lives at scratchpad `st2.log`
  (session-scoped — regenerate via `./build-linux/celeris selftest` if re-verifying).
  Real remaining speed lever for the RCWA solve is MKL (build-config), not the selftest.
- **Release prep DONE** (commits through `15c68d0`): `paper.bib` warning-free (grcwa 2020,
  meent arXiv:2406.12904 2024 filled); `CHANGELOG.md` `[Unreleased]` updated for v1.0
  (retarget→CPC, CLI split, paper.tex, --quick, CI g++-14 fix, selftest parallelization).
  Working tree clean except the pre-existing loose `README.md`/`ROADMAP.md`/`imgui.ini`.
- **v1.0.0 RELEASED + archived (2026-08-11, commit `c38503c`).** GitHub release `v1.0.0`
  cut via `gh release create` (triggers Zenodo). Zenodo minted: **version DOI
  `10.5281/zenodo.21891346`** (this exact code), **concept DOI `10.5281/zenodo.21891345`**
  (always latest). Wired in: `paper.tex` Program Summary (version+concept), README badge +
  `CITATION.cff` (concept DOI, version 1.0.0, date-released). CHANGELOG stamped `[1.0.0]`.
- **NEXT STEP = submission (HUMAN — needs Amarnath's accounts):**
  1. Submit to **CPC** via Elsevier Editorial Manager (upload `paper.tex` + `paper.bib` +
     `paper/figures/*.pdf`; solo corresponding author; `@ucf.edu` → APC waiver). At submission
     CPC creates the Mendeley Data / CPC Program Library deposit → fills the Program Summary
     "CPC Library link" field.
  2. Post the **arXiv** preprint (`physics.comp-ph` + cross-list `physics.optics`) same day —
     NEEDS A FACULTY ENDORSER for physics.comp-ph; line this up in advance.
  3. Fallback if CPC rejects → SoftwareX (minimal changes).
- Loose in working tree (pre-existing, untouched): `ROADMAP.md`, `imgui.ini`.

## Where we are in one line
CELERIS has a **content-complete first-paper draft** (`paper/paper.md`, 4 figures, verified
against grcwa + Stanford S4). The venue is now **Computer Physics Communications (CPC)**.
The remaining work is: convert the draft to Elsevier `elsarticle` LaTeX, tag a release +
Zenodo DOI, and submit. Software + validation are done and committed.

## THE VENUE DECISION (this changed twice — do not re-litigate)
- Goal: a first publication for grad-school recruitment. Author is an **undergraduate**
  (sophomore), so venue must allow **solo undergrad submission**.
- **JOSS** → dropped: 6-month public-repo clock (repo public ~2026-06-16 → not submittable
  until ~Dec 2026) + light for physics.
- **SciPost Physics Codebases** → **DEAD for us**: registration to submit is restricted to
  "professional academics (PhD students and above)" — an undergrad cannot be the submitting
  Contributor. (A PhD+ co-author could submit, but we're going solo-capable.)
- **CPC (Computer Physics Communications)** → **THE TARGET.** Elsevier, IF~7, the
  computational-physics software journal; **no academic-level gate** (undergrad can be solo
  corresponding author); fits the author's "computational physicist" identity.
- **SoftwareX** → **pre-planned FALLBACK** if CPC rejects (Elsevier, IF~3, lower novelty
  bar, no gate). Paper barely changes between CPC and SoftwareX.
- **APC:** UCF has an **Elsevier Read & Publish agreement** that waives OA APCs for UCF
  authors (auto via `@ucf.edu`; excludes society/Cell/Lancet; through 2028). So CPC/SoftwareX
  APC is effectively free for this author. CPC is also hybrid → subscription route is free too.
- A faculty **co-author is now OPTIONAL** (not required for CPC) but strategically valuable
  (champion + recommendation letter). Author leans solo for now.

## Author identity (baked into LICENSE/CITATION/CoC + git)
- Paper/citation author: **Amarnath Patel**, University of Central Florida, `am397421@ucf.edu`,
  ORCID **0009-0008-9460-082X**.
- Git commits authored as `jeebuscrossaint <thediamond270@gmail.com>` (matches existing
  history; set in repo-local git config).
- Repo: https://github.com/jeebuscrossaint/CELERIS (public).

## What shipped this session (all COMMITTED — persists)
- **Linux CPU-only build verified** (GCC 16 / CMake 4.4 / Ninja) — the whole plan rested on
  this; it builds + `selftest` passes (all 21 cases) + Python wheel installs + example runs.
- **3 real portability/perf bugs FIXED:** (1) pybind `solve_rcwa_2d` arg-count mismatch
  (broke Linux pip build); (2) `celeris_core` not `-fPIC` (broke the Python `.so` link);
  (3) MKL thread oversubscription — `eig.cpp` now pins `mkl_set_num_threads(1)` by default
  (MKL was ~5x SLOWER before; now **2.3x faster than AVX2** at M=8).
- **CUDA made cross-platform** (was Windows-only CMake): now builds on Linux (`/opt/cuda`,
  `libcudart_static.a` + `libculibos.a`), verified running on an RTX 4070.
- **Corrected a false claim:** README/ROADMAP said GPU far-field was "600-700x"; measured
  reality via `psfbench` is **~1-6x vs 16-core CPU** (kernel is memory-bound). Fixed everywhere.
- **S4 cross-check:** built Stanford S4 (phoebe-p fork) on Linux and matched CELERIS to
  **7e-8 (TE)** on an identical grating; three-way agreement CELERIS = grcwa = S4 = 0.93333.
- **Open-source/paper infra:** Apache-2.0 `LICENSE`, `CITATION.cff`, `CONTRIBUTING.md`,
  `CODE_OF_CONDUCT.md`, `CHANGELOG.md`, GitHub Actions CI (`.github/workflows/ci.yml`).
- **Paper:** `paper/paper.md` (full draft) + `paper/paper.bib`; `paper/reproduce_all.sh`;
  `paper/validation/s4_crosscheck.py`; `paper/figures/make_figures.py` + **4 committed PDFs**
  (fig_validation, fig_library, fig_psf, fig_achromatic). Competitor comparison table
  researched + verified. Honest limitations + AI-usage disclosure written.

## !!! SESSION-ISOLATION WARNING (read before continuing) !!!
The previous session's **scratchpad is GONE** in a new session (it lives under a
session-specific `/tmp/claude-.../<SESSION_ID>/scratchpad/`). That means the **built venv,
the installed `celeris`/`S4`/`matplotlib`/`mkl` wheels, and the S4 source build are NOT
available** to you. What PERSISTS: the whole repo (all commits above), and the
**system packages** installed via pacman: **openblas, suitesparse, boost** (system-wide).
To reproduce the Python/S4/figure work you must recreate a venv:
```
python -m venv .venv && . .venv/bin/activate
pip install .                       # the celeris module (needs network for Eigen/pybind)
pip install matplotlib numpy mkl-devel mkl-include
# S4 (optional, for the cross-check/validation figure): clone github.com/phoebe-p/S4,
#   patch NumPy-2: cast PyArray_ENABLEFLAGS((PyArrayObject*)Earr/Harr,...) in S4/main_python.c,
#   then: make S4_pyext BOOST_PREFIX=/usr BLAS_LIB=-lopenblas LAPACK_LIB=-lopenblas \
#         FFTW3_LIB=-lfftw3 CHOLMOD_LIB=-lcholmod CHOLMOD_INC=-I/usr/include/suitesparse \
#         PTHREAD_LIB=-lpthread   (see paper/validation/s4_crosscheck.py header)
```
The committed figures/scripts already capture the results, so you do NOT need to rerun any
of this unless regenerating figures or extending the validation.

## NEXT STEPS (in order) for the CPC submission
1. **Author reads `paper/paper.md`** and signs off on every claim (was the standing gate).
2. **Fill remaining paper `[TODO]`s:** AI-disclosure tool versions/date; optionally a
   Khorasaninejad-2016 efficiency figure (nice-to-have; pipeline already shown by Figs 2-4).
3. **Convert `paper/paper.md` → Elsevier `elsarticle` LaTeX** (`paper/paper.tex`), NOT SciPost.
   The CPC "Program Summary" block is already drafted at the top of `paper.md` — keep it.
   `paper/paper.bib` is BibTeX (verify the `[check]` competitor citations).
4. **Tag `v1.0` + mint a Zenodo DOI** for the exact commit; add the DOI to the paper + README.
5. **Submit to CPC** via Elsevier Editorial Manager (author is solo corresponding author;
   `@ucf.edu` for the APC waiver). Post the **arXiv preprint** (`physics.comp-ph` +
   `physics.optics`) the same day — needs a faculty endorsement for arXiv, so line that up.
6. If CPC desk-rejects/rejects → **SoftwareX** (fallback; minimal changes).

## Follow-on (future papers, NOT this one — keep scope off the SciPost paper)
- **Paper 2 (method, CPC or Optics Express):** the analytic **adjoint** (differentiable RCWA)
  — closes the one column CELERIS loses to fmmax/TORCWA/meent (autodiff). Derisk 1D first.
- **Arbitrary-polygon FFF** (Li normal-vector) to match S4's shape versatility.
- **GPU-batched library solve** (the real "make it fast" win; the current GPU only does
  far-field propagation, not the solve — MKL is the current solve-side lever).
- **Results papers** (Optics Express/Nanophotonics): achromat limit-pusher, wide-FOV doublet,
  inverse-designed device — co-authored, the source of letters.
- RULE: one code paper per venue; distinct content per paper; never re-publish the same paper.

## Honest framing to preserve (this is the credibility core)
- CELERIS is a **metalens design *pipeline*** (solver → library → design → analysis → GDSII),
  NOT the best raw solver. The GPU/autodiff-native tools (fmmax/TORCWA/meent) beat its
  *kernel*; CELERIS wins on the integrated, validated, fab-ready workflow none of them ship.
- Every number in the paper is measured and reproducible. Do NOT reintroduce hype (no
  "600-700x", no "beats Stanford"). The honest line is "matches S4 to 7e-8, from scratch,
  plus the whole pipeline." Do not fake commit dates, users, or results.
