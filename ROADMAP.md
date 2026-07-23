# CELERIS — Roadmap

**The goal is a paper *portfolio* out of one codebase — done the way that adds up
instead of backfiring.** The rule that makes it work: **one paper per distinct piece
of content.** The same "here is CELERIS" paper cannot go to four journals (that's
redundant publication — a red flag, not a flex). But the *code paper*, a *method
paper*, and several *results papers* are genuinely different content, so they are
genuinely different papers. That's how the codebase becomes a factory.

`[x]` done · `[ ]` todo · `[~]` partial.

> **Licensing (decided):** permissive open-source — recommend **Apache-2.0**
> (permissive + explicit patent grant; MIT/BSD-3 are simpler and equally acceptable).
> Supersedes the old "closed IP for a sale" stance; a permissive license still allows a
> later sale/relicense. A `LICENSE` file gates every venue below — do it first (§0).

---

## The paper portfolio (this is the whole strategy)

| # | Paper | Venue | Content (why it's distinct) | Author | Status |
|---|-------|-------|------------------------------|--------|--------|
| **P1** | **Code paper** | **SciPost Physics Codebases** | "Here is CELERIS, the end-to-end pipeline." | **Solo** | Primary near-term target |
| **P2** | **Method paper** | **CPC** | "A from-scratch **differentiable/adjoint** RCWA." | Solo / +advisor | **Contingent on the adjoint (§3) landing** |
| **R1** | Results | Optics Express / Nanophotonics | A single-etch achromat pushing the Chen-2018 group-delay limit | +advisor | Needs a new result |
| **R2** | Results | Optics Express / Optica | Wide-FOV doublet / aplanatic design | +advisor | Needs a new result |
| **R3** | Results | Optics Express / ACS Photonics | An inverse-designed device (rides on the adjoint) | +advisor | Needs adjoint + result |

**Read this table honestly:**
- **P1 (SciPost code paper) is the guaranteed, solo, near-term win.** No calendar gate,
  no fee, physicist-respected. **Submission-ready in ~2–4 weeks** of real work; arXiv
  preprint the same day; SciPost acceptance is a normal multi-month review after that.
- **P2 (CPC) only exists as a *separate* paper if the adjoint (§3) lands.** The adjoint
  is what makes CPC *distinct content* from P1 ("differentiable RCWA method" ≠ "here is
  the tool"). **If the adjoint doesn't land, do NOT also submit CPC** — it would overlap
  P1 and read as double-publishing. In that case the portfolio is P1 + results papers.
- **R1–R3 are where the volume legitimately lives**, but each needs an *actual new
  result*, not a reprint — that's real work, usually **co-authored with an advisor**
  (Eikenberry?), which is a feature: results papers are how the champion + comparative
  letters get built. The **solo flex stays P1.**

> **Reality check.** These are solid CV bricks that get cold emails **read and answered**
> — a peer-reviewed solo code paper + preprints in a subfield is most of what most
> undergrads never manage. They are not, by themselves, what makes a group fight over
> you. The paper earns the read; the summer, the vouching, and GRFP earn the fight.

---

## ⭐ Milestones

**P1-SUBMITTABLE (the near-term target)** = §1 (shared infra) + §2 (SciPost code paper:
docs, community files, benchmarks, `paper.md`). No calendar gate. Ship the **arXiv
preprint the same day** — that's the instant, citable win.

**P2-SUBMITTABLE** = §1 + the adjoint landing (§3d) + §3a–c (novelty framing, S4
validation, manuscript). **Contingent:** only pursue if the adjoint makes it genuinely
distinct from P1.

---

## 0. Decisions (make these FIRST)

- [ ] **License** → permissive; recommend **Apache-2.0** (add `LICENSE` + SPDX headers).
      Repo currently says "TBD." **Every venue is blocked until this file exists** —
      SciPost, CPC, arXiv-with-code, all of it.
- [ ] **Get an ORCID** (free, 5 min) — required/expected by SciPost and CPC.
- [ ] **Commit steadily from here on.** The whole history is one June-2026 burst; a
      clean, distributed commit history over the coming months reads as a maintained
      project (matters for reviewers judging "is this real software"). Every §1/§3 task
      = its own dated commit, not another sprint.

---

## 1. Shared infrastructure — EVERY paper needs this (HIGHEST LEVERAGE)

Today the repo is MSVC / AVX2 / Windows-centric with **no LICENSE, no CI, no test
suite**. This is the single most likely point of failure — reviewers are on Linux and
they *will* clone, build, and run it.

- [ ] **`LICENSE` file** (§0) + a short license/copyright header convention.
- [ ] **Linux build, clean, one command.** `cmake` + build on a fresh Ubuntu box with
      documented apt deps. Make portability *true and tested*, not "should work."
- [x] **CPU-only build is first-class.** CUDA must **not** be a build requirement — most
      reviewers have no GPU. VERIFIED on Linux (GCC 16 / CMake 4.4 / Ninja): CPU-only
      configure+build+`selftest` all green; the CLI and Python wheel build with no CUDA/MKL.
- [x] **CUDA builds on Linux too** (was Windows-only CMake). The GPU far-field kernel now
      compiles via `nvcc` and links the Linux static runtime; VERIFIED running on an
      RTX 4070 (`cuda::available()`=true). This keeps the GPU perf claim **reproducible on
      Linux** for the paper's benchmark — but GPU stays an *optional accelerator with CPU
      auto-fallback*, never a build requirement. NOTE for CI: hosted runners have **no
      GPU**, so CI stays CPU-only; GPU timing numbers are produced **locally** on the 4070
      and reported in the §3b benchmark (a compile-only CUDA CI job is optional, later).
- [ ] **CI (GitHub Actions)** — `.github/workflows/` does not exist yet; create it. Every
      push: Linux CPU-only build + `celeris selftest` + run the Python examples. A green
      CI badge is the cheapest credibility signal there is.
- [ ] **Regression/unit test suite.** Wrap the existing `selftest` locked numbers as
      CI-enforced tests + a `pytest` layer over the Python bindings asserting
      diffraction-limited FWHM / energy = 1 / grcwa (and later S4) agreement.
- [ ] **One-command reproduction.** A `paper/reproduce_all.(sh|py)` that regenerates every
      number/table/figure the papers cite. SciPost wants benchmarks reproducible; CPC
      effectively requires it; it makes every rebuttal trivial.
- [ ] **Tagged releases + `CHANGELOG.md`** — cut `v0.x` tags as milestones land (a formal
      release process; reviewers look for it).
- [ ] **Zenodo DOI** for the exact commit each paper describes; pin dep versions (Eigen,
      MKL optional, CUDA optional, pybind11, grcwa/S4 for validation).

## 2. P1 — SciPost Physics Codebases (the code paper — SOLO, near-term)

No calendar gate, no APC, physicist-respected. Once §1 is green this is mostly writing.
SciPost's acceptance checklist is the spec: **benchmarking tests, ≥1 detailed example
application, high-level code standards, a user guide, download/install/run docs, and a
demonstrable need for the community.**

### 2a. Software & repo (SciPost referee checks each)
- [ ] **Benchmarking tests** — your `selftest` + `validate` batteries, exposed and
      documented as the benchmark suite (energy conservation, 2D=1D, RCWA=TMM, grcwa/S4
      agreement, convergence vs M).
- [ ] **Detailed example application** — the **Khorasaninejad 2016 / Chen 2018
      reproductions** are exactly this: a worked, published-device design end-to-end.
      Foreground them as the flagship example.
- [ ] **Install/run docs** a stranger follows on Linux — clone → `pip install .`
      (scikit-build-core / `pyproject.toml` exist) → `import celeris` → first design;
      plus the native CLI/GUI build.
- [ ] **User guide** that contextualizes the tool and its added value (the pipeline gap).
- [ ] **API documentation** — Doxygen for the C++ core + docstrings on the pybind11
      bindings + a short "API overview." Proportionate, not exhaustive.
- [ ] **Community guidelines** — `CONTRIBUTING.md`, `CODE_OF_CONDUCT.md`, issue/PR
      templates, a stated support channel.
- [ ] **Evidence it's used for research** (specific, not aspirational): the reproductions
      themselves + areal-world example; bonus if ≥1 other person runs a design before
      submit (the old "lab-ready" goal — a labmate / NSBP contact, not a metalens expert).

### 2b. `paper.md` / manuscript sections (all expected)
- [ ] **Summary** — high-level, non-specialist; author + ORCID + UCF affiliation.
- [ ] **Statement of Need** — who it's for + the *pipeline gap* solvers don't fill.
- [ ] **State of the Field** — comparison vs S4/grcwa/TORCWA/fmmax/meent/RETICOLO. Build
      the feature table (rows = those + CELERIS; columns = 1D / 2D-vectorial / GPU /
      autodiff / arbitrary-shape / integrated design / analysis battery / achromatic+PB /
      GDSII export / published-device reproduction / GUI). You *win decisively* on the
      design+analysis+fab pipeline. **That asymmetry is the paper.** (Reused by P2.)
- [ ] **Software Design** — architecture + trade-offs (C++23 core, RCWA/FMM + Li
      factorization + Redheffer S-matrix, GPU-far-field-only, pybind11 layer) and why they
      matter for the application.
- [ ] **Research Impact / demonstrable need** — reproductions + example + any external use.
- [ ] **AI-usage disclosure** *(required by JOSS/SciPost norms — non-disclosure is an
      ethical breach).* State plainly CELERIS was developed with AI assistance (Claude
      Code): tools + versions, where used (code gen, tests, docs), and confirm **you
      reviewed, validated, and made all core design decisions.** Honest and fine — just
      mandatory.
- [ ] **References** — Moharam–Gaylord, Li, Liu–Fan/S4, Redheffer, Rumpf,
      Khorasaninejad 2016, Chen 2018 + the software archive DOI.

### 2c. Submit
- [ ] **Post the arXiv preprint** (`physics.optics` + `physics.comp-ph`) the day the draft
      is ready — the instant, citable, cold-email-able win.
- [ ] **Submit to SciPost Physics Codebases** the same week: repo public, `LICENSE` in,
      CI green, tagged release, benchmarks + example + docs in place. Then a normal
      multi-month refereeing cycle (respond to referees within ~2 weeks each round).

## 3. P2 — CPC (the METHOD paper: differentiable/adjoint RCWA) — CONTINGENT

CPC is a *separate* paper from P1 **only** if the adjoint lands — that's what makes it
"a new differentiable-RCWA method" instead of "here's the tool again." Build the adjoint
first (§3d); if it fights, drop P2 and keep P1 + results papers. Never let CPC become a
reprint of the SciPost paper.

### 3a. Novelty & positioning
- [ ] **Frame around the adjoint**: from-scratch RCWA with *analytic adjoint gradients* —
      directly answers "why not TORCWA?" (autodiff) with "I am too, from scratch, and I'm
      the whole pipeline." Reuse the P1 comparison table but now you also **win the
      autodiff column.**
- [ ] **Statement of need** + **kill "why not just use X"** for each competitor.

### 3b. External validation (reviewers WILL demand it)
- [ ] **Cross-check against S4 (Stanford `S⁴`)** — the reference everyone trusts. Match to
      ~1e-4 on a dielectric lamellar grating (published efficiencies), a metal grating, and
      a metalens meta-atom. ~1–2 days, high value. (Also strengthens P1's benchmarks.)
- [ ] **Reproduce a canonical published RCWA benchmark to the digit** (Moharam–Gaylord /
      Li gratings, already cited in `rcwa2d.cpp`). Table: CELERIS vs published vs S4.
- [ ] **Convergence figure** (promote `validate §1`): efficiency & energy-defect vs M.
- [x] **Honest timing benchmark — DONE, and it corrected a bad claim.** The old
      "600–700×" was NOT reproducible: `psfbench` (GPU vs CELERIS's own optimized 16-core
      CPU `compute_psf`, RTX 4070) measures **~1.1–5.8×** for far-field propagation, and
      the speedup *shrinks* with grid size (kernel is memory-bound). README/paper corrected
      to the measured numbers. The 600–700× was almost certainly vs a single-core or an
      early un-parallelized CPU path. REMAINING: benchmark the RCWA *solve* vs S4/grcwa on
      identical hardware (the solve, not propagation, is the real cost — see MKL).

### 3c. Manuscript
- [ ] **CPC "Program Summary" block** (title, licensing, languages, nature of problem,
      solution method, restrictions, running time, dependencies).
- [ ] **Method section** (RCWA + Li + Liu–Fan + Redheffer + Rumpf + **the adjoint
      derivation**) → **Results** (adjoint-driven inverse design + reproductions +
      S4/grcwa agreement + timing) → **Limitations** (§3e).
- [ ] **arXiv update** (new version / new preprint for the method paper) on submit.

### 3d. The adjoint (the gate for P2 existing at all — derisk EARLY)
- [ ] **Prototype on the 1D solver first** — check every analytic gradient against the
      existing central-difference values you already trust.
- [ ] **Then 2D** — hard parts: eigendecomposition backward pass + degenerate/near-
      degenerate eigenvalues. Budget accordingly.
- [ ] **Validate adjoint gradients vs finite differences** as a locked CI test.
- [ ] **Decision point:** 1D+2D land cleanly → write P2 (CPC). Otherwise → **no P2**;
      the adjoint becomes "future work" in P1 and you lean on results papers instead.

### 3e. Honest limitations (turn weakness into credibility — carry into every paper)
- [ ] CPU eigensolve; GPU accelerates far-field propagation only (exact numbers).
- [ ] Curved shapes staircase-converge slowly — confirmed in *both* CELERIS and grcwa.
- [ ] Achromatic multi-height designs imply a grayscale etch; single-etch PB is fabricable.
- [ ] Scope is metalenses/metasurfaces via RCWA — NOT FDTD or ray-tracing (a boundary).

## 4. R1–R3 — results-paper pipeline (optics journals, +advisor)

Each is a *new optical result produced with CELERIS*, distinct content, individually
more impressive to a physics committee than another software paper. **Each needs a
genuinely new result — that's real work, and it needs an advisor's steer on what's worth
publishing.** Placeholders until scoped with a PI:

- [ ] **R1 — Achromatic limit-pusher.** A single-etch achromat that extends the group-
      delay/aperture trade-off past Chen 2018, using the existing achromatic + PB
      machinery. Target: Optics Express / Nanophotonics.
- [ ] **R2 — Wide-FOV doublet / aplanatic.** The parked §P2 work: 2-surface corrector +
      lens, off-axis aberration correction. Classic design-paper format. Target: Optics
      Express / Optica.
- [ ] **R3 — Inverse-designed device.** A device showcase driven by the adjoint (§3d) —
      gated on the adjoint landing. Target: Optics Express / ACS Photonics.
- [ ] **Scope with an advisor first** — which of these is (a) novel enough to publish and
      (b) something a PI wants to co-author. This is the conversation that also starts the
      champion/letter relationship.

---

## Sequencing

1. **Week 1 — §0 decisions** (`LICENSE`, ORCID) + the §2b comparison table & statement of
   need. Cheap, pure thinking; reused by P1 *and* P2.
2. **Weeks 1–2 — §1 shared infra**: `LICENSE`, Linux build, CPU-only path, CI. Unblocks
   every reviewer and every paper at once.
3. **Weeks 2–3 — §1 test suite + one-command repro + Zenodo**; start §2 docs in parallel.
4. **Weeks 3–4 — §2 P1 code paper**: benchmarks + example + docs + `paper.md` →
   **arXiv preprint live → SUBMIT SciPost.** First win lands here (submit fast; acceptance
   is a normal multi-month cycle after).
5. **Weeks 2–8 (parallel) — §3d adjoint**, 1D then 2D. This is the fork that decides
   whether P2 (CPC) exists. Derisk early.
6. **Weeks 4–6 — §3b S4 cross-check + benchmark tables** (strengthens P1's benchmarks
   *and* seeds P2).
7. **If adjoint landed → Weeks 6–10 — §3a–c CPC method paper → SUBMIT CPC + arXiv.**
   If not → skip P2, go to R-pipeline.
8. **Ongoing — R1–R3** as new results materialize, scoped with an advisor.

---
---

# PRODUCT BACKLOG (PARKED)

> Everything below is the original product roadmap toward a Zemax/Lumerical competitor.
> It is **no longer the driver** — it is a record of what's built (paper-results material)
> and a quarry for future work. **Do not let it eat the semester.** The papers (§0–§4)
> come first. Items here that directly feed a paper are cross-referenced above.

The product goal posts were: **lab-ready** = usable by your lab for real designs;
**product-complete** = the Zemax bar (a buyer adopts it with no hand-holding). Scope
stays in the lane: **metalenses / metasurfaces via RCWA** — NOT general ray tracing,
NOT FDTD.

## P0. What's built (done) — the paper's results live here
- [x] **RCWA engine**: 1D (TE/TM, multilayer, S-matrix), 2D vectorial (P·Q), validated
      (energy=1, 2D=1D to 5e-12, RCWA=TMM, RCWA=EMT). Li/Liu–Fan inverse-rule 2D
      factorization (fixed the old Laurent solver that gave T≈0.07 where truth ≈0.96;
      transmittance spread over M 0.6→0.03, energy=1.000000 at every M).
- [x] **1D TE convergence fix** — half-space companion admittance sign (Y=−j·kz); the old
      +j·kz was the wrong branch for evanescent orders, so multi-order R/T split converged
      wrong (energy stayed conserved, masking it). Now subwavelength TE → 0.93333
      (|Δgrcwa|~4e-6); selftest **[5]** locks it. (CPC §3b material.)
- [x] **External cross-check vs grcwa** — square pillar, asymmetric rect (both pols),
      subwavelength TE/TM all match; baked into `selftest [8]`. (CPC §3b — *extend to S4*.)
- [x] **Meta-atom library**: rectangular pillars, fill sweep, 2D (fill_x, fill_y).
- [x] **Arbitrary meta-atom shapes** — 2D solver handles two-material shapes via sampled-
      grid Laurent factorization (`MetaShape` + `rasterize_eps` + `shape_operators`); CLI
      `design --shape circle|cross|ring`. Cross matches grcwa ~1e-3, energy conserved.
      HONEST: directional inverse-rule is unstable for non-separable shapes (→ Laurent,
      slower convergence); curved shapes staircase slowly in both CELERIS and grcwa.
      (CPC §3e limitation.) REMAINING: polygon/shape-aware GDS export; Li normal-vector FFF.
- [x] **Design**: hyperbolic focusing (amplitude-aware); polarization-multiplexed (bifocal).
- [x] **Pancharatnam–Berry (geometric-phase) design** — in-plane atom rotation, full 2×2
      Jones (`solve_jones`), HWP-atom finder, `pbdesign`. Geometric phase exact (~1e-15°),
      inherently diffraction-limited, capped by conversion (~0.64 TiO₂ square). Writes
      rotated-polygon GDS (`write_pb_gds`). Vortex/OAM & arbitrary-profile PB done.
- [x] **Arbitrary phase profiles** — focusing / vortex-OAM / deflector / axicon / freeform
      (`PhaseProfile` + `phase_profile_value()`), shared by both PB and propagation paths;
      `pbdesign --profile …` and `design --profile …`; freeform CGH via `--freeform-file`.
      selftest **[14]**.
- [x] **Achromatic / broadband design** — the high-value application, three paths:
      - Propagation-phase (`celeris achromatic`): dispersive fill×height library +
        two-objective (phase mod 2π **and** group delay) selection; drift 8.0→1.9 µm
        (4.2×). selftest **[15]**.
      - Single-etch (`--single-etch`): (phase, group-delay) plane spanned by **shape
        variety at one height** (square/circle/cross/ring × fill) — fabricable in one etch;
        drift 4.0→1.27 µm (3.2×). selftest **[15b]**.
      - Pancharatnam–Berry (`celeris pbachromatic`, modern recipe): PB rotation sets base
        phase *exactly*, a dispersive birefringent atom supplies group delay → decoupled
        objectives, base-phase RMS ~0. drift 4.64→1.90 µm (2.4×). selftest **[15c]**;
        rotated-rect GDS (`write_pb_rect_gds`). All three exposed in Python + a GUI
        Achromatic panel with a Library selector. (CPC §3e limitation: multi-height = grayscale etch.)
- [x] **Optimizers**: per-pillar Adam (inverse design), period×height (system). (FD-based —
      CPC §3d would replace with adjoint.)
- [x] **Analysis battery**: Strehl/FWHM/encircled, PSF, wavefront (OPD/Zernike/Maréchal),
      MTF, through-focus + caustic, chromatic, spot-vs-field, FOV, tolerance, birefringence,
      focal isolation.
- [x] **Field-resolved grids** (`analyze_field_grid` / `celeris fieldmap`) — Strehl + tan/sag
      FWHM across the full (θx,θy) field; selftest **[21]**.
- [x] **Efficiency breakdown** (`analyze_efficiency` / `celeris efficiency`) — per-order
      de_t/de_r, reflection, absorption = 1−R−T, useful-0th vs stray split; selftest **[20]**.
- [x] **Wide-FOV / off-axis (first cut)** — quadratic-phase lens (`design --profile quadratic`,
      `celeris widefov`): with an offset stop, quadratic stays coma-free where hyperbolic
      develops coma (±≥30° vs ±15.4° to Strehl 0.8); selftest **[17]**. (R2 seed.) REMAINING:
      doublet, aplanatic, GUI panel.
- [x] **Validation battery** (`celeris validate`) on real tabulated TiO₂ n,k (Siefke 2016):
      diffraction-limited focusing (phase Strehl 0.97, FWHM at λf/D), convergence tabulated.
- [x] **Published-device reproductions** (`celeris reproduce`) — **the P1/CPC headline**:
      - **Khorasaninejad 2016** (Science 352, 1190): NA=0.80 PB TiO₂ nanofins, transmitted-
        normalized conversion 96–99% (HWP quality), 2θ phase ~4–5° RMS, FWHM ~0.514·λ/NA;
        selftest **[16]** locks the 532 nm device.
      - **Chen 2018** (Nat. Nanotechnol. 13, 220): NA=0.20 D=26.4 µm broadband achromat,
        470–670 nm; reproduces the central physical limit — required group-delay span
        4.45 fs sits at the ~5 fs 600-nm-nanofin ceiling (diameter is GD-limited), GD
        objective flattens drift 5.0× in one etch; selftest **[18]**.
- [x] **Materials**: named registry (`by_name`/`catalog`) — analytic Sellmeier (SiO₂, BK7,
      Si₃N₄, Al₂O₃, GaN) + real tabulated n,k (TiO₂, c-Si, a-Si, Au, Ag, Al; CC0 from
      refractiveindex.info) in `data/`; `celeris materials`; selftest **[19]**.
- [x] **GPU**: far-field propagation kernel (~1–6× over 16-core CPU, measured — *propagation
      only*, §3b; the old 600–700× claim was unreproducible and has been corrected),
      CPU multicore library build.
- [x] **AVX2/FMA build** (default, ~2.7×) + **Intel MKL turbo** (opt-in, up to ~11× at
      high M; eigensolve is only ~15% — the cost is dense matmuls/inverses/S-matrix).
- [x] **I/O**: GDSII write (square + rectangular + rotated polygons), GDS read + in-app viewer.
- [x] **Python bindings (pybind11)** over celeris_core — full workflow, numpy maps;
      `pip install .` via scikit-build-core; 7 runnable examples. (P1 backbone.)
- [x] **GUI** (Dear ImGui, modularized celeris::gui): dockable workspace, all analyses,
      layout/library/polarization/achromatic panels, report export, project save/load, GPU
      indicator, dark mode, perf monitor.
- [x] **CLI**: design / polardesign / pbdesign / achromatic / pbachromatic / widefov /
      fieldmap / efficiency / materials / validate / reproduce / birefringence / selftest.
- [x] Honest negative findings recorded (GPU eigensolve, reference-LAPACK eigensolve).

## P1. Credibility & correctness (partly promoted to §1/§3b above)
- [x] Convergence studies documented (`validate §1`: metric vs M, vs sampling).
- [ ] Auto-convergence helper (suggest M for a target accuracy).
- [ ] Numerical edge cases hardened (Rayleigh anomalies, grazing orders, resonances).
- [ ] Regression/unit test suite + CI → **promoted to §1 (every paper needs it).**

## P2. Design quality & breadth
- [~] Full-2π libraries via height sweep (`design --auto-height`, `coverage()` metric) —
      TiO₂ square caps ~325–345°; closing the gap needs shape variety or higher index.
- [ ] Full Jones-matrix metasurfaces (arbitrary polarization transforms).
- [~] Wide-FOV → doublet (Arbabi 2016) + aplanatic; GUI panel. **(→ R2 results paper.)**
- [ ] Topology / freeform meta-atom inverse design (overlaps §3d adjoint / **R3**).
- [ ] Fab constraints baked into design (min feature, min gap, etch bias).

## P3. Productization
- [ ] Installer (WiX/NSIS) + code signing.
- [ ] One build, auto CPU/GPU at runtime → CPU-only path is **promoted to §1**; bundling
      is product-only.
- [ ] Example projects / templates gallery.
- [ ] User manual + tutorials + theory notes; **API docs (Doxygen) → promoted to §2 (P1).**
- [ ] Report export to PDF/HTML (today: txt + images).

## P4. Analysis & materials depth
- [ ] Polarization analysis maps (Stokes/Mueller, extinction-ratio map).
- [ ] Stray light / higher-order ghosts.
- [ ] Material editor in GUI (fit Sellmeier, anisotropic, thermo-optic); metals via
      analytic Lorentz–Drude (today tabulated).
- [ ] Fabrication-aware analysis (sidewall angle, corner rounding → performance).

## P5. GUI / UX polish
- [ ] ImPlot interactive charts (zoom/pan/hover).
- [ ] Parametric sweep / optimization UI, merit editor.
- [ ] Undo/redo, units management, validation, tooltips.
- [ ] 3D view (pillar array, field volume).
- [ ] Compare-designs side by side; recent files; drag-drop; DPI scaling.

## P6. I/O & formats
- [ ] GDS layers/datatypes + cell hierarchy (multi-material / multi-height masks).
- [ ] More export formats: DXF, OASIS, CIF; optional STL/3D.
- [ ] Import-and-analyze (GDS/phase-map → performance); foundry round-trip.
- [ ] Process/stack definition files.

## P7. Performance & scale
- [ ] GPU-batch the whole library/wavelength sweep in one launch.
- [ ] Out-of-core / mm-scale apertures (millions of pillars), memory streaming.
- [ ] Library/result caching; incremental recompute.
- [ ] Cancellable long runs; multi-GPU (later).

## P8. Ops / platform
- [x] **Linux build (CPU + CUDA) → DONE & verified (§1).** CPU-only and CUDA-on-Linux both
      build and run (RTX 4070); `selftest` green on GCC 16. CI workflow added. macOS optional.
- [ ] Versioning, release notes, changelog; crash reporting; auto-update.
- [ ] Licensing/activation if commercial → **superseded by the permissive-license decision.**
