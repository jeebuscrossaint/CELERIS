# CELERIS — Roadmap to "done"

The goal posts: **lab-ready** = good enough to hand to people in your lab and have
them use it for real designs. **product-complete** = the Zemax bar — a buyer can
adopt it with no hand-holding. This file is the full punch list to get from one to
the other. `[x]` = done, `[ ]` = todo, `[~]` = partial.

Scope stays in the lane: **metalenses / metasurfaces via RCWA**. Explicitly NOT
general ray tracing (Zemax's lane) or FDTD (Lumerical's lane).

---

## 0. Where we are today (done)
- [x] RCWA engine: 1D (TE/TM, multilayer, S-matrix), 2D vectorial (P·Q), validated
      (energy=1, 2D=1D to 5e-12, RCWA=TMM, RCWA=EMT)
- [x] Meta-atom library: rectangular pillars, fill sweep, 2D (fill_x, fill_y)
- [x] Design: hyperbolic focusing, amplitude-aware; polarization-multiplexed (bifocal)
- [x] Optimizers: per-pillar Adam (inverse design), period×height (system)
- [x] Analysis battery: Strehl/FWHM/encircled, PSF, wavefront (OPD/Zernike/Maréchal),
      MTF, through-focus + caustic, chromatic, spot-vs-field, FOV, tolerance,
      birefringence, focal isolation
- [x] GPU: far-field propagation kernel (600–700× over CPU), CPU multicore library build
- [x] Materials: Sellmeier / tabulated / constant, CSV import, built-ins
- [x] I/O: GDSII write (square + rectangular), GDS read + in-app viewer
- [x] GUI: dockable workspace, all analyses, layout/library/polarization panels,
      report export, project save/load, GPU indicator, live progress
- [x] CLI: design / polardesign / birefringence / selftest / benchmarks
- [x] Honest negative findings recorded (GPU eigensolve, reference-LAPACK eigensolve)

That's a strong, validated **engine + analysis + GUI**. What's missing is mostly
*credibility*, *design quality/breadth*, and *productization*.

---

## 1. Credibility & correctness  (makes it TRUSTED)
- [x] **Validation vs a published, fabricated metalens** — DONE. Two prongs:
      (a) `celeris validate` runs the credibility battery on REAL tabulated TiO₂ n,k
      (Siefke 2016 ALD, in `data/`) on fused silica — a real-TiO₂ focusing lens is
      **diffraction-limited** (phase Strehl 0.97, FWHM exactly at λf/D, NA 0.20), the
      meta-atom transmittance **converges** (Li factorization) and matches grcwa.
      (b) **`celeris reproduce`** reproduces the canonical visible device,
      **Khorasaninejad et al., Science 352, 1190 (2016)** — NA=0.80 PB-phase TiO₂
      nanofins, H=600 nm, designed at 405/532/660 nm (reported focusing efficiency
      86/73/66 %). Solves the PUBLISHED nanofin geometry (W,L,U from the paper) with
      the exact ALD-amorphous-TiO₂ n,k they deposited, and reports three efficiency
      numbers: transmittance (reflection-limited ~62–66 %), absolute spin-conversion
      (~60–66 %), and **transmitted-normalized conversion = HWP quality (96–99 %)** —
      near-ideal half-wave plates whose transmitted-normalized conversion correctly
      brackets the paper's focusing efficiency from ABOVE in all three; the absolute
      number sits lower only by the Fresnel reflection the paper's transmission-
      referenced efficiency normalizes out. RCWA-verifies the 2θ geometric phase
      (~4–5° RMS) and a diffraction-limited focal spot (FWHM ~ 0.514·λ/NA). HONEST
      FINDING: the 660 nm nanofin (fill_y=0.95, a 20 nm gap) needs M≥10 to converge
      (default bumped to 10). selftest **[16]** locks the 532 nm device. REMAINING
      (lower value): a published ACHROMATIC device (K2016 is monochromatic); shaped /
      higher-index meta-atoms to push absolute (reflection-limited) efficiency higher.
- [x] **Improved (Li/Liu–Fan) 2D Fourier factorization** — DONE. Replaced the basic
      Laurent 2D solver (which gave qualitatively WRONG transmittance, e.g. T≈0.07 where
      truth is ≈0.96, and broke energy beyond M≈8) with the inverse-rule / Liu–Fan
      eigenproblem (analytic directional ε_xx/ε_yy + ⟦1/ε⟧ for axis-aligned rectangles).
      `validate`'s transmittance spread over M now 0.6 → **0.03**, energy = 1.000000 at
      every M. Cross-checked against grcwa (see next): square pillar, asymmetric rect,
      and subwavelength grating all match.
- [x] Cross-check vs an external solver — DONE (grcwa, pip-installed, validated vs analytic
      TMM). CELERIS 2D matches grcwa on the square pillar, asymmetric rect (0.964/0.977 vs
      0.954/0.972, both pols), and subwavelength TE/TM. Baked into `selftest` [8].
- [x] **1D TE solver convergence** — FIXED. The 1D TE transmission for a subwavelength
      grating (Λ=0.3, λ=0.5) used to oscillate/diverge with M instead of converging to
      the grcwa-true 0.93333 (TM was always fine). ROOT CAUSE: a SIGN error in the
      half-space companion admittance for E_y. The layer modes propagate as
      exp(−k0·q·z) (forward wave e^{−i·kz·z}, companion ∂_z = −i·kz), so the matching
      half-space admittance must be **Y = −j·kz**, but the code used +j·kz. The +j·kz
      sign is absorbed harmlessly for the propagating order-0, but it is the wrong
      branch for the EVANESCENT orders, so the R/T split converged to the wrong value as
      more orders were kept (energy stayed conserved throughout — the boundary solve is
      structurally unitary — which masked the bug; the old selftest only checked energy).
      The diagnosis went deeper than the original note: the split was actually wrong for
      MULTI-order gratings too (not just TE oscillation), confirmed vs grcwa. Fix is a
      one-line sign flip in BOTH 1D paths (single-layer `rcwa1d.cpp` + S-matrix stack
      `rcwa_stack1d.cpp`); TM untouched. NOW: subwavelength TE → 0.93333 (|Δ grcwa|
      ~4e-6), multi-order Λ=1.0 matches grcwa exactly at λ=0.6 (0.79332 vs 0.79335) and
      λ=0.45 (0.89926 exact), normal & oblique, energy=1. selftest **[5]** now includes a
      1D-TE-converges-to-grcwa(0.93333) regression lock, and 1D TE finally agrees with
      the (already-correct) 2D solver in [8a]. Didn't affect the metalens design path
      (pure 2D), but closes a credibility gap a CREOL tester sweeping M would hit.
- [x] Convergence studies documented (metric vs harmonics M, vs sampling) — `celeris
      validate` §1 tabulates phase/transmittance/energy-conservation vs M.
- [ ] Auto-convergence helper (suggest M for a target accuracy)
- [ ] Numerical edge cases hardened (Rayleigh anomalies, grazing orders, resonances)
- [ ] Regression/unit test suite + CI (lock the validated numbers so nothing drifts)

## 2. Design quality & breadth  (makes it GOOD, not a toy)
- [~] **Full-2π libraries** via height/thickness sweep (multi-parameter) → lift Strehl
      from ~0.55 toward diffraction-limited. Biggest quality win. DONE so far:
      `design --auto-height` (and `optimize_height_for_2pi`) sweeps the etch depth,
      measures each height's effective phase coverage (new `coverage()` = 2π minus the
      largest gap on the circle — the honest "can I hit any target phase" metric, robust
      to phase wrapping where the old max-min `phase_span` was not) and mean transmittance,
      and picks the best coverage at the highest |t| (still a single etch). Demonstrated
      in `validate` [4]. HONEST FINDING: a TiO₂ *square* pillar at Λ=0.35/λ=0.532 caps
      coverage ~325–345°, not a clean 360° — closing the last gap with high |t| needs
      shape variety (cross/H, holes) or a higher-index material (REMAINING below).
- [~] **Achromatic / broadband** design (group-delay + dispersion engineering) — the
      highest-value application. DONE (first cut): `celeris achromatic` + the engine
      `design/achromatic.{hpp,cpp}`. Builds a **dispersive library** over a fill×height
      grid (so the (phase, group-delay) plane is COVERED — a single DOF traces only a
      1-D curve and can't set both), extracts each atom's group delay = dφ/dω (least-
      squares over the band, phase unwrapped), then a **two-objective selection** matches
      BOTH the base focusing phase (mod 2π) AND the radius-dependent group delay (GD error
      mapped to the band-edge phase error it causes, so `--gd-weight ~1` balances center
      vs edge). Verification reuses the library's STORED per-atom band response (no extra
      RCWA). RESULT (D=10,f=30,20% BW): adding the GD objective flattens chromatic focal
      drift **8.0µm → 1.9µm (4.2× tighter)** while base-phase RMS stays diffraction-limited
      (1.7°→4.7°). selftest **[15]** locks it (drift 2.4→0.6µm, still focusing). HONEST
      LIMITATIONS: (1) the achromatic design uses MULTIPLE pillar heights → a multi-level/
      grayscale etch (a single GDS layer encodes only in-plane footprints, not per-site
      depth); (2) the GD span (hence achromatic band×aperture) is bounded by the library —
      small apertures here have GD coverage »1, but a large aperture needs taller/coupled
      atoms; (3) low-Fresnel focal shift moves the absolute focus below target (both
      designs equally). **SINGLE-ETCH achromatic — DONE** (`celeris achromatic --single-etch`):
      the multi-height limitation (1) is the real fab blocker (grayscale lithography), so
      the fabricable variant spans the (phase, group-delay) plane with **SHAPE variety at
      ONE height** (square/circle/cross/ring × fill) instead of by depth — new
      `build_single_etch_library` + `MetaAtomSpec` + a spec-driven `build_dispersive_library_
      from_specs` (the fill×height builder is now a thin wrapper over it). RESULT (D=10,f=30,
      20% BW, ONE 1.1µm etch): drift **4.0µm → 1.27µm (3.2× tighter)**, GD RMS 1.88→0.79 fs,
      mean |t| 0.75 — nearly the multi-height 4.2× but fabricable in a single standard etch
      (the design leans on cross/ring shapes). HONEST: the etch DEPTH is the key knob — taller
      pillars accumulate more group delay (h=0.6→1.3×, h=1.1→3.2×); the GD span (hence
      achromatic aperture×BW) is still smaller than the multi-height grid. selftest **[15b]**
      locks it (genuinely single-height + GD-RMS reduced; gated on GD RMS not the noisy low-
      Fresnel drift). **Exposed in the Python bindings** — `build_dispersive_library`,
      `build_single_etch_library`, `build_dispersive_library_from_specs`, `design_achromatic_
      metalens`, `verify_achromatic_focus`, `to_metalens_design` + the Dispersive{Atom,Library}
      / MetaAtomSpec / AchromaticDesign / AchromaticFocalPoint classes; example
      `examples/python/06_achromatic.py` runs both library kinds (GD RMS 3.6×/2.4× flatter).
      **ACHROMATIC PANCHARATNAM–BERRY — DONE** (`celeris pbachromatic`, the modern recipe):
      the geometric (PB) phase sets the base profile EXACTLY by rotation while a DISPERSIVE
      birefringent atom (picked per radius from a fill_x×fill_y grid at ONE etch depth)
      supplies the group delay. This DECOUPLES the two objectives — unlike the propagation-
      phase achromat (which leaves a base-phase residual because each atom must hit both
      phase AND group delay), here the rotation hits the base phase for whatever atom is
      placed, so the atom is chosen PURELY for group delay (one objective) and the base-
      phase RMS is ~0 by construction. New module `design/pb_achromatic.{hpp,cpp}`
      (`DispersivePbAtom`/`DispersivePbLibrary`, `build_dispersive_pb_library` = two RCWA
      solves/atom/λ → spin-flip a_cross=(t_x−t_y)/2 + its group delay, FILTERS out non-
      birefringent near-square atoms whose a_cross≈0 gives a garbage group delay,
      `design_pb_achromatic_metalens`, `verify_pb_achromatic_focus`) + CLI + selftest [15c]
      + a per-site rotated-rectangle GDS writer `write_pb_rect_gds` (single mask layer).
      RESULT (D=10,f=30,20% BW, ONE 1.10µm etch): base-phase RMS ~5e-15° (exact), GD RMS
      1.12→0.18 fs, chromatic focal drift 4.64→1.90 µm (2.4× tighter), conversion cap 0.52.
      **Exposed in the Python bindings** — `build_dispersive_pb_library`,
      `design_pb_achromatic_metalens`, `verify_pb_achromatic_focus` + the Dispersive
      Pb{Atom,Library} / PbAchromaticDesign classes (rotation_deg / fill_x_map / fill_y_map
      → numpy 2D); example `examples/python/07_achromatic_pb.py` (GD RMS 6.4× flatter,
      base-phase RMS ~5e-15°). **Exposed in the GUI** — an **Achromatic** panel
      (`gui/`, `run_achromatic` worker) builds a fill×height dispersive library
      (10×5, 5-band) and designs standard (gd_weight=0) vs achromatic from the SAME
      library, showing the headline chromatic focal drift (e.g. 10.2→2.8 µm, 3.6×
      flatter at f=50/D=20), group-delay RMS, GD coverage, and an overlaid focus-vs-λ
      plot (standard orange vs achromatic green); writes GDS via `to_metalens_design`
      and honestly flags the multi-height (grayscale-etch) case. The panel now
      has a **Library** selector exposing ALL THREE paths (fill×height /
      single-etch `build_single_etch_library` / Pancharatnam-Berry
      `build_dispersive_pb_library`): `run_achromatic` branches on the kind, the
      PB path writes a rotated-rectangle GDS (`write_pb_rect_gds`), and the
      single-etch/PB paths add an Etch-height slider (the one-etch GD knob) +
      kind-aware labels. REMAINING: a published achromatic device reproduction
      (K2016 reproduced above is monochromatic).
- [~] More meta-atom shapes: circular, elliptical, cross/H, ring/hole. DONE: the 2D
      solver now handles arbitrary two-material shapes via a sampled-grid Fourier
      factorization (`MetaShape` + `RectCell2D::rasterize_eps` + `shape_operators` in
      rcwa2d.cpp); CLI `design --shape circle|cross|ring [--shape-param]`. Uses the
      stable LAURENT factorization (the same scheme grcwa uses for grid layers): the
      cross matches grcwa to ~1e-3 and energy is conserved. HONEST FINDINGS: (1) the
      faster *directional inverse-rule* factorization is numerically UNSTABLE for
      non-separable shapes (diverges above ~M=10), so shapes use Laurent (slower
      convergence, needs higher M); Li's normal-vector FFF would recover speed (TODO).
      (2) Curved shapes (circle/ellipse) converge slowly in BOTH CELERIS and grcwa
      (staircase of a curved interface) — confirmed independently. (3) At Λ=0.35/λ=0.532
      a CROSS clears the 330° coverage target at ~half the pillar height of a square
      (h 0.56 vs 0.94, equal Strehl); circle ≈ square. REMAINING: shape-aware (polygon)
      GDS export — today the GDS writes square footprints regardless of `--shape`.
- [x] **Rotated pillars → Pancharatnam–Berry (geometric) phase** for circular
      polarization optics — DONE. `design --shape`/birefringence gave the building
      blocks; added in-plane atom **rotation** to the 2D solver (`RectCell2D::rotation_rad`;
      a rotated rectangle is non-separable so it takes the sampled-grid Laurent path),
      the full **2×2 Jones matrix** (`solve_jones`, two RCWA solves), an **HWP-atom
      finder** (max spin-flip conversion `|t_x−t_y|²/4`), and a PB design path
      (`design_pb_metalens` + CLI `pbdesign`): one fixed half-wave-plate atom rotated
      per site stamps the focusing phase as 2θ on the cross-circular output. Geometric
      phase is EXACT (design phase error ~1e-15°) and amplitude is uniform → inherently
      diffraction-limited, capped only by conversion efficiency (~0.64 for a TiO₂ square
      atom). RCWA-VERIFIED: solving the rotated atom confirms the spin-flip phase tracks
      −2θ to ~5° RMS (D=20, M=6). Writes a **rotated-pillar GDS** (`write_pb_gds`, true
      rotated polygons — closes part of the shape-aware GDS gap). HONEST FINDING: the
      ~5–7° residual is a method/discretization mismatch (θ=0 uses the analytic separable
      path; rotated angles use grid-Laurent at finite resolution) — shrinks with grid/M.
      Vortex/OAM & arbitrary-profile PB are DONE (any φ(r) flows through the same
      rotation map). **Achromatic PB — DONE** (`celeris pbachromatic`, the modern
      recipe; see the §2 achromatic item).
- [x] Arbitrary phase profiles: beam deflector, axicon, vortex/OAM, freeform wavefront,
      hologram/CGH. DONE. The profile machinery was lifted into a neutral `PhaseProfile`
      (focusing | vortex/OAM | deflector | axicon | **freeform**) + `phase_profile_value()`
      in `design/phase_profile.hpp` — the rotation map θ=−handedness·φ/2 (PB path) and the
      library lookup (propagation path) are both profile-agnostic, so the SAME profiles now
      drive BOTH design paths. CLI: `pbdesign --profile …` (geometric phase, exact) AND
      `design --profile …` (propagation phase) share `--charge`, `--deflect-deg`,
      `--deflect-azimuth`, `--axicon-deg`, and the **freeform** map via
      `--freeform-file <grid.txt> --freeform-extent <µm>` (a whitespace phase grid → an
      arbitrary loaded CGH/hologram, bilinearly sampled) + `--recon-z`. Per-profile optical
      proofs are now a SHARED helper (`profile_optical_proof`): focusing → on-axis z-peak;
      vortex → focal-plane donut (null + ring); deflector 10° → beam lands at ~10°; axicon →
      extended on-axis line focus; freeform → reconstruction at z. selftest [14] locks the
      freeform reproduction (a linear ramp is bilinear-exact, max |Δφ|~1e-13; focusing exact
      at nodes). HONEST: the propagation path carries the library's finite phase error
      (~34° for a deflector) + non-uniform |t|; the PB path stamps the same profiles exactly
      (RMS ~1e-14°). REMAINING (lower value): same profiles in the Python bindings + GUI.
- [ ] Full Jones-matrix metasurfaces (arbitrary polarization transforms)
- [~] Wide-FOV / off-axis designs (quadratic phase, doublets, aplanatic). DONE
      (first cut): a **quadratic-phase** wide-FOV lens. New `PhaseProfileKind::Quadratic`
      (phi = -k*r^2/(2f)) on the propagation-phase path (`design --profile quadratic`),
      and `celeris widefov` — designs a hyperbolic AND a quadratic lens for the same f/D
      from the SAME library and sweeps the field angle through an aperture STOP placed in
      front of the metasurface. KEY PHYSICS: under tilt the quadratic phase is the same
      parabola simply RECENTERED at x0 = f*sin(theta), so its focus just translates and
      stays sharp; the hyperbolic lens develops coma. The wide-FOV benefit appears ONLY
      with the offset stop (default: front focal plane, distance = f): each field angle
      then samples a DECENTERED low-NA patch of the larger lens — coma-free for the
      quadratic (vertex-centered parabola), coma-laden for the hyperbolic. New
      `analyze_wide_fov` (offset-stop, decentered-patch propagation) alongside the
      existing `analyze_field_of_view` (stop-at-lens; the right model for a normal lens's
      usable field, where quadratic shows NO benefit — an honest finding). RESULT
      (f=30, lens D=70, stop D=20 at 30µm): hyperbolic FOV-to-Strehl-0.8 ±15.4° vs
      quadratic ≥±30° (sweep cap); to-0.5 ±21.5° vs ≥±30°. HONEST COSTS: resolution is
      set by the (small) stop ~lambda*f/D_stop not the full lens, and the focal surface
      is curved (Petzval). selftest **[17]** locks it (at 30°: hyperbolic 0.27 vs
      quadratic 0.96). Exposed in the Python bindings (PhaseProfileKind.Quadratic +
      PhaseProfile.quadratic factory). REMAINING: doublet (2-surface corrector + lens,
      Arbabi 2016) and aplanatic designs; a GUI panel.
- [ ] Topology / freeform meta-atom inverse design (research-grade)
- [ ] Fab constraints baked into design (min feature, min gap, etch bias)

## 3. Productization  (makes it USABLE / SHIPPABLE)
- [x] **Python bindings (pybind11)** over celeris_core — scripting + numpy/matplotlib;
      how this market actually works. DONE: `_celeris` extension (pybind11 v3.0.4,
      `bindings/python/celeris_module.cpp`) exposes the full workflow — Material
      (+ builtins), 1D/2D RCWA solvers, UnitCellLibrary / build_unit_cell_library /
      optimize_height_for_2pi, design_metalens, analyze_focus, compute_psf (maps
      come back as numpy arrays). Opt-in `-DCELERIS_BUILD_PYTHON=ON` stages an
      importable package at `build/python/celeris` (MKL/CUDA DLLs bundled beside the
      .pyd); `pip install .` works via scikit-build-core (`pyproject.toml`). Four
      runnable examples in `examples/python/` (focusing lens, RCWA convergence,
      real-material dispersion, full-2π height sweep). Verified: import + end-to-end
      design reproduces the engine's numbers (diffraction-limited FWHM, energy=1).
- [ ] **Installer** (WiX/NSIS) + code signing → a thing people download, not a dev exe
- [ ] **One build, auto CPU/GPU** at runtime (bundle CUDA runtime, fall back to CPU if
      no device) — today GPU is a separate build
- [ ] Example projects / templates gallery (focusing lens, deflector, polarization
      splitter, …) so a new user starts from something
- [ ] User manual + tutorials + theory notes; API docs (Doxygen)
- [ ] Report export to PDF/HTML (today: txt + images)

## 4. Analysis & materials depth  (Zemax parity, mostly incremental)
- [ ] Field-resolved grids (PSF/MTF/Strehl across the full field, not just a row)
- [ ] Efficiency breakdown (per diffraction order, reflection, absorption budget)
- [ ] Polarization analysis maps (Stokes/Mueller, extinction-ratio map)
- [ ] Stray light / higher-order ghosts
- [~] Bundled real n,k library — TiO₂ (Siefke 2016, ALD amorphous) shipped in
      `data/TiO2_Siefke.csv`; still need a-Si, c-Si, GaN, SiN, Al₂O₃, Au/Ag/Al
- [ ] Material editor in GUI (fit Sellmeier to data, anisotropic, thermo-optic)
- [ ] Fabrication-aware analysis (sidewall angle, corner rounding → performance)

## 5. GUI / UX polish  (makes it feel like a product)
- [x] Dark mode + themed native title bar (immersive dark), commit in title, app icon
- [x] Internal performance monitor (FPS + frame-time graph, process CPU%/memory, device)
- [x] GUI modularized (split the monolith into celeris::gui: theme/textures/app_state/workers)
- [ ] ImPlot interactive charts (zoom/pan/hover) — replace static plots
- [ ] Parametric sweep / optimization UI (vary param → live metric plot), merit editor
- [ ] Undo/redo, units management, input validation, tooltips/contextual help
- [ ] 3D view (pillar array, field volume)
- [ ] Compare-designs side by side; recent files; drag-drop; layout persistence; DPI scaling

## 6. I/O & formats  (interoperability)
- [ ] GDS layers/datatypes + cell hierarchy (multi-material / multi-height masks)
- [ ] More export formats: DXF, OASIS, CIF; optional STL/3D
- [ ] Import-and-analyze (GDS/phase-map → performance), round-trip with foundries
- [ ] Process/stack definition files

## 7. Performance & scale
- [x] **AVX2/FMA build (default)** — `CELERIS_AVX2=ON` adds `/arch:AVX2`; MSVC's SSE2
      x64 baseline left Eigen's dense complex matmuls badly under-vectorized. ~2.7x
      free speedup (a 2D solve at M=10 went 44s -> 16s), no dependency.
- [x] **Intel MKL turbo (opt-in)** — `CELERIS_USE_MKL=ON` routes ALL of Eigen's dense
      BLAS/LAPACK through MKL (`EIGEN_USE_MKL_ALL`) *plus* the zgeev shim. KEY FINDING:
      the eigensolve is only ~15% of a solve — the real cost is the dense matmuls /
      inverses in operator assembly + the S-matrix recursion, which MKL multithreads.
      On top of AVX2 this is up to ~11x over the SSE2 baseline at high M (M=10: 44s->3.9s,
      M=12: ~120s->10.7s). Needs MKL on PATH (`pip install mkl-devel mkl-include`);
      trades the single-self-contained-binary property, hence opt-in. (The old "only
      MKL beats Eigen *for zgeev*" note was right about zgeev but missed that zgeev
      wasn't the bottleneck.)
- [ ] GPU-batch the whole library/wavelength sweep in one launch
- [ ] Out-of-core / mm-scale apertures (millions of pillars), memory streaming
- [ ] Library/result caching; incremental recompute (only what changed)
- [ ] Cancellable long runs; multi-GPU (later)

## 8. Ops / platform
- [ ] Linux build + CI (engine is already portable); macOS optional
- [ ] Versioning, release notes, changelog; crash reporting (opt-in); auto-update
- [ ] Licensing/activation if commercial

---

## ⭐ LAB-READY milestone (the "yo try my shit" line)
You do NOT need all of the above to share with labmates. Lab-ready =
**trusted + good designs + installable + a couple examples**:
1. **#1 validation case** (one published device reproduced) — they'll believe it
2. **Full-2π library** (#2) — so their designs are actually sharp, not capped at 0.55
3. **Installer or Python bindings** (#3, pick one) — so they can run it without your toolchain
4. **3–4 example projects + a short quickstart** (#3) — so they're not staring at a blank app
5. Real n,k material library (#4) — so they design in their actual materials

Hit those five and it's lab-shareable. Everything else is the road from lab tool
to sellable product.

## Rough sequencing
1. **Validation case** (#1) — cheap, software-only, unlocks trust
2. **Full-2π libraries** (#2) — biggest quality jump
3. **Python bindings + examples** (#3) — adoption
4. **Achromatic + PB-phase** (#2) — the differentiated, high-value designs
5. **Installer + auto CPU/GPU + docs** (#3/#8) — shippable
6. Breadth (#4/#5/#6) and scale (#7) as demand dictates

idk some notes by gemini

Conversation with Gemini

what does this look like to you

Based on the file image_83d105.png, this looks like a specialized software application interface designed for metalens simulation and design.

It appears to be an optical engineering tool used to model, optimize, and analyze metasurfaces—specifically looking at how subwavelength nano-pillars interact with light.

Here is a breakdown of what the different panels show:
1. Design Parameters (Left Panel)

This panel lets you configure the physical and optical properties of the metalens:

    Dimensions: It specifies a focal length of 50.0 µm, an aperture of 20.0 µm, a wavelength of 0.532 µm (green light), a period of 0.200 µm, and a pillar height of 0.47 µm.

    Materials: The pillar material is set to Titanium Dioxide (TiO2​, index n≈2.40) on an N-BK7 glass substrate.

    Simulation & Output: There are options for running the design (likely using RCWA, or Rigorous Coupled-Wave Analysis, as indicated by the "RCWA harmonics" field) and exporting the layout to GDSII format (metalens.gds), which is standard for semiconductor fabrication.

2. Design Summary (Middle Top Panel)

This gives a snapshot of the resulting performance characteristics:

    Array Size: 100×100 pillars (10,000 total).

    Numerical Aperture (NA): 0.20.

    Phase Coverage: 332 degrees.

    Errors & Efficiency: Shows an RMS phase error of 9.4 degrees and a mean transmission efficiency of 76.3%.

3. Chromatic Performance Graph (Right Panel)

    This plot tracks Focal length (um) vs wavelength.

    The downward-sloping blue line clearly illustrates chromatic aberration, demonstrating how the focal length decreases as the wavelength increases—a classic characteristic of diffractive optical elements and simple metalenses unless they are explicitly designed to be achromatic.

4. Status Log (Bottom Panel)

    The log shows that a simulation or calculation just finished ("Done.").

    It notes that the propagation was hardware-accelerated using a local GPU: Propagation: GPU (NVIDIA GeForce RTX 4070 Laptop GPU).

what are 10 things its missing to be a proper metalens design program

While the interface in the screenshot covers the absolute basics—defining a single unit cell, calculating phase coverage, and running a simple 2D or 3D propagation—it is quite barebones compared to commercial-grade optical design suites (like Ansys Lumerical, Synopsys LightTools/RSoft, or Comsol).

To be a robust, industry-standard metalens design program, it is missing several critical features:
1. Robust Achromatic Design Suite

The graph on the right shows a massive drop in focal length as wavelength increases. A proper program needs tools to design Achromatic Metalenses (AMLs). This requires solving for both phase and group delay / group delay dispersion (GDD) across a bandwidth, rather than just matching a single target phase at one wavelength.
2. Complex Unit-Cell Geometry Library

The software seems limited to simple "pillars" (likely cylinders or squares). Modern metasurfaces use intricate geometries to manipulate light properties. A proper tool needs a built-in library for:

    Holes, rings, and cross-shaped resonators.

    "Molecules" (multipillars per unit cell).

    Freeform topology optimization shapes.

3. Full Vectorial polarization Control

The interface has a "Polarization" tab, but a proper program needs deep support for birefringent unit cells (e.g., rectangular or elliptical pillars). It should automatically calculate and optimize the full Jones Matrix or Mueller Matrix for every cell to handle arbitrary polarizations, asymmetric focusing, or generating orbital angular momentum (OAM) beams.
4. Advanced Optimization Algorithms

"Optimize Design (period x height)" is incredibly basic. A professional program requires:

    Adjoint-based inverse design (Topology Optimization) to let the computer figure out the exact nanostructure shapes for a specific target efficiency.

    Global optimization algorithms (Genetic Algorithms, Particle Swarm) for unit-cell multiplexing.

5. Multi-Angle / Wide Field-of-View (FoV) Engineering

Metalenses notoriously suffer from severe off-axis aberrations (coma and astigmatism). A proper software needs tools to analyze and optimize the phase profile for oblique incident angles, or to design multi-layer metalens doublets/triplets to correct for wide FoV.
6. Rigorous Material Dispersion Engine

The user input requires a static refractive index (index n = 2.40). Real materials experience dispersion (index changes drastically with wavelength) and have a complex refractive index (n~=n+ik) where k represents material absorption/loss. The program needs a built-in material database (e.g., Palik data) that accurately accounts for real-world material losses across the spectrum.
7. Full-Wave Solvers (FDTD/FEM) Co-Simulation

While RCWA (Rigorous Coupled-Wave Analysis) is great for fast, periodic unit-cell libraries, it can struggle with large-angle scattering or edge-effects at the boundaries of the metalens where the period assumption breaks down. A complete suite needs integrated FDTD (Finite-Difference Time-Domain) or FEM (Finite Element Method) solvers to simulate the entire macro-scale metalens accurately.
8. Ray-Tracing Export & Macro-System Integration

A metalens is rarely used entirely on its own; it's usually part of a larger camera or sensor system. A proper program needs to seamlessly export the metalens as a phase mask or holographic surface into ray-tracing software like Zemax OpticStudio or Code V for system-level camera design.
9. Fabrication Constraint Checking (DFM)

Nano-lithography (E-beam or DUV photolithography) has physical limits. A proper design tool needs a Design for Manufacturing (DFM) checker to flag features that can't actually be built, such as:

    Aspect ratios that are too high (pillars will collapse).

    Gaps between pillars that are smaller than the minimum lithographic resolution (side-wall merging).

10. Thermal and Mechanical Stress Analysis (STOP)

Metasurfaces are often deployed in rugged environments (automotive sensors, satellites, phones). Professional software includes Structural, Thermal, and Optical Performance (STOP) analysis to model how thermal expansion changes the pillar spacing and refractive index, and how mechanical stress warps the substrate, degrading focus.

Is there a specific feature or workflow from this list you are trying to implement or evaluate for your current project?


Gemini is AI and can make mistakes.
