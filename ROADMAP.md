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
- [ ] **Validation vs a published, fabricated metalens** — reproduce a paper's
      measured efficiency/Strehl/FWHM (e.g. Khorasaninejad 2016 TiO₂). THE #1 item.
- [ ] Convergence studies documented (metric vs harmonics M, vs sampling) so users
      know how to trust/tune accuracy
- [ ] Cross-check vs an external solver (S4 / grcwa / Lumerical) on 1–2 benchmarks
- [ ] Auto-convergence helper (suggest M for a target accuracy)
- [ ] Numerical edge cases hardened (Rayleigh anomalies, grazing orders, resonances)
- [ ] Regression/unit test suite + CI (lock the validated numbers so nothing drifts)

## 2. Design quality & breadth  (makes it GOOD, not a toy)
- [ ] **Full-2π libraries** via height/thickness sweep (multi-parameter) → lift Strehl
      from ~0.55 toward diffraction-limited. Biggest quality win.
- [ ] **Achromatic / broadband** design (group-delay + dispersion engineering) — the
      highest-value application
- [ ] More meta-atom shapes: circular, elliptical, cross/H (needs non-separable
      Fourier or Li factorization)
- [ ] **Rotated pillars → Pancharatnam–Berry (geometric) phase** for circular
      polarization optics
- [ ] Arbitrary phase profiles: beam deflector, axicon, vortex/OAM, freeform wavefront,
      hologram/CGH
- [ ] Full Jones-matrix metasurfaces (arbitrary polarization transforms)
- [ ] Wide-FOV / off-axis designs (quadratic phase, doublets, aplanatic)
- [ ] Topology / freeform meta-atom inverse design (research-grade)
- [ ] Fab constraints baked into design (min feature, min gap, etch bias)

## 3. Productization  (makes it USABLE / SHIPPABLE)
- [ ] **Python bindings (pybind11)** over celeris_core — scripting + numpy/matplotlib;
      how this market actually works
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
- [ ] Bundled real n,k library (TiO₂, a-Si, c-Si, GaN, SiN, Al₂O₃, Au/Ag/Al)
- [ ] Material editor in GUI (fit Sellmeier to data, anisotropic, thermo-optic)
- [ ] Fabrication-aware analysis (sidewall angle, corner rounding → performance)

## 5. GUI / UX polish  (makes it feel like a product)
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
- [ ] GPU-batch the whole library/wavelength sweep in one launch
- [ ] Out-of-core / mm-scale apertures (millions of pillars), memory streaming
- [ ] Library/result caching; incremental recompute (only what changed)
- [ ] Cancellable long runs; multi-GPU (later)
- [ ] (Optional) Intel MKL eigensolve plug-in — only real lever left for the library build

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
