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
- [~] **Validation vs a published, fabricated metalens** — `celeris validate` runs the
      credibility battery on REAL tabulated TiO₂ n,k (Siefke 2016 ALD, in `data/`)
      on fused silica. DONE: a real-TiO₂ focusing lens is **diffraction-limited**
      — phase Strehl 0.97, FWHM exactly at λf/D, NA 0.20 — and the meta-atom
      transmittance now **converges** (Li factorization landed) and matches an external
      solver (grcwa). Absolute efficiency is trustworthy. REMAINING for a full
      paper reproduction: a specific published device + (for K2016) the PB-phase
      design path; and full-2π libraries to push the focusing efficiency higher.
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
- [ ] **1D TE solver convergence** — NEW finding: the 1D solver's TE transmission for a
      mid-period grating (Λ=0.3, λ=0.5) oscillates with M (0.973→0.887) instead of
      converging to the grcwa-true 0.933 (1D TM is fine). Pre-existing; exposed by the
      grcwa cross-check. Doesn't affect the metalens design path (pure 2D), but fix for
      full credibility.
- [x] Convergence studies documented (metric vs harmonics M, vs sampling) — `celeris
      validate` §1 tabulates phase/transmittance/energy-conservation vs M.
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
