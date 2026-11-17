# R1 metric: fraction of the achromatic limit reached by a single-etch library

Issue #23. This note fixes the definitions every R1 sweep reports. The formulas are
implemented in [`metric.py`](metric.py); `python3 research/r1/metric.py` runs the
numeric self-checks and prints the table in §6.

Units: µm, fs, rad/fs; c = 0.299792458 µm/fs.

## 1. What CELERIS measures

For each meta-atom, CELERIS fits a group delay τ = dφ/dω across the band. This is the
least-squares slope of the unwrapped transmission phase: arg a_cross for
Pancharatnam–Berry (PB) atoms, arg t for propagation-phase atoms. A library's
**delay span** is

  Δτ_lib = max τ − min τ over the atoms kept by the quality filters (#20).

Δτ_lib is the only library property the limit formulas need.

## 2. The lens needs a delay span; diameter follows from it

An achromatic lens of focal length f and radius R must delay its centre ray relative
to its edge ray by the path difference over c (Presutti & Monticone, Eq. 4):

  ΔT(R, f) = (√(R² + f²) − f) / c.   (1)

CELERIS computes this as `required_gd_span_fs`. Setting ΔT = Δτ_lib and inverting gives
the largest achromatic diameter the library supports:

  D_lib = 2f √((1 + cΔτ_lib / f)² − 1).   (2)

**Eq. 2 is Buisine et al. Eq. 15** with Δβ₁ = Δτ_lib.

### Buisine Eq. 7 is the same identity

Buisine Eq. 7 reads D = 2h(Δn − Δn₀)·NA / (1 − √(1 − NA²)). Write x = cΔτ and
NA = R/√(R² + f²) (lens in air). From Eq. 1, √(R² + f²) − f = x, so
1 − √(1 − NA²) = x/√(R² + f²) and NA/(1 − √(1 − NA²)) = R/x. Then

  2x · R/x = 2R = D.

So Eq. 7 and Eq. 15 are the same geometric relation:
- Eq. 7 is parametrised by NA, Eq. 15 by f;
- the library enters Eq. 7 as the optical-path spread x = h(Δn − Δn₀) and Eq. 15 as cΔβ₁;
- the identity holds exactly when h(Δn − Δn₀) = cΔτ.

The paper does not reconcile the two equations, and `metric.py` checks the equality
numerically. **Neither equation is a bound by itself.** Each converts a delay span into
a diameter. The physics is in how large Δτ can be.

Open: Buisine define Δn as the "maximal variation of effective refractive index over
the bandwidth and achromatic subset" and give no formula for Δn₀. We have not read
their supplement verbatim, so we do not claim that their Δn equals cΔτ/h exactly. R1
uses Δτ, which is what CELERIS measures and what Eqs. 1–2 need.

## 3. The bound on Δτ (Presutti & Monticone 2020)

Any linear time-invariant structure obeys a delay–bandwidth product bound
(P&M Eq. 3):

  ΔT · Δω ≤ κ.

For a lossless dielectric waveguide-like structure of thickness L, Tucker's κ is
(P&M Eq. 8):

  κ = 2π (L/λ_c)(n_max − n_min) = ω_c L Δn / c.

With L = h (the single etch depth), n_max the pillar index and n_min the background
index, the bound on the library's delay span is:

  Δτ_PM = (h Δn / c) · (ω_c / Δω).   (3)

The first factor is the **non-dispersive delay line**:

  Δτ_DL = h(n_max − n_min)/c.   (4)

This is the delay difference between crossing h entirely in pillar and entirely in
background. The second factor is the inverse fractional bandwidth. So for a band with
fractional bandwidth Δω/ω_c:

  Δτ_PM = Δτ_DL / (Δω/ω_c).

A plain delay line, whose delay does not depend on frequency, gets only the fraction
Δω/ω_c of the bound. The rest is available only to frequency-dependent (slow-light or
resonant) delay, and that costs efficiency. This is the trade-off R1 maps.

Alternatives in P&M:
- the single-resonator κ = 2 (Eq. 7), for thin resonant atoms;
- Miller's κ = (π/√3)(L/λ_c)η_max (Eq. 9), which strictly needs L ≫ λ_c.

R1 uses Tucker's κ as the primary bound, because its assumption is a guided-wave
structure with Δω/ω_c < 1, which matches a pillar. Miller's bound is a cross-check.

## 4. Definitions R1 reports

| symbol | definition | meaning |
|---|---|---|
| Δτ_lib | Eq. §1 | delay span the filtered library supplies |
| η_DL | Δτ_lib / Δτ_DL | fraction of the delay-line scale; **can exceed 1** (guided-mode dispersion) |
| η_PM | Δτ_lib / Δτ_PM = η_DL · Δω/ω_c | fraction of the Presutti–Monticone (Tucker) bound; **must be ≤ 1** |
| D_lib | Eq. 2 with Δτ_lib | largest achromatic diameter at a given f |
| D_PM | Eq. 2 with Δτ_PM | the bound's diameter at the same f |

Conventions, fixed so results are comparable:

- **ω_c = (ω_min + ω_max)/2**, the band centre in frequency. It is not the ω of the
  wavelength midpoint. For 470–670 nm, λ_c = 552.5 nm and Δω/ω_c = 0.3509.
- **n_max = n_pillar(λ_c)** and **n_min = n_background(λ_c)**, the phase indices as in
  Tucker's formula. For dispersive materials, report a group-index variant later as a
  sensitivity check; it is not the headline.
- **Lens in air** (n_b = 1 in P&M's NA = n_b sin(arctan(R/f))).
- **Sanity rule: η_PM > 1 is a red flag, not a result.** It means an aliased or slipped
  group delay (#20, #21), or an atom outside Tucker's model (a resonator), and must be
  investigated before it is reported.
- **Efficiency is a second axis, not a fixed cut.** Report η_PM as a curve against the
  amplitude threshold (mean |a_cross| for PB atoms, mean |t| for propagation-phase
  atoms). Buisine keep T ≥ 0.9 and drop resonant atoms; CELERIS's default is mean
  |a_cross| ≥ 0.30. The span is set by a few atoms at the edge of the distribution, and
  at H = 1 µm the group delays cluster near 4 and 9 fs. So one cut would hide the
  trade-off.

## 5. What the bound assumes

- The bound is for a diffraction-limited lens (Strehl 1) and a time-invariant, linear,
  passive structure. P&M do not claim it is tight. They report that all published
  metalenses obey the Tucker and Miller versions.
- Eq. 1 neglects geometric aberrations and higher-order dispersion (group-delay
  dispersion, GDD) over wide bands. Buisine flag GDD as necessary for large bandwidths.
  R1's full-lens check (verify_*_focus) is locally periodic, not full-wave.

## 6. Numbers for the libraries that already exist

Band 470–670 nm. f = 64.67 µm (the Chen 2018 lens: NA 0.20, D 26.4 µm, required delay
4.45 fs). Constant-index rows use n = 2.40. The ALD row uses Siefke TiO₂ at λ_c
(n = 2.434). Spans are from PR #22 (2026-11-16).

| library | Δτ_lib (fs) | Δτ_DL | Δτ_PM | η_DL | η_PM | D_lib (µm) | D_PM (µm) |
|---|---|---|---|---|---|---|---|
| pilot n=2.4, H=1.2 µm (7 pts) | 5.97 | 5.60 | 15.97 | 1.065 | 0.374 | 30.6 | 50.7 |
| fig_achromatic n=2.4, H=0.6 µm (7 pts) | 2.58 | 2.80 | 7.99 | 0.921 | 0.323 | 20.1 | 35.5 |
| selftest [18] n=2.4, H=0.6 µm (3 pts) | 2.21 | 2.80 | 7.99 | 0.789 | 0.277 | 18.6 | 35.5 |
| `reproduce --device chen2018` (ALD), H=0.6 µm | 3.11 | 2.87 | 8.18 | 1.084 | 0.380 | 22.0 | 36.0 |

Reading:
- CELERIS's simple rectangles land at **η_DL ≈ 0.8–1.1**, about the delay-line scale,
  and **η_PM ≈ 0.28–0.38**, close to Δω/ω_c = 0.35. This is what Eq. 3 predicts for
  non-resonant atoms.
- For comparison, the published Chen 2018 device needs 4.45 fs at H = 0.6 µm, which is
  η_PM ≈ 0.54 and η_DL ≈ 1.55. This uses Chen's published lens requirement, not a
  CELERIS-computed library. Their coupled, resonant nanofins go past the delay line.
  That is the expected route above η_PM ≈ Δω/ω_c, and R1 should measure what it costs
  in efficiency.
- The 3-point row is the least reliable (#21). R1 sweeps use ≥ 7 band points and spot
  checks against 21 points.

## 7. Sources and how they were checked

- Buisine et al., "Beyond Group Delay: Fundamental Limits of Achromatic Metasurfaces",
  arXiv:2610.08047 (Oct 2026). Eqs. 7 and 15 and the library details were read verbatim
  from the HTML main text; the supplement was not.
- Presutti & Monticone, "Focusing on Bandwidth: Achromatic Metalens Limits",
  arXiv:2001.10899 / Optica 7, 624 (2020). Eqs. 3, 4, 6–9 were read verbatim (about 74%
  of the HTML), and Eqs. 11–12 were checked to follow from 6–9. Table 1 was not
  verified.
