# P4 — Certified surrogate models for the RCWA forward map

**Status:** proposal / not started. **Sequencing: strictly AFTER the P1 JOSS submission is
out the door, and after (or alongside the tail of) the §3d adjoint.** Nothing here touches
`paper/` or the in-flight submission. If you are reading this while the JOSS submission is
still unsent, close this file and go do that instead.

---

## 1. Why this exists

The novelty risk on P1 is already known and recorded: *"another RCWA solver."* The §3d
adjoint mitigates it, but be honest about what the adjoint actually is — **catch-up**. It
closes the one column CELERIS currently loses to fmmax / TORCWA / meent (autodiff). Landing
it means CELERIS finally *matches* the field on differentiability. It does not put CELERIS
anywhere the field is not.

A **learned surrogate for the forward map** is the first thing in this project's plan that
the autodiff-native tools do not ship. And CELERIS is unusually well-positioned for it,
because the expensive part of surrogate work is not the model — it is a trustworthy ground-truth
generator, and this repo already has one:

- a solver validated three ways (CELERIS = grcwa = S4 = 0.93333, TE residual 7e-8),
- 55 locked-tolerance assertions gating every release,
- a meta-atom library sweep that *already* generates exactly the dataset a surrogate needs,
- GPU throughput for bulk data generation,
- and a documented culture of retracting overstated numbers (the 600–700x → 1–6x fix).

That last one is the actual asset. See §3.

## 2. The two tiers (do Tier 1 first; Tier 2 may never happen)

**Tier 1 — the meta-atom library map.** Learn
`(λ, pillar width/radius, height, material) → complex transmission t`.
Low-dimensional input, smooth output, and the library sweep is the dominant cost in a design
run. Tractable in weeks, not months. This is the whole of the near-term plan.

**Tier 2 — freeform / patterned unit cells.** Input becomes a *pattern*, not a handful of
scalars, so this needs operator learning (FNO / DeepONet class) and a much larger dataset.
Only consider it if Tier 1 clears its gates and the arbitrary-polygon FFF work (§3-adjacent)
has landed. Do not start here.

## 3. The differentiator — and an honest novelty audit

**Be clear-eyed: "neural network surrogate for metasurface design" is a crowded literature.**
There are many papers training nets on meta-atom sweeps. Shipping one more is not a
contribution and would deserve the same desk-reject the P1 framing was built to avoid.

What is *not* crowded, and what this project is uniquely equipped to do:

> **A surrogate with a distribution-free error guarantee against a machine-precision
> reference, wired to an automatic fallback to the exact solver when the guarantee does not
> hold — so the pipeline never silently returns a wrong device.**

Concretely: use **conformal prediction** to calibrate a per-query error bound on the
surrogate's complex transmission. That yields a statement of the form *"with ≥1−α coverage,
|t_surrogate − t_exact| ≤ ε,"* with no distributional assumptions. Queries whose predicted
bound exceeds the design tolerance get routed to the exact RCWA solve. The surrogate becomes
an *accelerator with a certificate*, not a replacement you have to trust on faith.

This is the same philosophy as the locked-tolerance selftest, applied to a statistical
component — which is why it fits this repo and not a generic ML-photonics paper. It is also
the honest answer to the question every referee will ask: *"how do I know the fast path
didn't lie to me?"* Nobody in the surrogate-metasurface literature answers that properly.

## 4. Why the adjoint comes first (real dependency, not just ordering)

- The exact adjoint gives **gradients of the true forward map**, which are (a) a far stronger
  training/validation signal than function values alone and (b) the reference for checking
  that surrogate gradients are usable for design, not just its values.
- It enables the payoff architecture: **surrogate for global exploration → exact adjoint for
  local refinement.** Multi-fidelity / trust-region optimization, with the conformal gate
  deciding when to spend an exact solve. The surrogate explores cheaply; the adjoint
  guarantees the final device is real.
- Without the adjoint, the surrogate is a lookup-table accelerator. With it, it is a design
  method.

## 5. Baselines that must be beaten (the honesty gate)

Most ML-for-photonics work skips this. Do not.

1. **Cubic spline / RBF interpolation over the existing library grid.** The Tier-1 map is
   smooth and low-dimensional — classical interpolation is the real baseline, not "no
   surrogate." **Kill criterion: if a spline matches the net inside tolerance at lower total
   cost, Tier 1 has no reason to exist.** Report that result honestly and pivot to Tier 2 or
   drop the whole line.
2. **The exact solver with MKL.** Report the **break-even point** including data generation
   and training cost. For a single one-off design the surrogate is a net loss and the paper
   must say so; it wins across repeated design campaigns / large spectral sweeps. That
   break-even curve is itself a publishable, useful result.
3. **Ablate the certificate.** Show what fraction of queries the conformal gate routes to the
   exact solver, and what the achieved coverage is versus nominal 1−α.

## 6. Build plan

- [ ] **Data.** Script a reproducible library-sweep dump (params → complex t) from the
      existing sweep machinery. Fixed seed, versioned, hash-stamped. Held-out test set drawn
      from a *different* region of parameter space than the training grid, not a random split —
      random splits over a dense grid massively overstate accuracy.
- [ ] **Baseline first.** Spline/RBF interpolator + error table. Do this **before** any
      network, so the kill criterion is evaluated honestly rather than retroactively.
- [ ] **Model.** Small MLP with Fourier features on `(λ, geometry)`; predict real/imag parts
      (not amplitude/phase — phase wrapping will wreck the loss). Keep it small enough to
      ship; this is not a scale exercise.
- [ ] **Conformal calibration.** Split-conformal on a held-out calibration set; report
      empirical coverage vs nominal across parameter regions. Flag where coverage degrades.
- [ ] **Fallback gate.** `solve_or_surrogate(..., tol)` in the design path: surrogate when the
      certified bound clears `tol`, exact solve otherwise. Log the routing fraction.
- [ ] **Locked test cases.** Add selftest assertions in the repo idiom: (a) surrogate error on
      the held-out set below a locked bound; (b) empirical conformal coverage ≥ nominal on a
      fixed seed; (c) an end-to-end design built through the gated path matching the
      exact-path device inside tolerance. Physics gate stays untouched.
- [ ] **Hybrid loop** (needs §3d): surrogate-explored → adjoint-refined design, versus
      adjoint-only from a cold start. Compare exact-solve count to reach equal device quality.
      This is the headline number.

## 7. Paper placement

| # | Paper | Venue | Content | Gate |
|---|-------|-------|---------|------|
| **P4** | Certified surrogate + hybrid design | CPC / Optics Express / ACS Photonics | "A calibrated surrogate with distribution-free error guarantees and exact-solver fallback, and the multi-fidelity design loop it enables." | Tier 1 clears §5 kill criterion; adjoint landed for the hybrid result |

Distinct content from P1 (the pipeline) and P2 (the adjoint method) — the existing
one-paper-per-distinct-content rule holds. If Tier 1 dies on the spline baseline, **there is
no P4**; say so and move on, exactly as §3d is already gated.

## 8. Why this is worth doing beyond the paper

CELERIS is currently, on paper, an optics tool. This line converts it into the substrate for
computational/statistical-physics work: a validated forward model, differentiated, accelerated,
with calibrated inference on top. That is a portable methodological identity rather than a
subfield commitment — the same machinery transfers to any domain with an expensive forward
model. It is also the one addition here that is simultaneously a novelty fix for the project
and direct evidence for a compute/ML-physics research trajectory.

## 9. Scope guard

- Do **not** add ML claims, hedges, or "future work" edits to the submitted P1 paper. It is
  submission-ready; leave it alone.
- Torch/JAX must be an **optional** Python-side dependency. The C++ engine, the CPU-only
  build, and CI stay ML-free — reviewers clone on Linux with no GPU and no ML stack.
- The exact solver remains the source of truth. The surrogate is never the default path
  without a certificate.
