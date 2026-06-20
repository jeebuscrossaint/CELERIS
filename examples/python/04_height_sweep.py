"""Full-2pi search: sweep pillar height to maximize phase coverage at high |t|.

A single etch depth can cap phase coverage short of 360 deg (or hit it only at a
low-transmittance height), which limits the Strehl. optimize_height_for_2pi
sweeps the height, measures each height's coverage + mean transmittance, and
picks the best. Plots coverage and transmittance vs height.

    PYTHONPATH=build/python python examples/python/04_height_sweep.py
"""
import os
import celeris as cel

WAVELENGTH, PERIOD = 0.532, 0.35
air, glass = cel.materials.air(), cel.materials.fused_silica()
tio2 = cel.Material.constant(2.40 + 0.0j, "TiO2")

res = cel.optimize_height_for_2pi(
    pillar=tio2, background=air, incident=air, substrate=glass,
    period_um=PERIOD, wavelength_um=WAVELENGTH,
    thick_lo=0.30, thick_hi=1.10, n_heights=12,
    fill_min=0.1, fill_max=0.9, fill_samples=24, M=8,
    coverage_target_deg=330.0)

print(f"best height = {res.best_thickness_um:.3f} um  "
      f"coverage = {res.coverage_deg:.0f} deg  "
      f"mean transmittance = {res.mean_transmittance:.3f}  "
      f"(target {'reached' if res.reached_target else 'NOT reached'})")
print(f"\n{'height(um)':>10}  {'coverage(deg)':>13}  {'mean T':>8}")
for e in res.sweep:
    mark = "  <-- best" if abs(e.thickness_um - res.best_thickness_um) < 1e-9 else ""
    print(f"{e.thickness_um:>10.3f}  {e.coverage_deg:>13.0f}  {e.mean_transmittance:>8.3f}{mark}")

try:
    import matplotlib.pyplot as plt
    h = [e.thickness_um for e in res.sweep]
    cov = [e.coverage_deg for e in res.sweep]
    tr = [e.mean_transmittance for e in res.sweep]
    fig, ax1 = plt.subplots(figsize=(7, 4.5))
    ax1.plot(h, cov, "o-", color="C0", label="coverage (deg)")
    ax1.axhline(res.coverage_target_deg, ls=":", color="C0", alpha=0.6)
    ax1.axvline(res.best_thickness_um, ls="--", color="0.4")
    ax1.set(xlabel="pillar height (um)", ylabel="phase coverage (deg)")
    ax2 = ax1.twinx()
    ax2.plot(h, tr, "s--", color="C3", label="mean transmittance")
    ax2.set_ylabel("mean transmittance |t|^2")
    ax1.set_title("Full-2pi height sweep — TiO2 square pillar")
    fig.tight_layout()
    out = os.path.join(os.path.dirname(__file__), "04_height_sweep.png")
    fig.savefig(out, dpi=130)
    print("\nsaved", out)
except ImportError:
    print("(install matplotlib to render the figure)")
