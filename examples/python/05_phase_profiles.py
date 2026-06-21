"""Arbitrary phase profiles + freeform holograms, on both design paths.

The SAME PhaseProfile drives two physically different metasurfaces:
  * propagation-phase (`design_metalens`): vary the pillar SIZE to hit phi(x,y)
    via a library lookup -> a finite residual phase error;
  * geometric-phase / Pancharatnam-Berry (`design_pb_metalens`): rotate one fixed
    half-wave-plate atom by theta = -handedness*phi/2 -> phi is EXACT.

This demo designs a vortex, a deflector, and a freeform (CGH) hologram on the
propagation path, then a PB focusing lens + PB vortex (exact geometric phase),
printing the design fidelity of each.

Run (from the repo root, after building with -DCELERIS_BUILD_PYTHON=ON):
    PYTHONPATH=build/python python examples/python/05_phase_profiles.py
"""
import math
import os
import celeris as cel

# --- shared parameters ------------------------------------------------------
WAVELENGTH = 0.532   # um
PERIOD = 0.35        # um lattice pitch
HEIGHT = 0.60        # um pillar height
FOCAL = 50.0         # um
DIAMETER = 14.0      # um aperture (kept small so the PB sweep stays quick)

air = cel.materials.air()
glass = cel.materials.fused_silica()
tio2 = cel.Material.constant(2.40 + 0.0j, "TiO2")

# === propagation-phase path =================================================
lib = cel.build_unit_cell_library(
    pillar=tio2, background=air, incident=air, substrate=glass,
    period_um=PERIOD, wavelength_um=WAVELENGTH, thickness_um=HEIGHT,
    fill_min=0.08, fill_max=0.92, n_samples=32, M=8)
print(f"library: {len(lib.fill)} pillars, "
      f"coverage = {math.degrees(lib.coverage()):.0f} deg\n")

print("=== propagation-phase path (design_metalens) ===")
profiles = {
    "focusing":  cel.PhaseProfile.focusing(FOCAL),
    "vortex l=2": cel.PhaseProfile.vortex(charge=2, focal_length_um=FOCAL),
    "deflector 10deg": cel.PhaseProfile.deflector(deflect_deg=10.0),
    "axicon 5deg": cel.PhaseProfile.axicon(axicon_deg=5.0),
}
for name, prof in profiles.items():
    d = cel.design_metalens(lib, prof, diameter_um=DIAMETER)
    print(f"  {name:<16} {d.n_cells}x{d.n_cells}  "
          f"RMS phase err {d.rms_phase_error_deg:5.1f} deg  mean|t| {d.mean_amplitude:.3f}")

# --- freeform / CGH: build an arbitrary phi(x,y) in numpy -------------------
# A linear ramp == a beam deflector; here it reproduces the analytic deflector,
# proving the freeform path. Any phi(x,y) array (a CGH) works the same way.
try:
    import numpy as np
    n = 64
    xs = np.linspace(-DIAMETER / 2, DIAMETER / 2, n)
    X, Y = np.meshgrid(xs, xs)
    k = 2 * math.pi / WAVELENGTH
    phase = k * math.sin(math.radians(10.0)) * X          # 10 deg deflector ramp
    prof = cel.PhaseProfile.freeform(phase, extent_um=DIAMETER)
    d = cel.design_metalens(lib, prof, diameter_um=DIAMETER)
    print(f"  {'freeform(ramp)':<16} {d.n_cells}x{d.n_cells}  "
          f"RMS phase err {d.rms_phase_error_deg:5.1f} deg  mean|t| {d.mean_amplitude:.3f}"
          "   <- matches the analytic deflector")
except ImportError:
    np = None
    print("  (install numpy for the freeform demo)")

# === geometric-phase / Pancharatnam-Berry path ==============================
# Find the best half-wave-plate atom, then rotate it per site (M=6: dev-fast).
print("\n=== geometric-phase path (design_pb_metalens) ===")
atom = cel.find_hwp_atom(
    pillar=tio2, background=air, incident=air, substrate=glass,
    period_um=PERIOD, wavelength_um=WAVELENGTH, thickness_um=HEIGHT,
    fill_min=0.1, fill_max=0.9, n_samples=12, M=6)
print(f"  HWP atom: fill ({atom.fill_x:.2f}, {atom.fill_y:.2f}), "
      f"retardance {atom.retardance_deg:.0f} deg, "
      f"conversion eff {atom.conversion_efficiency:.3f}")

for name, prof in [("PB focusing", cel.PhaseProfile.focusing(FOCAL)),
                   ("PB vortex l=1", cel.PhaseProfile.vortex(charge=1))]:
    pb = cel.design_pb_metalens(atom, PERIOD, WAVELENGTH, prof, diameter_um=DIAMETER)
    print(f"  {name:<14} {pb.n_cells}x{pb.n_cells}  "
          f"RMS phase err {pb.rms_phase_error_deg:.2e} deg (exact)  "
          f"conv eff {pb.conversion_efficiency:.3f}")

# --- figure: compare the propagation-path designs ---------------------------
if np is not None:
    try:
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, 4, figsize=(16, 4))
        for ax, (name, prof) in zip(axes, profiles.items()):
            d = cel.design_metalens(lib, prof, diameter_um=DIAMETER)
            im = ax.imshow(d.fill_map, cmap="viridis", origin="lower")
            ax.set(title=name, xlabel="cell i", ylabel="cell j")
            fig.colorbar(im, ax=ax, fraction=0.046)
        fig.suptitle("propagation-phase pillar maps for different phase profiles")
        fig.tight_layout()
        out = os.path.join(os.path.dirname(__file__), "05_phase_profiles.png")
        fig.savefig(out, dpi=120)
        print("\nsaved", out)
    except ImportError:
        print("\n(install matplotlib to render the figure)")
