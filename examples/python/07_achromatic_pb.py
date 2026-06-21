"""Achromatic Pancharatnam-Berry metalens -- geometric phase + dispersion.

This is the MODERN achromatic recipe. A standard metalens fixes the phase only
at the design wavelength, so its focal length drifts with lambda. To focus a
whole band at one plane each site must match BOTH the base focusing phase AND the
radius-dependent GROUP DELAY (dphi/domega).

The propagation-phase achromat (see 06_achromatic.py) makes ONE atom hit both
terms -- a two-objective fit that leaves a base-phase residual because the library
is discrete. The PB achromat DECOUPLES them:

  * the GEOMETRIC (PB) phase -- rotating a birefringent atom by theta imprints a
    phase -handedness*2*theta on the spin-flipped output that is EXACT and
    wavelength-INDEPENDENT, so the rotation hits the base phase for WHATEVER atom
    is placed (base-phase RMS ~0 by construction);
  * the atom is then chosen PURELY for its group delay (ONE objective) from a
    dispersive birefringent library at ONE etch depth (rotated rectangles of
    varying footprint -> a single fabrication step, no grayscale).

So the achromatic limit here is ONLY the library's group-delay coverage. Passing
gd_weight=0 reproduces a dispersion-blind "standard" PB lens from the same library
-- the apples-to-apples baseline. verify_pb_achromatic_focus then reports the
rigorous f(lambda) of each from the library's stored per-atom band response.

This is a small/fast configuration (D=6 um) so it runs in ~1-2 min. Scale up
diameter / fill-samples / harmonics for a real design.

Run (from the repo root, after building with -DCELERIS_BUILD_PYTHON=ON):
    PYTHONPATH=build/python python examples/python/07_achromatic_pb.py
"""
import os
import celeris as cel

# --- design target + band ---------------------------------------------------
# As in 06_achromatic.py, the headline figure of merit is the GROUP-DELAY RMS
# (the lever the design directly minimizes, and what the engine's selftest gates
# on) -- robust at any aperture. The chromatic focal-drift number is supplementary
# and noisy at small apertures (low Fresnel number); it cleans up as D grows.
FOCAL = 30.0         # um
DIAMETER = 10.0      # um aperture (small -> fast)
LAMBDA0 = 0.532      # um center wavelength (green)
FRAC_BW = 0.20       # 20% fractional bandwidth
N_BAND = 7           # band samples
PERIOD = 0.35        # um lattice pitch
# Birefringent atoms accumulate group delay with optical path, so a taller pillar
# widens the GD span each footprint can supply (h~1.1um at Lambda=0.35 ~ AR 4,
# fabricable). The etch DEPTH is the dominant single-etch knob.
HEIGHT = 1.10        # um single etch depth
M = 6                # RCWA harmonics (dev-fast; raise for accuracy)
N_FILLS = 12         # fill_x x fill_y grid side

air = cel.materials.air()
glass = cel.materials.bk7()
tio2 = cel.Material.constant(2.40 + 0.0j, "TiO2")

# Band samples, ascending, centered on LAMBDA0.
lam_lo = LAMBDA0 * (1.0 - 0.5 * FRAC_BW)
lam_hi = LAMBDA0 * (1.0 + 0.5 * FRAC_BW)
band = [lam_lo + (lam_hi - lam_lo) * j / (N_BAND - 1) for j in range(N_BAND)]

print(f"achromatic PB target: f={FOCAL} um  D={DIAMETER} um  "
      f"band=[{lam_lo:.3f},{lam_hi:.3f}] um ({FRAC_BW*100:.0f}% BW, {N_BAND} samples)")
print(f"  SINGLE etch depth {HEIGHT} um, RCP illumination, {N_FILLS}x{N_FILLS} fill grid\n")

# --- dispersive BIREFRINGENT library (one etch depth) -----------------------
# Two RCWA solves/atom/wavelength -> spin-flip a_cross + its group delay. The
# builder FILTERS out near-square (non-birefringent) atoms whose a_cross ~ 0 give
# a garbage group delay that would pollute the selection.
print("building dispersive birefringent library (two solves/atom/wavelength)...")
lib = cel.build_dispersive_pb_library(
    pillar=tio2, background=air, incident=air, substrate=glass,
    period_um=PERIOD, band_wavelengths_um=band, center_wavelength_um=LAMBDA0,
    fill_min=0.10, fill_max=0.90, n_fills=N_FILLS, thickness_um=HEIGHT, M=M)
span = lib.gd_max_fs - lib.gd_min_fs
print(f"  kept {len(lib.atoms)} birefringent atoms (of {lib.n_fill**2} grid points), "
      f"group-delay span {span:.2f} fs\n")

# --- standard (chromatic) vs achromatic, from the SAME library --------------
# gd_weight=0 = best-conversion atom everywhere = a STANDARD PB lens (baseline);
# gd_weight=1 = vary the atom per radius to match the group delay.
std = cel.design_pb_achromatic_metalens(lib, FOCAL, DIAMETER, handedness=1, gd_weight=0.0)
ach = cel.design_pb_achromatic_metalens(lib, FOCAL, DIAMETER, handedness=1, gd_weight=1.0)

print(f"  standard   (gd_weight 0): base-phase RMS {std.rms_phase_error_deg:.2e} deg, "
      f"GD RMS {std.rms_group_delay_error_fs:.2f} fs, mean|a_cross| {std.mean_amplitude:.3f}, "
      f"conv {std.mean_conversion:.3f}")
print(f"  achromatic (gd_weight 1): base-phase RMS {ach.rms_phase_error_deg:.2e} deg, "
      f"GD RMS {ach.rms_group_delay_error_fs:.2f} fs, mean|a_cross| {ach.mean_amplitude:.3f}, "
      f"conv {ach.mean_conversion:.3f}")
print("  base-phase RMS ~0 for BOTH: the rotation stamps the base phase EXACTLY "
      "(geometric phase) -- the atom is free to chase group delay.")
print(f"  GD budget: required {ach.required_gd_span_fs:.2f} fs, "
      f"available {ach.available_gd_span_fs:.2f} fs -> coverage {ach.gd_coverage:.2f}")
if ach.gd_coverage < 1.0:
    print("  HONEST: GD coverage < 1 -> achromatic only over a reduced aperture/bandwidth "
          "(taller/coupled atoms widen the span).")

# The ROBUST signal that dispersion engineering worked is that the group-delay
# objective REDUCES the group-delay RMS (the lever it directly minimizes).
gd_factor = (std.rms_group_delay_error_fs / ach.rms_group_delay_error_fs
             if ach.rms_group_delay_error_fs > 0 else float("inf"))
print(f"  >> group-delay RMS: {std.rms_group_delay_error_fs:.2f} fs -> "
      f"{ach.rms_group_delay_error_fs:.2f} fs  ({gd_factor:.1f}x flatter dispersion)\n")

# --- rigorous f(lambda) of both (no new RCWA solves) ------------------------
cs = cel.verify_pb_achromatic_focus(lib, std)
ca = cel.verify_pb_achromatic_focus(lib, ach)
drift_std = max(p.focal_length_um for p in cs) - min(p.focal_length_um for p in cs)
drift_ach = max(p.focal_length_um for p in ca) - min(p.focal_length_um for p in ca)
print(f"  {'lambda(nm)':>10}  {'standard f(um)':>16}  {'achromatic f(um)':>16}")
for s, a in zip(cs, ca):
    print(f"  {s.wavelength_um*1000:>10.0f}  {s.focal_length_um:>16.2f}  "
          f"{a.focal_length_um:>16.2f}")
print(f"  chromatic focal drift (supplementary, noisy at small D): "
      f"standard {drift_std:.2f} um -> achromatic {drift_ach:.2f} um")

# --- figure: chromatic focal curves -----------------------------------------
try:
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 5))
    wl = [p.wavelength_um * 1000 for p in cs]
    ax.plot(wl, [p.focal_length_um for p in cs], "-o", color="C3", label="standard PB")
    ax.plot(wl, [p.focal_length_um for p in ca], "-o", color="C0", label="achromatic PB")
    ax.axhline(FOCAL, color="k", lw=0.8, ls=":", label="target f")
    ax.set(xlabel="wavelength (nm)", ylabel="focal length (um)",
           title="Achromatic Pancharatnam-Berry: chromatic focal shift")
    ax.legend(fontsize=9)
    fig.tight_layout()
    out = os.path.join(os.path.dirname(__file__), "07_achromatic_pb.png")
    fig.savefig(out, dpi=120)
    print("\nsaved", out)
except ImportError:
    print("\n(install matplotlib to render the focal-shift figure)")
