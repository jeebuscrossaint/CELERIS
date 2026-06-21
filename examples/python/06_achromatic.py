"""Achromatic (broadband) metalens design via dispersion engineering.

A standard metalens fixes the phase only at the design wavelength, so its focal
length drifts as f(lambda) ~ f0*lambda0/lambda -- the #1 metalens limitation. To
focus a whole BAND at one plane, each site must match BOTH the base focusing
phase (mod 2pi) AND the radius-dependent GROUP DELAY (dphi/domega). A single
geometric DOF traces only a 1-D curve in the (phase, group-delay) plane, so the
library must span 2-D. Two ways:

  * fill x height grid  (build_dispersive_library)   -> multi-level / grayscale
    etch (the per-site depth is not encoded in one GDS layer);
  * SHAPE variety at ONE height (build_single_etch_library) -> fabricable in a
    single lithography step (the honest tradeoff: a smaller group-delay span).

Both feed the SAME two-objective selection (design_achromatic_metalens). Passing
gd_weight=0 reproduces a dispersion-blind "standard" design from the same library
-- the apples-to-apples baseline. verify_achromatic_focus then reports the
rigorous f(lambda) of each from the library's stored per-atom band response.

This is a small/fast configuration (D=6 um) so it runs in ~1-2 min. Scale up
diameter / fill-samples / harmonics for a real design.

Run (from the repo root, after building with -DCELERIS_BUILD_PYTHON=ON):
    PYTHONPATH=build/python python examples/python/06_achromatic.py
"""
import os
import celeris as cel

# --- design target + band ---------------------------------------------------
# Small aperture keeps the demo fast (~1-2 min). The headline achromatic figure
# of merit here is the GROUP-DELAY RMS (the lever the design directly minimizes,
# and what the engine's own selftest gates on) -- it is robust at any aperture.
# The chromatic focal-drift number is supplementary and NOISY at small apertures
# (low Fresnel number -> the on-axis focus shifts below target and the peak
# finder can lock onto a sidelobe at isolated wavelengths); it cleans up as D
# grows. Raise DIAMETER (and FOCAL to keep NA sane) to push it into that regime.
FOCAL = 20.0         # um
DIAMETER = 6.0       # um aperture (small -> fast)
LAMBDA0 = 0.532      # um center wavelength (green)
FRAC_BW = 0.20       # 20% fractional bandwidth
N_BAND = 7           # band samples
PERIOD = 0.35        # um lattice pitch
M = 6                # RCWA harmonics (dev-fast; raise for accuracy)

air = cel.materials.air()
glass = cel.materials.fused_silica()
tio2 = cel.Material.constant(2.40 + 0.0j, "TiO2")

# Band samples, ascending, centered on LAMBDA0.
lam_lo = LAMBDA0 * (1.0 - 0.5 * FRAC_BW)
lam_hi = LAMBDA0 * (1.0 + 0.5 * FRAC_BW)
band = [lam_lo + (lam_hi - lam_lo) * j / (N_BAND - 1) for j in range(N_BAND)]

print(f"achromatic target: f={FOCAL} um  D={DIAMETER} um  "
      f"band=[{lam_lo:.3f},{lam_hi:.3f}] um ({FRAC_BW*100:.0f}% BW, {N_BAND} samples)\n")


def demo(name, lib):
    """Design standard vs achromatic from one library and compare drift."""
    span = lib.gd_max_fs - lib.gd_min_fs
    print(f"=== {name} ===")
    print(f"  library: {len(lib.atoms)} atoms, group-delay span {span:.2f} fs")

    std = cel.design_achromatic_metalens(lib, FOCAL, DIAMETER, gd_weight=0.0)   # baseline
    ach = cel.design_achromatic_metalens(lib, FOCAL, DIAMETER, gd_weight=1.0)   # achromatic
    print(f"  standard   (gd_weight 0): base-phase RMS {std.rms_phase_error_deg:5.1f} deg, "
          f"GD RMS {std.rms_group_delay_error_fs:.2f} fs, mean|t| {std.mean_amplitude:.3f}")
    print(f"  achromatic (gd_weight 1): base-phase RMS {ach.rms_phase_error_deg:5.1f} deg, "
          f"GD RMS {ach.rms_group_delay_error_fs:.2f} fs, mean|t| {ach.mean_amplitude:.3f}")
    print(f"  GD budget: required {ach.required_gd_span_fs:.2f} fs, "
          f"available {ach.available_gd_span_fs:.2f} fs -> coverage {ach.gd_coverage:.2f}")
    print(f"  single height? {ach.single_height}  "
          f"(heights [{ach.min_height_um:.2f}, {ach.max_height_um:.2f}] um)")

    # The ROBUST signal that dispersion engineering worked is that the group-delay
    # objective REDUCES the group-delay RMS (the lever it directly minimizes).
    gd_factor = (std.rms_group_delay_error_fs / ach.rms_group_delay_error_fs
                 if ach.rms_group_delay_error_fs > 0 else float("inf"))
    print(f"  >> group-delay RMS: {std.rms_group_delay_error_fs:.2f} fs -> "
          f"{ach.rms_group_delay_error_fs:.2f} fs  ({gd_factor:.1f}x flatter dispersion)")

    # Rigorous f(lambda) of both -- flat = achromatic. No new RCWA solves. This
    # focal-drift number is supplementary and noisy at small apertures (see the
    # header note); GD RMS above is the metric to trust here.
    cs = cel.verify_achromatic_focus(lib, std)
    ca = cel.verify_achromatic_focus(lib, ach)
    drift_std = max(p.focal_length_um for p in cs) - min(p.focal_length_um for p in cs)
    drift_ach = max(p.focal_length_um for p in ca) - min(p.focal_length_um for p in ca)
    print(f"  {'lambda(nm)':>10}  {'standard f(um)':>16}  {'achromatic f(um)':>16}")
    for s, a in zip(cs, ca):
        print(f"  {s.wavelength_um*1000:>10.0f}  {s.focal_length_um:>16.2f}  "
              f"{a.focal_length_um:>16.2f}")
    print(f"  chromatic focal drift (supplementary, noisy at small D): "
          f"standard {drift_std:.2f} um -> achromatic {drift_ach:.2f} um\n")
    return std, ach, cs, ca


# === 1. fill x height library (multi-etch) ==================================
hx_lib = cel.build_dispersive_library(
    pillar=tio2, background=air, incident=air, substrate=glass,
    period_um=PERIOD, band_wavelengths_um=band, center_wavelength_um=LAMBDA0,
    fill_min=0.08, fill_max=0.92, n_fills=10, thick_lo=0.40, thick_hi=1.40,
    n_heights=6, M=M)
hx_std, hx_ach, hx_cs, hx_ca = demo("fill x height (multi-etch / grayscale)", hx_lib)

# === 2. single-etch library (shape diversity at one height) =================
# Taller pillar widens the group-delay span each shape can supply (h~1.1um at
# Lambda=0.35 ~ AR 4, fabricable) -- the key single-etch knob.
se_lib = cel.build_single_etch_library(
    pillar=tio2, background=air, incident=air, substrate=glass,
    period_um=PERIOD, band_wavelengths_um=band, center_wavelength_um=LAMBDA0,
    fill_min=0.08, fill_max=0.92, n_fills=10, thickness_um=1.10, M=M)
se_std, se_ach, se_cs, se_ca = demo("single-etch (shape-diverse, one 1.10um etch)", se_lib)

# Show which shapes the single-etch achromatic design leaned on -- shape
# diversity (not depth) is what supplied the group delay.
shape_names = {cel.MetaShape.Rectangle: "square", cel.MetaShape.Ellipse: "circle",
               cel.MetaShape.Cross: "cross", cel.MetaShape.Ring: "ring"}
mix = {n: 0 for n in shape_names.values()}
for q in se_ach.atom_index:
    mix[shape_names[se_lib.atoms[q].shape]] += 1
print("single-etch shape mix (achromatic): "
      + ", ".join(f"{v} {k}" for k, v in mix.items()))

# --- figure: chromatic focal curves (reuse the verifies computed above) -----
try:
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 5))
    for (cs, ca), label, ls in (((hx_cs, hx_ca), "fill x height", "-"),
                                ((se_cs, se_ca), "single-etch", "--")):
        wl = [p.wavelength_um * 1000 for p in cs]
        ax.plot(wl, [p.focal_length_um for p in cs], ls, color="C3",
                label=f"standard ({label})")
        ax.plot(wl, [p.focal_length_um for p in ca], ls, color="C0",
                label=f"achromatic ({label})")
    ax.axhline(FOCAL, color="k", lw=0.8, ls=":", label="target f")
    ax.set(xlabel="wavelength (nm)", ylabel="focal length (um)",
           title="Chromatic focal shift: standard vs achromatic")
    ax.legend(fontsize=8)
    fig.tight_layout()
    out = os.path.join(os.path.dirname(__file__), "06_achromatic.png")
    fig.savefig(out, dpi=120)
    print("\nsaved", out)
except ImportError:
    print("\n(install matplotlib to render the focal-shift figure)")
