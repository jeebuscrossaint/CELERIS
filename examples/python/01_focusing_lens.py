"""Design a TiO2-on-glass focusing metalens and analyze its focus.

The headline workflow: build a phase library by sweeping pillar fill, design a
hyperbolic focusing lens, then propagate to the focal plane and report the
standard optical metrics. Saves a PSF + pillar-map figure next to this script.

Run (from the repo root, after building with -DCELERIS_BUILD_PYTHON=ON):
    PYTHONPATH=build/python python examples/python/01_focusing_lens.py
"""
import math
import os
import celeris as cel

# --- design parameters (a green-light, NA~0.2 lens) -------------------------
WAVELENGTH = 0.532   # um
PERIOD = 0.35        # um lattice pitch (subwavelength)
HEIGHT = 0.60        # um pillar height
FOCAL = 50.0         # um
DIAMETER = 20.0      # um aperture  -> NA ~ 0.20

# --- materials --------------------------------------------------------------
air = cel.materials.air()
glass = cel.materials.fused_silica()
tio2 = cel.Material.constant(2.40 + 0.0j, "TiO2")   # ~visible TiO2 index

# --- phase library: sweep square-pillar fill --------------------------------
lib = cel.build_unit_cell_library(
    pillar=tio2, background=air, incident=air, substrate=glass,
    period_um=PERIOD, wavelength_um=WAVELENGTH, thickness_um=HEIGHT,
    fill_min=0.08, fill_max=0.92, n_samples=32, M=8)
print(f"library: {len(lib.fill)} pillars, "
      f"phase coverage = {math.degrees(lib.coverage()):.0f} deg, "
      f"mean |t| = {sum(lib.amplitude) / len(lib.amplitude):.3f}")

# --- design + analyze -------------------------------------------------------
lens = cel.design_metalens(lib, focal_length_um=FOCAL, diameter_um=DIAMETER)
print(f"design: {lens.n_cells}x{lens.n_cells} pillars, "
      f"RMS phase error = {lens.rms_phase_error_deg:.1f} deg, "
      f"mean |t| = {lens.mean_amplitude:.3f}")

foc = cel.analyze_focus(lens, lib, FOCAL, WAVELENGTH, DIAMETER)
print(f"focus: Strehl = {foc.strehl:.3f}  (phase-only {foc.phase_strehl:.3f}),  "
      f"FWHM = {foc.fwhm_um:.3f} um  (diffraction limit {foc.diffraction_limit_um:.3f}),  "
      f"encircled energy = {foc.encircled_energy:.3f}")

# --- PSF map ----------------------------------------------------------------
psf = cel.compute_psf(lens, lib, FOCAL, WAVELENGTH, DIAMETER, n=128, half_window_um=4.0)

try:
    import numpy as np
    import matplotlib.pyplot as plt

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))
    ext = [-psf.half_window_um, psf.half_window_um] * 2
    im = ax1.imshow(np.log10(psf.intensity / psf.intensity.max() + 1e-6),
                    extent=ext, cmap="inferno", origin="lower")
    ax1.set(title=f"focal-plane PSF (log)\nFWHM {foc.fwhm_um:.2f} um, Strehl {foc.strehl:.2f}",
            xlabel="x (um)", ylabel="y (um)")
    fig.colorbar(im, ax=ax1, label="log10 |E|^2 (norm)")

    fm = ax2.imshow(lens.fill_map, cmap="viridis", origin="lower")
    ax2.set(title=f"pillar fill map ({lens.n_cells}x{lens.n_cells})",
            xlabel="cell i", ylabel="cell j")
    fig.colorbar(fm, ax=ax2, label="fill fraction")
    fig.tight_layout()
    out = os.path.join(os.path.dirname(__file__), "01_focusing_lens.png")
    fig.savefig(out, dpi=130)
    print("saved", out)
except ImportError:
    print("(install matplotlib + numpy to render the figure)")
