"""Use real tabulated material data (TiO2, Siefke 2016) in a design.

Loads the bundled n,k CSV with numpy, builds a dispersive Material, and shows
how the achievable phase coverage shifts with wavelength across the visible band
when you account for real dispersion + loss (k). This is the honest material
path: no hardcoded index, real refractiveindex.info data.

    PYTHONPATH=build/python python examples/python/03_material_dispersion.py
"""
import os
import celeris as cel

REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
CSV = os.path.join(REPO, "data", "TiO2_Siefke.csv")

import numpy as np
wl, n, k = np.loadtxt(CSV, comments="#", unpack=True)   # columns: wavelength_um n k
print(f"loaded {len(wl)} points from {os.path.basename(CSV)}: "
      f"{wl.min():.3f}-{wl.max():.3f} um, n@0.53um ~ {np.interp(0.53, wl, n):.3f}")

tio2 = cel.Material.tabulated(list(wl), list(n), list(k), "TiO2 (Siefke 2016)")
air, glass = cel.materials.air(), cel.materials.fused_silica()

# How does phase coverage / transmittance vary with wavelength (fixed geometry)?
PERIOD, HEIGHT = 0.35, 0.60
import math
print(f"\n{'lambda(um)':>10}  {'n+ik':>16}  {'coverage(deg)':>13}  {'mean|t|':>8}")
for lam in (0.45, 0.50, 0.532, 0.60, 0.65):
    idx = tio2.index(lam)
    lib = cel.build_unit_cell_library(tio2, air, air, glass, PERIOD, lam, HEIGHT,
                                      0.1, 0.9, 20, 8)
    mean_t = sum(lib.amplitude) / len(lib.amplitude)
    print(f"{lam:>10.3f}  {idx.real:6.3f}+{idx.imag:6.4f}i  "
          f"{math.degrees(lib.coverage()):>13.0f}  {mean_t:>8.3f}")

try:
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(wl * 1000, n, "C0-", label="n")
    ax.set(xlabel="wavelength (nm)", ylabel="n", title="TiO2 dispersion (Siefke 2016, ALD amorphous)")
    ax2 = ax.twinx()
    ax2.plot(wl * 1000, k, "C3-", label="k")
    ax2.set_ylabel("k (extinction)")
    ax.axvspan(380, 700, color="0.9", zorder=0)
    fig.tight_layout()
    out = os.path.join(os.path.dirname(__file__), "03_material_dispersion.png")
    fig.savefig(out, dpi=130)
    print("\nsaved", out)
except ImportError:
    print("(install matplotlib to render the figure)")
