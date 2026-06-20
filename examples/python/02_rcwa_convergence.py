"""RCWA convergence: zeroth-order transmittance and energy vs harmonic count M.

A standard sanity check a CREOL tester will run: does the solver converge, and
is energy conserved? Sweeps M for a TiO2 square pillar (the 2D vectorial solver,
post-Li factorization) and prints/plots T0 and R+T vs M. Expect T0 to settle and
sum_de = 1.000000 at every M.

    PYTHONPATH=build/python python examples/python/02_rcwa_convergence.py
"""
import os
import celeris as cel

WAVELENGTH, PERIOD, HEIGHT, FILL = 0.532, 0.35, 0.60, 0.5
air, glass = cel.materials.air(), cel.materials.fused_silica()
tio2 = cel.Material.constant(2.40 + 0.0j, "TiO2")

cell = cel.RectCell2D(tio2, air, FILL, FILL, HEIGHT)        # square pillar
stack = cel.Rcwa2DStack(PERIOD, PERIOD, [cell])

Ms = list(range(2, 15))
t0, energy = [], []
print(f"{'M':>3}  {'T0':>9}  {'R+T':>10}")
for M in Ms:
    res = cel.solve_rcwa_2d(air, stack, glass, WAVELENGTH, 0.0, 0.0,
                            1.0 + 0j, 0.0 + 0j, M, M)
    t0.append(res.de_t0)
    energy.append(res.sum_de)
    print(f"{M:>3}  {res.de_t0:9.5f}  {res.sum_de:10.7f}")

spread = max(t0[len(t0) // 2:]) - min(t0[len(t0) // 2:])
print(f"\nT0 spread over the upper half of the M-sweep: {spread:.4f} (small = converged)")

try:
    import matplotlib.pyplot as plt
    fig, ax1 = plt.subplots(figsize=(7, 4.5))
    ax1.plot(Ms, t0, "o-", color="C0", label="zeroth-order T")
    ax1.set(xlabel="harmonic half-count M", ylabel="zeroth-order transmittance")
    ax2 = ax1.twinx()
    ax2.plot(Ms, energy, "s--", color="C3", label="R + T")
    ax2.set_ylabel("energy R + T")
    ax2.set_ylim(0.999, 1.001)
    ax1.set_title("2D RCWA convergence — TiO2 square pillar")
    lines = ax1.get_lines() + ax2.get_lines()
    ax1.legend(lines, [l.get_label() for l in lines], loc="center right")
    fig.tight_layout()
    out = os.path.join(os.path.dirname(__file__), "02_rcwa_convergence.png")
    fig.savefig(out, dpi=130)
    print("saved", out)
except ImportError:
    print("(install matplotlib to render the figure)")
