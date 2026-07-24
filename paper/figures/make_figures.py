#!/usr/bin/env python3
"""Generate the figures for the CELERIS paper from live CELERIS runs.

Outputs (vector PDF, for LaTeX) into this directory:
  fig_validation.pdf  -- CELERIS vs Stanford S4 convergence + agreement (external validation)
  fig_library.pdf     -- meta-atom phase library (transmission phase & |t| vs fill)
  fig_psf.pdf         -- diffraction-limited focal spot (PSF map + line cut)

Requires: celeris (pip install . from repo root), numpy, matplotlib. The S4 panel is
optional -- if S4 is not importable, that panel shows CELERIS-only convergence.

Run:  python paper/figures/make_figures.py
"""
import math
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import celeris as cel

HERE = os.path.dirname(os.path.abspath(__file__))
# Colorblind-safe (Wong) palette.
BLUE, VERM, GREEN, GRAY = "#0072B2", "#D55E00", "#009E73", "#555555"
plt.rcParams.update({"font.size": 11, "axes.grid": True, "grid.alpha": 0.3,
                     "figure.dpi": 130, "savefig.bbox": "tight"})


def fig_validation():
    """CELERIS vs S4 on the identical freestanding n=1.5 grating (selftest [5] case)."""
    Lam, fill, thick, wl = 0.3, 0.5, 0.5, 0.5
    glass = cel.Material.constant(complex(1.5, 0.0)); air = cel.materials.air()
    Ms = [4, 6, 8, 10, 14, 20, 30, 40]

    def cel_t0(M, pol):
        g = cel.BinaryGrating1D(glass, air, Lam, fill, thick)
        r = cel.solve_rcwa_1d(air, g, air, wl, 0.0, M, pol)
        return r.de_t[len(r.orders) // 2]

    cel_te = [cel_t0(M, cel.Pol.TE) for M in Ms]
    cel_tm = [cel_t0(M, cel.Pol.TM) for M in Ms]

    s4_te = s4_tm = None
    try:
        import S4
        def s4_t0(M, s, p):
            S = S4.New(Lattice=Lam, NumBasis=M)
            S.SetMaterial(Name='V', Epsilon=1.0); S.SetMaterial(Name='G', Epsilon=2.25)
            S.AddLayer(Name='t', Thickness=0, Material='V')
            S.AddLayer(Name='g', Thickness=thick, Material='V')
            S.SetRegionRectangle(Layer='g', Material='G', Center=(0, 0), Angle=0,
                                 Halfwidths=(fill * Lam / 2, 0))
            S.AddLayer(Name='b', Thickness=0, Material='V')
            S.SetFrequency(1.0 / wl)
            S.SetExcitationPlanewave(IncidenceAngles=(0, 0), sAmplitude=s, pAmplitude=p)
            ft, _ = S.GetPowerFlux(Layer='t', zOffset=0)
            fb, _ = S.GetPowerFlux(Layer='b', zOffset=0)
            return (fb / ft).real
        s4_te = [s4_t0(M, 1, 0) for M in Ms]
        s4_tm = [s4_t0(M, 0, 1) for M in Ms]
    except Exception as e:  # noqa: BLE001
        print("  [fig_validation] S4 not available, CELERIS-only panel:", e)

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(10, 4))
    axL.plot(Ms, cel_te, "o-", color=BLUE, label="CELERIS (TE)")
    axL.plot(Ms, cel_tm, "s-", color=VERM, label="CELERIS (TM)")
    if s4_te:
        axL.plot(Ms, s4_te, "x--", color=BLUE, mew=2, label="S$^4$ (TE)")
        axL.plot(Ms, s4_tm, "+--", color=VERM, mew=2, label="S$^4$ (TM)")
    axL.axhline(0.93333, color=GRAY, ls=":", lw=1, label="grcwa (TE) 0.93333")
    axL.set(xlabel="Fourier order $M$", ylabel="0th-order transmittance $T_0$",
            title="Convergence vs external solvers")
    axL.legend(fontsize=8)

    if s4_te:
        axR.semilogy(Ms, np.abs(np.array(cel_te) - np.array(s4_te)), "o-", color=BLUE,
                     label="|CELERIS $-$ S$^4$| (TE)")
        axR.semilogy(Ms, np.abs(np.array(cel_tm) - np.array(s4_tm)), "s-", color=VERM,
                     label="|CELERIS $-$ S$^4$| (TM)")
        axR.set(xlabel="Fourier order $M$", ylabel="|$T_0$ difference|",
                title="CELERIS $-$ S$^4$ agreement")
        axR.legend(fontsize=9)
    else:
        axR.axis("off")
    fig.tight_layout()
    out = os.path.join(HERE, "fig_validation.pdf"); fig.savefig(out); plt.close(fig)
    print("  wrote", out)


def _tio2_library():
    air = cel.materials.air(); glass = cel.materials.fused_silica()
    tio2 = cel.Material.constant(2.40 + 0.0j, "TiO2")
    return cel.build_unit_cell_library(
        pillar=tio2, background=air, incident=air, substrate=glass,
        period_um=0.35, wavelength_um=0.532, thickness_um=0.60,
        fill_min=0.08, fill_max=0.92, n_samples=40, M=8), 0.532


def fig_library():
    lib, wl = _tio2_library()
    fill = np.array(lib.fill)
    phase = np.degrees(np.array(lib.phase)) % 360.0
    amp = np.array(lib.amplitude)

    fig, ax1 = plt.subplots(figsize=(6.2, 4.2))
    ax1.plot(fill, phase, "o-", color=BLUE, label="phase")
    ax1.set(xlabel="fill fraction", ylabel="transmission phase (deg)",
            title=f"TiO$_2$ meta-atom library ($\\lambda$={wl*1e3:.0f} nm, $\\Lambda$=350 nm)")
    ax1.set_ylim(0, 360)
    ax2 = ax1.twinx()
    ax2.plot(fill, amp, "s--", color=VERM, label="|t|")
    ax2.set_ylabel("transmission |t|", color=VERM); ax2.set_ylim(0, 1.05)
    ax2.tick_params(axis="y", labelcolor=VERM); ax2.grid(False)
    ax1.text(0.05, 0.9, f"coverage = {math.degrees(lib.coverage()):.0f}$^\\circ$",
             transform=ax1.transAxes, fontsize=10)
    fig.tight_layout()
    out = os.path.join(HERE, "fig_library.pdf"); fig.savefig(out); plt.close(fig)
    print("  wrote", out)


def fig_psf():
    lib, wl = _tio2_library()
    focal, diameter = 50.0, 20.0
    lens = cel.design_metalens(lib, focal_length_um=focal, diameter_um=diameter)
    foc = cel.analyze_focus(lens, lib, focal, wl, diameter)
    psf = cel.compute_psf(lens, lib, focal, wl, diameter, n=201, half_window_um=4.0)
    I = np.array(psf.intensity); I /= I.max()
    W = psf.half_window_um
    x = np.linspace(-W, W, I.shape[0])

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(10, 4.2))
    im = axA.imshow(np.log10(I + 1e-6), extent=[-W, W, -W, W], cmap="inferno", origin="lower")
    axA.set(title=f"focal-plane PSF (log)  Strehl {foc.strehl:.2f}",
            xlabel="x ($\\mu$m)", ylabel="y ($\\mu$m)")
    fig.colorbar(im, ax=axA, label="$\\log_{10}|E|^2$ (norm)")

    axB.plot(x, I[I.shape[0] // 2, :], color=BLUE, label="CELERIS focus")
    axB.axhline(0.5, color=GRAY, ls=":", lw=1)
    axB.set(title=f"central line cut\nFWHM {foc.fwhm_um:.3f} $\\mu$m "
                  f"(diffraction limit {foc.diffraction_limit_um:.3f} $\\mu$m)",
            xlabel="x ($\\mu$m)", ylabel="intensity (norm)")
    axB.legend(fontsize=9)
    fig.tight_layout()
    out = os.path.join(HERE, "fig_psf.pdf"); fig.savefig(out); plt.close(fig)
    print("  wrote", out, f"[Strehl {foc.strehl:.3f}, FWHM {foc.fwhm_um:.3f} vs DL {foc.diffraction_limit_um:.3f}]")


def fig_achromatic():
    """Chen-2018-style broadband achromatic PB metalens: chromatic focal shift,
    standard vs group-delay-engineered, from one single-etch birefringent library."""
    air = cel.materials.air(); glass = cel.materials.bk7()
    tio2 = cel.Material.constant(2.40 + 0.0j, "TiO2")
    focal, diameter, lam0, bw, nband, period, height, M, nfill = \
        30.0, 10.0, 0.532, 0.20, 7, 0.35, 1.10, 6, 12
    lo, hi = lam0 * (1 - 0.5 * bw), lam0 * (1 + 0.5 * bw)
    band = [lo + (hi - lo) * j / (nband - 1) for j in range(nband)]
    lib = cel.build_dispersive_pb_library(
        pillar=tio2, background=air, incident=air, substrate=glass,
        period_um=period, band_wavelengths_um=band, center_wavelength_um=lam0,
        fill_min=0.10, fill_max=0.90, n_fills=nfill, thickness_um=height, M=M)
    std = cel.design_pb_achromatic_metalens(lib, focal, diameter, handedness=1, gd_weight=0.0)
    ach = cel.design_pb_achromatic_metalens(lib, focal, diameter, handedness=1, gd_weight=1.0)
    cs = cel.verify_pb_achromatic_focus(lib, std)
    ca = cel.verify_pb_achromatic_focus(lib, ach)
    wl = [p.wavelength_um * 1000 for p in cs]
    drift_s = max(p.focal_length_um for p in cs) - min(p.focal_length_um for p in cs)
    drift_a = max(p.focal_length_um for p in ca) - min(p.focal_length_um for p in ca)

    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    ax.plot(wl, [p.focal_length_um for p in cs], "-o", color=VERM,
            label=f"standard PB (drift {drift_s:.1f} $\\mu$m)")
    ax.plot(wl, [p.focal_length_um for p in ca], "-o", color=BLUE,
            label=f"achromatic PB (drift {drift_a:.1f} $\\mu$m)")
    ax.axhline(focal, color=GRAY, ls=":", lw=1, label="target $f$")
    ax.set(xlabel="wavelength (nm)", ylabel="focal length ($\\mu$m)",
           title="Broadband achromatic metalens (Chen 2018 recipe):\n"
                 "single-etch group-delay engineering flattens the focal shift")
    ax.legend(fontsize=9)
    fig.tight_layout()
    out = os.path.join(HERE, "fig_achromatic.pdf"); fig.savefig(out); plt.close(fig)
    print("  wrote", out, f"[focal drift {drift_s:.2f} -> {drift_a:.2f} um]")


if __name__ == "__main__":
    print("Generating CELERIS paper figures...")
    fig_validation()
    fig_library()
    fig_psf()
    fig_achromatic()
    print("done.")
