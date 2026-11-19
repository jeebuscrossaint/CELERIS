#!/usr/bin/env python3
"""Strehl-based achievable aperture D_S for an R1 library (issue #33). PROTOTYPE.

The span-based D_lib (metric.md) assumes every delay in the span is usable at any
phase. Propagation-phase libraries couple phase and group delay, so this script
measures what a design actually reaches. Per band wavelength it computes the on-axis
Strehl at the design focus from the chosen atoms' stored band response:

    S(lambda) = |sum_s t_s(lambda) exp(i k sag_s)|^2 / (sum_s |t_s(lambda)|)^2

The constant phase and group-delay offsets drop out of S. Atom selection uses the
CELERIS picker (a numpy copy, checked against design_achromatic_metalens) and scans
both offsets. D_S is the largest aperture whose band-minimum S stays >= 0.8.

Status: prop libraries only. PB needs DispersivePbAtom.cross bound to Python.

Run from the repo root:
    python research/r1/coverage.py tio2 sio2 1.0
"""
import argparse
import math
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import metric  # noqa: E402
import sweep  # noqa: E402

C = metric.C_UM_PER_FS
DIAMETERS_UM = (4, 8, 12, 16, 20, 26.4, 32, 40)
STREHL_CUT = 0.8  # Marechal criterion


class PropLibrary:
    """Arrays of one propagation-phase library: phases, delays, band response."""

    def __init__(self, lib, period_um: float):
        a = lib.atoms
        self.period = period_um
        self.ph0 = np.array([x.phase0_rad for x in a])
        self.gd = np.array([x.group_delay_fs for x in a])
        self.mamp = np.array([x.mean_amplitude for x in a])
        self.T = np.array([[x.amplitude[j] * np.exp(1j * x.phase_rad[j])
                            for j in range(len(lib.wavelengths_um))] for x in a])
        self.lam = np.array(lib.wavelengths_um)
        w = 2 * np.pi * C / self.lam
        self.dw_half = 0.5 * abs(w[0] - w[-1])
        self.k0 = 2 * np.pi / lib.center_wavelength_um
        self.gd_min, self.gd_max = lib.gd_min_fs, lib.gd_max_fs

    def sites(self, D: float):
        """Radii of the square-grid sites (as in the C++ design) and the in-aperture mask."""
        n = max(1, round(D / self.period))
        c = (n - 1) / 2
        x = (np.arange(n) - c) * self.period
        X, Y = np.meshgrid(x, x)
        r = np.hypot(X, Y).ravel()
        return r, r <= D / 2

    def pick(self, sag, gd_ref, phi_off, gd_weight=1.0, amp_weight=0.25):
        """Per-site atom index; same cost as design_achromatic_metalens."""
        tphi = -self.k0 * sag + phi_off
        tgd = -sag / C + gd_ref
        dphi = np.angle(np.exp(1j * (self.ph0[None, :] - tphi[:, None])))
        dgd = (self.gd[None, :] - tgd[:, None]) * self.dw_half
        cost = dphi**2 + gd_weight * dgd**2 + amp_weight * (1 - self.mamp[None, :])**2
        return np.argmin(cost, axis=1)

    def strehl(self, idx, sag, inside):
        """On-axis Strehl at the design focus, one value per band wavelength."""
        t = self.T[idx][inside]
        k = 2 * np.pi / self.lam
        field = (t * np.exp(1j * np.outer(sag[inside], k))).sum(0)
        return np.abs(field)**2 / np.abs(t).sum(0)**2

    def best_strehl(self, D, focal, gd_weight, n_gd=13, n_phi=8):
        """Band Strehl of the best pick over the offset scan (by band-minimum S)."""
        r, inside = self.sites(D)
        sag = np.sqrt(r * r + focal * focal) - focal
        req = (math.hypot(D / 2, focal) - focal) / C
        best = None
        for gref in np.linspace(self.gd_min, self.gd_max + req, n_gd):
            for po in np.linspace(0, 2 * np.pi, n_phi, endpoint=False):
                s = self.strehl(self.pick(sag, gref, po, gd_weight), sag, inside)
                if best is None or s.min() > best.min():
                    best = s
        return req, best


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("material")
    ap.add_argument("substrate")
    ap.add_argument("height_um", type=float)
    ap.add_argument("--config", default="full", choices=sorted(sweep.CONFIGS))
    args = ap.parse_args()

    import celeris as cel
    cfg = sweep.CONFIGS[args.config]
    lam_min, lam_max = cfg["band"]
    _, wc, _ = metric.band(lam_min, lam_max)
    lam_c = 2 * math.pi * C / wc
    pillar = cel.materials.by_name(args.material)
    substrate = cel.materials.by_name(args.substrate)
    P = sweep.choose_period_um(substrate, lam_min, cfg["period_margin"])
    fmin, fmax = sweep.fill_range(P, cfg["min_feature_um"], cfg["min_gap_um"])
    band = sweep.wavelengths(lam_min, lam_max, cfg["n_band"])
    lib, _ = sweep.build(cel, "prop", pillar, substrate, P, band, lam_c, fmin, fmax,
                         cfg["n_fill"], args.height_um, cfg["M"])
    L = PropLibrary(lib, P)
    F = sweep.FOCAL_UM

    # The numpy picker must reproduce the C++ one at the C++ offsets.
    for D in (10.0, 25.0):
        r, _ = L.sites(D)
        sag = np.sqrt(r * r + F * F) - F
        req = (math.hypot(D / 2, F) - F) / C
        gref = 0.5 * (L.gd_min + L.gd_max) + 0.5 * req
        d = cel.design_achromatic_metalens(lib, F, D)
        if not np.array_equal(L.pick(sag, gref, 0.0), np.array(d.atom_index)):
            print(f"picker mismatch vs design_achromatic_metalens at D={D}", file=sys.stderr)
            return 1

    span = L.gd.max() - L.gd.min()
    print(f"{args.material}/{args.substrate} H={args.height_um} P={P} atoms={len(lib.atoms)} "
          f"span={span:.2f} fs  D_lib={metric.achievable_diameter_um(span, F):.1f} um")
    for gdw, label in ((0.0, "standard"), (1.0, "achromat")):
        d_s = None
        for D in DIAMETERS_UM:
            req, s = L.best_strehl(D, F, gdw)
            if s.min() >= STREHL_CUT:
                d_s = D
            print(f"  {label:8s} D {D:5.1f} um  required {req:5.2f} fs  S_min {s.min():.3f}")
        print(f"  {label}: D_S(0.8) >= {d_s} um (grid {DIAMETERS_UM})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
