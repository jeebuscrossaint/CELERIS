#!/usr/bin/env python3
"""Strehl-based achievable aperture D_S for R1 libraries (issue #33).

The span-based D_lib (metric.md) assumes every delay in the span is usable at any
phase. This script measures what a design actually reaches. Per band wavelength it
computes the on-axis Strehl at the design focus from the chosen atoms' stored band
response:

    S(lambda) = |sum_s u_s(lambda) exp(i k sag_s)|^2 / (sum_s |u_s(lambda)|)^2

where u_s is the field a site contributes:

- prop: u_s = t_s(lambda), the atom's transmission.
- pb:   u_s = a_cross,s(lambda) exp(i phi_geo,s), with the rotation phase
        phi_geo,s = -k0 sag_s - arg a_cross,s(lambda0) that the PB design applies.

A constant phase or group-delay offset drops out of S. Atom selection uses the
CELERIS picker (a numpy copy, checked against design_achromatic_metalens or
design_pb_achromatic_metalens) and scans the free offsets: group-delay reference
for both kinds, and the base-phase offset for prop (PB hits any base phase by
rotation). D_S is the largest aperture whose band-minimum S stays >= 0.8.

Run from the repo root. One library, printed:
    python research/r1/coverage.py one tio2 sio2 1.0 --kind pb
Every library of a sweep config, one JSON each (skipped if it exists):
    python research/r1/coverage.py run --config full
    python research/r1/coverage.py run --config full --M 10 --only "^prop-tio2"
    python research/r1/coverage.py report --config full
"""
import argparse
import json
import math
import re
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import metric  # noqa: E402
import sweep  # noqa: E402

C = metric.C_UM_PER_FS
DIAMETERS_UM = (4, 8, 12, 16, 20, 26.4, 32, 40)
STREHL_CUT = 0.8  # Marechal criterion
WEIGHTS = ((0.0, "standard"), (1.0, "achromat"))  # gd_weight of the picker


class Library:
    """Arrays of one library: base phase, delay, mean amplitude, band response."""

    def __init__(self, lib, period_um: float, kind: str):
        a = lib.atoms
        self.kind = kind
        self.period = period_um
        self.ph0 = np.array([x.phase0_rad for x in a])
        self.gd = np.array([x.group_delay_fs for x in a])
        self.mamp = np.array([x.mean_amplitude for x in a])
        n_band = len(lib.wavelengths_um)
        if kind == "prop":
            self.U = np.array([[x.amplitude[j] * np.exp(1j * x.phase_rad[j])
                                for j in range(n_band)] for x in a])
        else:
            self.U = np.array([list(x.cross) for x in a])
        self.lam = np.array(lib.wavelengths_um)
        w = 2 * np.pi * C / self.lam
        self.dw_half = 0.5 * abs(w[0] - w[-1])
        self.k0 = 2 * np.pi / lib.center_wavelength_um
        self.gd_min, self.gd_max = lib.gd_min_fs, lib.gd_max_fs

    def sites(self, D: float):
        """Radii of the square-grid sites (as in the C++ designs) and the in-aperture mask."""
        n = max(1, round(D / self.period))
        c = (n - 1) / 2
        x = (np.arange(n) - c) * self.period
        X, Y = np.meshgrid(x, x)
        r = np.hypot(X, Y).ravel()
        return r, r <= D / 2

    def pick(self, sag, gd_ref, phi_off, gd_weight=1.0, amp_weight=0.25):
        """Per-site atom index; same cost as the C++ design for this kind."""
        tgd = -sag / C + gd_ref
        dgd = (self.gd[None, :] - tgd[:, None]) * self.dw_half
        cost = gd_weight * dgd**2 + amp_weight * (1 - self.mamp[None, :])**2
        if self.kind == "prop":
            tphi = -self.k0 * sag + phi_off
            cost = cost + np.angle(np.exp(1j * (self.ph0[None, :] - tphi[:, None])))**2
        return np.argmin(cost, axis=1)

    def strehl(self, idx, sag, inside):
        """On-axis Strehl at the design focus, one value per band wavelength."""
        u = self.U[idx]
        if self.kind == "pb":
            u = u * np.exp(1j * (-self.k0 * sag - self.ph0[idx]))[:, None]
        u, s = u[inside], sag[inside]
        k = 2 * np.pi / self.lam
        field = (u * np.exp(1j * np.outer(s, k))).sum(0)
        return np.abs(field)**2 / np.abs(u).sum(0)**2

    def best_strehl(self, D, focal, gd_weight, n_gd=13, n_phi=8):
        """Band Strehl of the best pick over the offset scan (by band-minimum S)."""
        r, inside = self.sites(D)
        sag = np.sqrt(r * r + focal * focal) - focal
        req = (math.hypot(D / 2, focal) - focal) / C
        phis = np.linspace(0, 2 * np.pi, n_phi, endpoint=False) if self.kind == "prop" else [0.0]
        best = None
        for gref in np.linspace(self.gd_min, self.gd_max + req, n_gd):
            for po in phis:
                s = self.strehl(self.pick(sag, gref, po, gd_weight), sag, inside)
                if best is None or s.min() > best.min():
                    best = s
        return req, best


def check_picker(cel, L: Library, lib, focal: float) -> None:
    """The numpy picker must reproduce the C++ one at the C++ offsets."""
    for D in (10.0, 25.0):
        r, _ = L.sites(D)
        sag = np.sqrt(r * r + focal * focal) - focal
        req = (math.hypot(D / 2, focal) - focal) / C
        gref = 0.5 * (L.gd_min + L.gd_max) + 0.5 * req
        if L.kind == "prop":
            d = cel.design_achromatic_metalens(lib, focal, D)
        else:
            d = cel.design_pb_achromatic_metalens(lib, focal, D)
        if not np.array_equal(L.pick(sag, gref, 0.0), np.array(d.atom_index)):
            raise RuntimeError(f"{L.kind} picker mismatch vs the C++ design at D={D}")


def coverage(cel, cfg: dict, entry: tuple, n_fill=None) -> dict:
    """Build one library and return its Strehl-vs-diameter table and D_S."""
    mat, sub, height, kind = entry
    lam_min, lam_max = cfg["band"]
    _, wc, _ = metric.band(lam_min, lam_max)
    lam_c = 2 * math.pi * C / wc
    pillar, substrate = cel.materials.by_name(mat), cel.materials.by_name(sub)
    P = sweep.choose_period_um(substrate, lam_min, cfg["period_margin"])
    fmin, fmax = sweep.fill_range(P, cfg["min_feature_um"], cfg["min_gap_um"])
    band = sweep.wavelengths(lam_min, lam_max, cfg["n_band"])

    t0 = time.time()
    lib, _ = sweep.build(cel, kind, pillar, substrate, P, band, lam_c, fmin, fmax,
                         n_fill or cfg["n_fill"], height, cfg["M"])
    L = Library(lib, P, kind)
    F = sweep.FOCAL_UM
    check_picker(cel, L, lib, F)

    span = float(L.gd.max() - L.gd.min())
    out = dict(material=mat, substrate=sub, height_um=height, kind=kind, M=cfg["M"],
               n_band=cfg["n_band"], n_fill=n_fill or cfg["n_fill"], period_um=P,
               focal_um=F, strehl_cut=STREHL_CUT, n_atoms=len(lib.atoms), span_fs=span,
               d_lib_um=metric.achievable_diameter_um(span, F), designs={})
    for gdw, label in WEIGHTS:
        rows, d_s = [], None
        for D in DIAMETERS_UM:
            req, s = L.best_strehl(D, F, gdw)
            rows.append(dict(diameter_um=D, required_fs=req, strehl_min=float(s.min()),
                             strehl=[float(x) for x in s]))
            if s.min() >= STREHL_CUT:
                d_s = D
        out["designs"][label] = dict(gd_weight=gdw, d_s_um=d_s, rows=rows)
    out["total_s"] = time.time() - t0
    return out


def print_one(res: dict) -> None:
    print(f"{res['kind']} {res['material']}/{res['substrate']} H={res['height_um']} "
          f"P={res['period_um']} M={res['M']} atoms={res['n_atoms']} "
          f"span={res['span_fs']:.2f} fs  D_lib={res['d_lib_um']:.1f} um")
    for label, d in res["designs"].items():
        for row in d["rows"]:
            print(f"  {label:8s} D {row['diameter_um']:5.1f} um  required "
                  f"{row['required_fs']:5.2f} fs  S_min {row['strehl_min']:.3f}")
        print(f"  {label}: D_S({STREHL_CUT}) >= {d['d_s_um']} um (grid {DIAMETERS_UM})")


def config_of(args) -> tuple[dict, Path]:
    """The sweep config (with any --M override) and its coverage result directory."""
    cfg = dict(sweep.CONFIGS[args.config])
    name = f"coverage-{args.config}"
    if args.M:
        cfg["M"] = args.M
        name += f"-M{args.M}"
    return cfg, Path(args.out) / name


def cmd_one(args) -> int:
    import celeris as cel
    cfg, _ = config_of(args)
    print_one(coverage(cel, cfg, (args.material, args.substrate, args.height_um, args.kind),
                       n_fill=args.n_fill))
    return 0


def cmd_run(args) -> int:
    import celeris as cel
    if not hasattr(cel.DispersivePbAtom, "cross"):
        print("this celeris build does not expose DispersivePbAtom.cross; rebuild from main",
              file=sys.stderr)
        return 2
    cfg, out = config_of(args)
    out.mkdir(parents=True, exist_ok=True)
    for entry in cfg["libraries"]:
        key = sweep.key_of(entry)
        if args.only and not re.search(args.only, key):
            continue
        path = out / (key + ".json")
        if path.exists() and not args.force:
            print(f"skip {path.name} (exists)")
            continue
        res = coverage(cel, cfg, entry)
        path.write_text(json.dumps(res, indent=1))
        ds = {k: v["d_s_um"] for k, v in res["designs"].items()}
        print(f"{path.name}: span {res['span_fs']:.2f} fs, D_lib {res['d_lib_um']:.1f} um, "
              f"D_S {ds}  [{res['total_s']:.0f} s]", flush=True)
    return 0


def cmd_report(args) -> int:
    _, out = config_of(args)
    files = sorted(out.glob("*.json"))
    if not files:
        print(f"no results in {out}")
        return 1
    print(f"{'library':32s} {'span fs':>8s} {'D_lib um':>9s} {'D_S std':>8s} {'D_S achr':>9s}")
    for f in files:
        r = json.loads(f.read_text())
        d = r["designs"]
        print(f"{f.stem:32s} {r['span_fs']:8.2f} {r['d_lib_um']:9.1f} "
              f"{str(d['standard']['d_s_um']):>8s} {str(d['achromat']['d_s_um']):>9s}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--config", default="full", choices=sorted(sweep.CONFIGS))
    common.add_argument("--M", type=int, default=None,
                        help="override the config's Fourier order")
    common.add_argument("--out", default=str(HERE / "out"))

    p = sub.add_parser("one", parents=[common], help="one library, printed")
    p.add_argument("material")
    p.add_argument("substrate")
    p.add_argument("height_um", type=float)
    p.add_argument("--kind", default="prop", choices=("prop", "pb"))
    p.add_argument("--n-fill", type=int, default=None,
                   help="override the config's fill samples (sampling-density control)")
    p.set_defaults(fn=cmd_one)

    p = sub.add_parser("run", parents=[common], help="every library of a config, to JSON")
    p.add_argument("--only", default=None, help="regex on the library key")
    p.add_argument("--force", action="store_true", help="recompute existing results")
    p.set_defaults(fn=cmd_run)

    p = sub.add_parser("report", parents=[common], help="table of D_S per library")
    p.set_defaults(fn=cmd_report)

    args = ap.parse_args()
    return args.fn(args)


if __name__ == "__main__":
    sys.exit(main())
