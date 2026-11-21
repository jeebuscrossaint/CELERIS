"""R1 sweep driver: delay span and fraction of the achromatic limit, per library.

For every (material, substrate, height, library kind) in a config this script:

1. picks the lattice period so the cell is subwavelength in the substrate at the
   shortest band wavelength (no propagating substrate orders; issue #30),
2. limits the fill range by an assumed minimum feature and minimum gap,
3. builds the dispersive library with CELERIS (propagation-phase squares + circles,
   or a Pancharatnam-Berry rectangle grid) with N band points,
4. reports the delay span as a curve against an efficiency threshold (metric.md §4),
   with eta_DL, eta_PM and D_lib from metric.py,
5. rebuilds the library with 2N-1 band points and reports how much each span moves
   (band-sampling convergence; issue #21 showed coarse bands alias).

Usage (from the repo root, with the celeris Python package installed):

    python research/r1/sweep.py run --config smoke
    python research/r1/sweep.py run --config pilot
    python research/r1/sweep.py report --config pilot

Each library is written to research/r1/out/<config>/<key>.json and skipped on a
re-run if the file exists, so a long sweep can be interrupted and resumed.

Results are only valid with the refined band unwrap (#21, PR #29); the script
refuses to run against a celeris build without `n_refined_samples`.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import metric  # noqa: E402  (same directory)

# Lens used to turn a delay span into a diameter: the Chen-2018 achromat
# (NA 0.20, D 26.4 um -> f = 64.67 um). Fixed so every row is comparable.
FOCAL_UM = 13.2 / math.tan(math.asin(0.20))

# Efficiency thresholds for the span-vs-efficiency curve. For PB atoms the
# amplitude is mean |a_cross|; for propagation-phase atoms it is mean |t|.
THRESHOLDS = (0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90)

# Configs. Fabrication limits are ASSUMED round numbers, not taken from a process;
# they are recorded in every output file and can be changed per config.
CONFIGS = {
    # Plumbing check: seconds, not physics.
    "smoke": dict(
        band=(0.470, 0.670), n_band=5, M=3, n_fill=4, min_feature_um=0.06,
        min_gap_um=0.06, period_margin=0.95,
        libraries=[("tio2", "sio2", 0.6, "prop"), ("tio2", "sio2", 0.6, "pb")],
    ),
    # First real look: TiO2 on fused silica, two heights, both kinds.
    "pilot": dict(
        band=(0.470, 0.670), n_band=7, M=6, n_fill=10, min_feature_um=0.06,
        min_gap_um=0.06, period_margin=0.95,
        libraries=[(m, s, h, k) for (m, s) in [("tio2", "sio2")]
                   for h in (0.6, 1.0) for k in ("prop", "pb")],
    ),
    # The study grid (materials x heights x kinds). Heights stop where the aspect
    # ratio height / min_feature passes ~25.
    "full": dict(
        band=(0.470, 0.670), n_band=7, M=6, n_fill=12, min_feature_um=0.06,
        min_gap_um=0.06, period_margin=0.95,
        libraries=[(m, s, h, k)
                   for (m, s) in [("tio2", "sio2"), ("si3n4", "sio2"), ("gan", "al2o3")]
                   for h in (0.4, 0.6, 0.8, 1.0, 1.2, 1.5)
                   for k in ("prop", "pb")],
    ),
    # Fourier-order (M) convergence check (#34): run with --M 8 / --M 10 and compare
    # against "full" (M = 6). H = 1.0 per material and kind, plus the worst
    # band-convergence case (TiO2 PB H = 0.8).
    "mcheck": dict(
        band=(0.470, 0.670), n_band=7, M=6, n_fill=12, min_feature_um=0.06,
        min_gap_um=0.06, period_margin=0.95,
        libraries=[(m, s, 1.0, k)
                   for (m, s) in [("tio2", "sio2"), ("si3n4", "sio2"), ("gan", "al2o3")]
                   for k in ("prop", "pb")] + [("tio2", "sio2", 0.8, "pb")],
    ),
}


def wavelengths(lo: float, hi: float, n: int) -> list[float]:
    """n wavelengths evenly spaced in lambda over [lo, hi] (um)."""
    return [lo + (hi - lo) * i / (n - 1) for i in range(n)]


def choose_period_um(substrate, lam_min: float, margin: float) -> float:
    """Largest 10-nm-step period with no propagating substrate or air order at lam_min.

    A first substrate order propagates when lambda < P * n_sub (normal incidence),
    and an air order when lambda < P. Take margin * lam_min / max(n_sub, 1) and
    round DOWN to 10 nm.
    """
    n_sub = complex(substrate.index(lam_min)).real
    p = margin * lam_min / max(n_sub, 1.0)
    return math.floor(p * 100.0 + 1e-9) / 100.0


def fill_range(period_um: float, min_feature_um: float, min_gap_um: float) -> tuple[float, float]:
    """(fill_min, fill_max) allowed by the minimum feature and the minimum gap."""
    return min_feature_um / period_um, 1.0 - min_gap_um / period_um


def build(cel, kind: str, pillar, substrate, period: float, band: list[float],
          lam_c: float, fmin: float, fmax: float, n_fill: int, height: float, M: int):
    """Build one library; returns (library, per-atom amplitude used for the curve)."""
    air = cel.materials.air()
    if kind == "prop":
        fills = [fmin + (fmax - fmin) * i / (n_fill - 1) for i in range(n_fill)]
        specs = [cel.MetaAtomSpec(shape, f, f, height, 0.5)
                 for f in fills for shape in (cel.MetaShape.Rectangle, cel.MetaShape.Ellipse)]
        lib = cel.build_dispersive_library_from_specs(
            pillar, air, air, substrate, period, band, lam_c, specs, M)
        return lib, [a.mean_amplitude for a in lib.atoms]
    if kind == "pb":
        lib = cel.build_dispersive_pb_library(
            pillar=pillar, background=air, incident=air, substrate=substrate,
            period_um=period, band_wavelengths_um=band, center_wavelength_um=lam_c,
            fill_min=fmin, fill_max=fmax, n_fills=n_fill, thickness_um=height, M=M)
        return lib, [a.mean_amplitude for a in lib.atoms]
    raise ValueError(f"unknown library kind {kind!r}")


def atom_key(kind: str, a) -> tuple:
    """Geometry key used to match atoms between the N and 2N-1 builds."""
    if kind == "prop":
        return (str(a.shape), round(a.fill_x, 9))
    return (round(a.fill_x, 9), round(a.fill_y, 9))


def span_curve(gd: list[float], amp: list[float]) -> list[dict]:
    """Delay span among atoms with amplitude >= each threshold (None if < 2 atoms)."""
    out = []
    for t in THRESHOLDS:
        sel = [g for g, a in zip(gd, amp) if a >= t]
        out.append(dict(threshold=t, n_atoms=len(sel),
                        span_fs=(max(sel) - min(sel)) if len(sel) >= 2 else None))
    return out


def run_library(cel, cfg: dict, entry: tuple, converge: bool) -> dict:
    mat, sub, height, kind = entry
    lam_min, lam_max = cfg["band"]
    _, wc, _ = metric.band(lam_min, lam_max)
    lam_c = 2.0 * math.pi * metric.C_UM_PER_FS / wc
    pillar, substrate = cel.materials.by_name(mat), cel.materials.by_name(sub)
    period = choose_period_um(substrate, lam_min, cfg["period_margin"])
    fmin, fmax = fill_range(period, cfg["min_feature_um"], cfg["min_gap_um"])
    n_max = complex(pillar.index(lam_c)).real
    n_sub = complex(substrate.index(lam_min)).real

    t0 = time.time()
    band = wavelengths(lam_min, lam_max, cfg["n_band"])
    lib, amp = build(cel, kind, pillar, substrate, period, band, lam_c, fmin, fmax,
                     cfg["n_fill"], height, cfg["M"])
    gd = [a.group_delay_fs for a in lib.atoms]
    curve = span_curve(gd, amp)
    t_build = time.time() - t0

    conv = None
    if converge:
        band2 = wavelengths(lam_min, lam_max, 2 * cfg["n_band"] - 1)
        lib2, _ = build(cel, kind, pillar, substrate, period, band2, lam_c, fmin, fmax,
                        cfg["n_fill"], height, cfg["M"])
        gd2 = {atom_key(kind, a): a.group_delay_fs for a in lib2.atoms}
        # Span of the SAME kept atoms with their 2N-1 delays (an atom the finer
        # library drops is reported, not silently matched).
        matched = [gd2.get(atom_key(kind, a)) for a in lib.atoms]
        missing = sum(g is None for g in matched)
        diffs = [abs(g2 - g) for g, g2 in zip(gd, matched) if g2 is not None]
        curve2 = span_curve([g2 if g2 is not None else g for g, g2 in zip(gd, matched)], amp)
        conv = dict(n_band=len(band2), max_atom_gd_change_fs=max(diffs) if diffs else None,
                    atoms_dropped_by_finer_band=missing,
                    span_change_fs=[None if (c1["span_fs"] is None or c2["span_fs"] is None)
                                    else c2["span_fs"] - c1["span_fs"]
                                    for c1, c2 in zip(curve, curve2)])

    for c in curve:
        s = c["span_fs"]
        if s is None:
            continue
        c["eta_dl"] = metric.eta_dl(s, height, n_max, 1.0)
        c["eta_pm"] = metric.eta_pm(s, height, n_max, 1.0, lam_min, lam_max)
        c["d_lib_um"] = metric.achievable_diameter_um(s, FOCAL_UM)

    return dict(
        material=mat, substrate=sub, height_um=height, kind=kind,
        band_um=[lam_min, lam_max], n_band=cfg["n_band"], M=cfg["M"],
        period_um=period, period_over_cutoff=period * n_sub / lam_min,
        fill_range=[fmin, fmax], min_feature_um=cfg["min_feature_um"],
        min_gap_um=cfg["min_gap_um"], n_fill=cfg["n_fill"],
        n_max=n_max, lam_c_um=lam_c, focal_um=FOCAL_UM,
        delay_line_fs=metric.delay_line_fs(height, n_max, 1.0),
        tucker_bound_fs=metric.tucker_bound_fs(height, n_max, 1.0, lam_min, lam_max),
        n_atoms=len(lib.atoms), refined_samples=sum(a.n_refined_samples for a in lib.atoms),
        atoms=[dict(gd_fs=a.group_delay_fs, amp=m, fill_x=a.fill_x, fill_y=a.fill_y,
                    fit_rms_rad=a.phase_fit_rms_rad, n_refined=a.n_refined_samples)
               for a, m in zip(lib.atoms, amp)],
        curve=curve, convergence=conv, build_s=t_build, total_s=time.time() - t0,
    )


def key_of(entry: tuple) -> str:
    mat, sub, h, kind = entry
    return f"{kind}-{mat}-on-{sub}-H{h:.2f}"


def config_of(args) -> tuple[dict, Path]:
    """The config (with any --M override) and its result directory."""
    cfg = dict(CONFIGS[args.config])
    name = args.config
    if getattr(args, "M", None):
        cfg["M"] = args.M
        name = f"{args.config}-M{args.M}"
    return cfg, Path(args.out) / name


def cmd_run(args) -> int:
    import celeris as cel
    if not hasattr(cel.DispersiveAtom, "n_refined_samples"):
        print("this celeris build has no refined band unwrap (#21); rebuild from main",
              file=sys.stderr)
        return 2
    cfg, out = config_of(args)
    out.mkdir(parents=True, exist_ok=True)
    for entry in cfg["libraries"]:
        path = out / (key_of(entry) + ".json")
        if path.exists() and not args.force:
            print(f"skip {path.name} (exists)")
            continue
        res = run_library(cel, cfg, entry, converge=not args.no_converge)
        path.write_text(json.dumps(res, indent=1))
        c0 = res["curve"][0]
        print(f"{path.name}: P={res['period_um']*1000:.0f} nm, {res['n_atoms']} atoms, "
              f"span@{c0['threshold']:.2f} = {c0['span_fs'] or float('nan'):.2f} fs, "
              f"eta_PM {c0.get('eta_pm', float('nan')):.3f}  [{res['total_s']:.0f} s]",
              flush=True)
    return 0


def cmd_report(args) -> int:
    _, out = config_of(args)
    rows = [json.loads(p.read_text()) for p in sorted(out.glob("*.json"))]
    if not rows:
        print(f"no results in {out}")
        return 1
    print(f"{'library':32s} {'P nm':>5s} {'atoms':>5s} {'thr':>4s} {'span fs':>8s} "
          f"{'eta_DL':>7s} {'eta_PM':>7s} {'D_lib um':>8s} {'dspan(2N-1)':>11s}")
    for r in rows:
        name = key_of((r["material"], r["substrate"], r["height_um"], r["kind"]))
        conv = r["convergence"]
        for i, c in enumerate(r["curve"]):
            if c["span_fs"] is None:
                continue
            ds = conv["span_change_fs"][i] if conv else None
            flag = "  <-- eta_PM > 1: investigate" if c.get("eta_pm", 0) > 1 else ""
            print(f"{name if i == 0 else '':32s} {r['period_um']*1000:5.0f} {c['n_atoms']:5d} "
                  f"{c['threshold']:4.2f} {c['span_fs']:8.2f} {c['eta_dl']:7.3f} "
                  f"{c['eta_pm']:7.3f} {c['d_lib_um']:8.1f} "
                  f"{'' if ds is None else f'{ds:+.2f}':>11s}{flag}")
    return 0


def cmd_compare(args) -> int:
    """Span change of each library in result dir B against the same library in A."""
    a_dir, b_dir = Path(args.out) / args.a, Path(args.out) / args.b
    rows = 0
    print(f"{'library':32s} {'thr':>4s} {args.a + ' fs':>12s} {args.b + ' fs':>14s} {'change':>7s}")
    for pb in sorted(b_dir.glob("*.json")):
        pa = a_dir / pb.name
        if not pa.exists():
            print(f"{pb.stem:32s} (not in {args.a})")
            continue
        ra, rb = json.loads(pa.read_text()), json.loads(pb.read_text())
        for i, (ca, cb) in enumerate(zip(ra["curve"], rb["curve"])):
            if ca["span_fs"] is None or cb["span_fs"] is None:
                continue
            print(f"{pb.stem if i == 0 else '':32s} {ca['threshold']:4.2f} "
                  f"{ca['span_fs']:12.2f} {cb['span_fs']:14.2f} "
                  f"{cb['span_fs'] - ca['span_fs']:+7.2f}")
            rows += 1
    return 0 if rows else 1


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("run", "report"):
        p = sub.add_parser(name)
        p.add_argument("--config", choices=sorted(CONFIGS), default="smoke")
        p.add_argument("--out", default=str(HERE / "out"))
        p.add_argument("--M", type=int, default=None,
                       help="override the config's Fourier order (results go to <config>-M<M>)")
        if name == "run":
            p.add_argument("--force", action="store_true", help="recompute existing results")
            p.add_argument("--no-converge", action="store_true",
                           help="skip the 2N-1 band-point rebuild")
    p = sub.add_parser("compare", help="span change between two result directories")
    p.add_argument("a", help="reference result directory name, e.g. full")
    p.add_argument("b", help="result directory to compare, e.g. mcheck-M8")
    p.add_argument("--out", default=str(HERE / "out"))
    args = ap.parse_args()
    return {"run": cmd_run, "report": cmd_report, "compare": cmd_compare}[args.cmd](args)


if __name__ == "__main__":
    sys.exit(main())
