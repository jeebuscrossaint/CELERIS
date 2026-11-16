"""R1 metric: how much of the achromatic limit a single-etch meta-atom library reaches.

The derivation and the reasoning behind each definition are in metric.md (same
directory). This module only implements the formulas so that every R1 script
uses the same definitions. Running it executes numeric self-checks and prints
the metric for the libraries that already exist in the repo.

Units: lengths in micrometres, times in femtoseconds, angular frequencies in
rad/fs. With these, omega = 2*pi*c/lambda with c = 0.299792458 um/fs.

Pure Python (stdlib only): no RCWA is run here.
"""
from __future__ import annotations

import math
from pathlib import Path

C_UM_PER_FS = 0.299792458  # speed of light [um/fs]


# --------------------------------------------------------------------------- band
def omega(lam_um: float) -> float:
    """Angular frequency [rad/fs] of a vacuum wavelength [um]."""
    return 2.0 * math.pi * C_UM_PER_FS / lam_um


def band(lam_min_um: float, lam_max_um: float) -> tuple[float, float, float]:
    """Return (delta_omega, omega_c, fractional bandwidth) for a wavelength band.

    omega_c is the centre in FREQUENCY, (omega_min + omega_max) / 2, which is
    what the delay-bandwidth bound uses. It is not the omega of the
    wavelength-midpoint.
    """
    w_hi, w_lo = omega(lam_min_um), omega(lam_max_um)
    dw = w_hi - w_lo
    wc = 0.5 * (w_hi + w_lo)
    return dw, wc, dw / wc


# ----------------------------------------------------------------- lens geometry
def required_delay_fs(diameter_um: float, focal_um: float) -> float:
    """Edge-minus-centre group delay an achromatic lens needs [fs].

    Presutti & Monticone Eq. 4: dT = (sqrt(R^2 + F^2) - F) / c. Same quantity as
    `required_gd_span_fs` in CELERIS's design_pb_achromatic_metalens.
    """
    r = 0.5 * diameter_um
    return (math.hypot(r, focal_um) - focal_um) / C_UM_PER_FS


def achievable_diameter_um(delay_span_fs: float, focal_um: float) -> float:
    """Largest achromatic diameter a library with this group-delay span supports [um].

    Inverse of required_delay_fs. Identical to Buisine et al. Eq. 15,
    D = 2 f sqrt([1 + (c/f) d_beta1]^2 - 1), with d_beta1 = delay span.
    """
    x = C_UM_PER_FS * delay_span_fs
    return 2.0 * focal_um * math.sqrt((1.0 + x / focal_um) ** 2 - 1.0)


def buisine_eq7_diameter_um(path_spread_um: float, na: float) -> float:
    """Buisine et al. Eq. 7, D = 2 h(dn - dn0) NA / (1 - sqrt(1 - NA^2)).

    `path_spread_um` is h(dn - dn0). metric.md shows this equals Eq. 15 once
    h(dn - dn0) = c * delay span and NA = R / sqrt(R^2 + f^2) (lens in air).
    """
    return 2.0 * path_spread_um * na / (1.0 - math.sqrt(1.0 - na * na))


# ------------------------------------------------------------------- the limits
def delay_line_fs(height_um: float, n_max: float, n_min: float = 1.0) -> float:
    """Non-dispersive delay-line delay spread h (n_max - n_min) / c [fs].

    The delay difference between light crossing height h entirely in the
    pillar material and entirely in the background. It is a reference scale,
    not a strict bound: guided-mode dispersion can exceed it (eta_dl > 1).
    """
    return height_um * (n_max - n_min) / C_UM_PER_FS


def tucker_bound_fs(height_um: float, n_max: float, n_min: float,
                    lam_min_um: float, lam_max_um: float) -> float:
    """Presutti & Monticone delay span bound with Tucker's kappa [fs].

    Eq. 3 dT * d_omega <= kappa with Eq. 8 kappa = 2 pi (L / lambda_c)(n_max - n_min)
    = omega_c L dn / c, so dT_max = delay_line * (omega_c / d_omega). Valid for
    fractional bandwidth < 1 (Tucker's assumption).
    """
    dw, wc, fbw = band(lam_min_um, lam_max_um)
    if fbw >= 1.0:
        raise ValueError(f"Tucker bound needs fractional bandwidth < 1 (got {fbw:.3f})")
    return delay_line_fs(height_um, n_max, n_min) * wc / dw


def eta_dl(delay_span_fs: float, height_um: float, n_max: float, n_min: float = 1.0) -> float:
    """Fraction of the delay-line scale a library reaches (can exceed 1)."""
    return delay_span_fs / delay_line_fs(height_um, n_max, n_min)


def eta_pm(delay_span_fs: float, height_um: float, n_max: float, n_min: float,
           lam_min_um: float, lam_max_um: float) -> float:
    """Fraction of the Presutti-Monticone (Tucker) bound a library reaches.

    Must be <= 1. A value above 1 means a numerical error (e.g. an aliased or
    slipped group delay, see issues #20/#21) or a meta-atom outside Tucker's
    model (a resonator). Treat it as a red flag, not a result.
    """
    return delay_span_fs / tucker_bound_fs(height_um, n_max, n_min, lam_min_um, lam_max_um)


# ---------------------------------------------------------------- material data
def tabulated_n(csv_path: Path, lam_um: float) -> float:
    """Linear interpolation of n from a CELERIS data/*.csv (lambda_um n k)."""
    rows = []
    for line in csv_path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.replace(",", " ").split()
        try:
            rows.append((float(parts[0]), float(parts[1])))
        except (ValueError, IndexError):
            continue  # header line
    rows.sort()
    for (l0, n0), (l1, n1) in zip(rows, rows[1:]):
        if l0 <= lam_um <= l1:
            return n0 + (n1 - n0) * (lam_um - l0) / (l1 - l0)
    raise ValueError(f"{lam_um} um outside {csv_path.name}")


# ------------------------------------------------------------------ self-checks
def _self_check() -> None:
    # Eq. 7 == Eq. 15 (== inverse of P&M Eq. 4) for a range of geometries.
    for f in (20.0, 64.67, 500.0):
        for dtau in (0.5, 2.6, 6.0, 20.0):
            d15 = achievable_diameter_um(dtau, f)
            r = 0.5 * d15
            na = r / math.hypot(r, f)
            d7 = buisine_eq7_diameter_um(C_UM_PER_FS * dtau, na)
            assert abs(d7 - d15) < 1e-9 * d15, (f, dtau, d7, d15)
            assert abs(required_delay_fs(d15, f) - dtau) < 1e-9 * dtau
    # Chen 2018: NA 0.20, D 26.4 um -> ~4.45 fs (CELERIS selftest [18]).
    f = 13.2 / math.tan(math.asin(0.20))
    assert abs(required_delay_fs(26.4, f) - 4.45) < 0.01
    # Tucker bound = delay line / fractional bandwidth.
    _, _, fbw = band(0.470, 0.670)
    assert abs(tucker_bound_fs(1, 2.4, 1, 0.47, 0.67) * fbw - delay_line_fs(1, 2.4, 1)) < 1e-12


def _report() -> None:
    lo, hi = 0.470, 0.670
    dw, wc, fbw = band(lo, hi)
    lam_c = 2 * math.pi * C_UM_PER_FS / wc
    print(f"band {lo*1e3:.0f}-{hi*1e3:.0f} nm: d_omega {dw:.4f} rad/fs, omega_c {wc:.4f} rad/fs "
          f"(lambda_c {lam_c*1e3:.1f} nm), fractional bandwidth {fbw:.4f}")
    f_chen = 13.2 / math.tan(math.asin(0.20))
    print(f"Chen 2018 lens (D 26.4 um, f {f_chen:.2f} um): required delay "
          f"{required_delay_fs(26.4, f_chen):.2f} fs\n")

    n_ald = tabulated_n(Path(__file__).resolve().parents[2] / "data" / "TiO2_Siefke.csv", lam_c)
    # (label, measured span [fs], height [um], n_pillar used for the limits)
    # Spans: PR #22 / issue #20 (2026-11-16). n=2.4 rows are constant-index models.
    libs = [
        ("pilot n=2.4 H=1.2 (7 pts, filtered)", 5.97, 1.2, 2.40),
        ("fig_achromatic n=2.4 H=0.6 (7 pts)", 2.58, 0.6, 2.40),
        ("selftest [18] n=2.4 H=0.6 (3 pts)", 2.21, 0.6, 2.40),
        (f"reproduce chen2018 ALD H=0.6 (n(lc)={n_ald:.3f})", 3.11, 0.6, n_ald),
    ]
    print(f"{'library':48s} {'span':>6s} {'DL':>6s} {'PM':>6s} {'eta_DL':>7s} {'eta_PM':>7s} "
          f"{'D_lib':>7s} {'D_PM':>7s}  (fs, fs, fs, -, -, um, um at f={f_chen:.1f} um)")
    for label, span, h, n in libs:
        dl = delay_line_fs(h, n)
        pm = tucker_bound_fs(h, n, 1.0, lo, hi)
        print(f"{label:48s} {span:6.2f} {dl:6.2f} {pm:6.2f} {eta_dl(span, h, n):7.3f} "
              f"{eta_pm(span, h, n, 1.0, lo, hi):7.3f} {achievable_diameter_um(span, f_chen):7.1f} "
              f"{achievable_diameter_um(pm, f_chen):7.1f}")


if __name__ == "__main__":
    _self_check()
    print("self-checks passed\n")
    _report()
