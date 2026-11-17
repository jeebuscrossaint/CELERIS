"""Check the Python material registry against the shipped data (issue #25).

Run from the repo root (so data/ is found):  python tests/python/check_materials.py
Plain asserts, no pytest needed. Exits non-zero on the first failure.
"""
from pathlib import Path

import celeris as cel

m = cel.materials
REPO = Path(__file__).resolve().parents[2]

# Every canonical entry is listed, and every tabulated file is found.
names = {i.name for i in m.catalog()}
for want in ("air", "sio2", "bk7", "si3n4", "al2o3", "gan", "tio2", "a-si", "c-si"):
    assert want in names, f"{want} missing from catalog: {sorted(names)}"
for info in m.catalog():
    if info.tabulated:
        assert info.available, f"{info.data_file} not found"
        assert m.find_data_file(info.data_file), info.data_file

# by_name('tio2') reproduces data/TiO2_Siefke.csv exactly at a table row.
rows = [l.split() for l in (REPO / "data" / "TiO2_Siefke.csv").read_text().splitlines()
        if l.strip() and not l.startswith("#")]
lam, n, k = (float(x) for x in rows[len(rows) // 2])
got = m.by_name("TiO2").index(lam)
assert abs(got.real - n) < 1e-9 and abs(got.imag - k) < 1e-12, (lam, got, n, k)

# Aliases / case / '-' vs '_' resolve to the same entry.
assert m.by_name("A_Si").index(0.62) == m.by_name("a-si").index(0.62)

# Analytic models: published values at 550 nm (lossless).
sap = m.sapphire().index(0.55)
gan = m.gallium_nitride().index(0.55)
assert abs(sap.real - 1.7704) < 2e-3 and sap.imag == 0.0, sap  # Malitson-Dodge n_o
assert 2.35 < gan.real < 2.45 and gan.imag == 0.0, gan        # alpha-GaN n_o ~2.41

# Unknown names raise.
try:
    m.by_name("unobtainium")
except RuntimeError:
    pass
else:
    raise AssertionError("by_name accepted an unknown material")

print(f"material registry OK: {len(names)} entries; tio2 n({lam:.4f})={got.real:.6f}, "
      f"sapphire {sap.real:.4f}, GaN {gan.real:.4f}")
