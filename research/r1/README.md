# R1: how much of the achromatic limit a single-etch library reaches

This directory holds the scripts for the R1 study (ROADMAP §4): the delay span that a
single-etch meta-atom library supplies, as a fraction of the Presutti–Monticone bound.

| file | what it is |
|---|---|
| [`metric.md`](metric.md) | Definitions: delay span, η_DL, η_PM, D_lib, and why. |
| [`metric.py`](metric.py) | The formulas (stdlib only). `python3 metric.py` runs self-checks. |
| [`sweep.py`](sweep.py) | The sweep driver: builds libraries with CELERIS and applies the metric. |

## Running the sweep

```sh
pip install .                                   # the celeris Python package, from the repo root
python research/r1/sweep.py run --config smoke  # plumbing check, seconds
python research/r1/sweep.py run --config pilot  # TiO2 on fused silica, two heights
python research/r1/sweep.py report --config pilot
```

Results go to `research/r1/out/<config>/` (one JSON per library; not committed).
A re-run skips libraries that already have a result, so a long sweep can be resumed.

## Choices the driver makes (and why)

- **Period.** The largest 10 nm step below 0.95 · λ_min / n_sub. This keeps every
  substrate diffraction order evanescent across the band. A 400 nm cell on glass
  is *not* subwavelength at 470 nm: up to 25 % of the incident power goes into
  substrate orders (issue #30).
- **Fill range.** Limited by an assumed minimum feature and minimum gap (60 nm by
  default). These are round numbers, not values from a specific process, and every
  output file records them.
- **Libraries.** `prop` uses square and circular pillars, propagation phase, and mean |t|
  as the efficiency. `pb` uses a rectangle (fill_x × fill_y) grid with Pancharatnam–Berry
  phase and mean |a_cross| as the efficiency. Both use CELERIS's default quality filters
  (#20, #27).
- **Efficiency is an axis.** The span is reported at each threshold in
  `THRESHOLDS` (metric.md §4), not at one cut.
- **Band-sampling convergence.** Each library is rebuilt with 2N − 1 band points, and
  the change in span is reported. Coarse bands alias resonant atoms (#21). Treat a large
  change the way an unconverged M would be treated.
- **Fourier-order convergence.** `--M` overrides the config's Fourier order and writes
  to `out/<config>-M<M>/`. The `mcheck` config holds the libraries for this check (#34):
  ```sh
  python research/r1/sweep.py run --config mcheck --M 8 --no-converge
  python research/r1/sweep.py run --config mcheck --M 10 --no-converge --only "^prop-|pb-tio2-on-sio2-H1.00"
  python research/r1/sweep.py compare full mcheck-M8
  ```
- **Red flag.** η_PM > 1 is marked in the report. It means a numerical problem or a
  resonator outside Tucker's model, never a result to quote as is.

The driver needs the refined band unwrap (#21). It exits with an error on an older
build.
