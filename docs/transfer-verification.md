# Transfer export validation — 3 October 2026

The [schema and physics limits](transfer-export.md) describe what these checks
establish. The kerr-sbi production binary and recorded research artifacts were
not modified. All renders used a new local binary and fresh scratch directories.

## Commands and simulation identity

```sh
cargo test --workspace
cargo test -p nullgeo-cli --no-default-features --test transfer_cli
cargo clippy --workspace --all-targets -- -D warnings
cargo fmt --all -- --check
cargo build --release -p nullgeo-cli --target-dir target/transfer-v1
/Users/jameswirth/VSCodeProjects/kerr-sbi/.venv/bin/python \
  examples/verify_transfer.py --binary target/transfer-v1/release/nullgeo \
  --output scratch/transfer-verification-final
/Users/jameswirth/VSCodeProjects/kerr-sbi/.venv/bin/python -B \
  examples/transfer_visibility.py scratch/transfer-verification-final/grid \
  --q 2.8 --g-power 4
```

All checks passed, including a repeated CLI integration test and Clippy after the
final validation change. Python used NumPy 2.5.3. The consumer virtual environment
was only used to execute Python with its existing NumPy installation.

* Binary: `target/transfer-v1/release/nullgeo` (release, parallel).
* SHA-256: `e972cea4c2e264affb7e770e6367d3a0e944b3aa27d8401aef24a59fc7cbfccf`.
* Build ID: `transfer-v1-1791022257915363000`.
* Base revision: `c5a5cca50c8670dad3ee252ef06a3239d388a1c9`, **dirty**, with this implementation.
* Compiler: `rustc 1.89.0 (29483883e 2025-08-04)`, `aarch64-apple-darwin`.
* Machine-readable results: `scratch/transfer-verification-final/verification.json`.
  Scratch is ignored by Git; the results below are the tracked summary.

## Reconstruction and coverage

Scene: 32×24 pixels, Kerr mass=1 and spin=0.7, camera at Cartesian radius 85,
24-degree horizontal FOV, disk [8,18], tol=1e-8, escape radius=340,
max_steps=20000. Original emission: q=2, g_power=3. Independent ordinary renders
used (q,g_power)=(2.8,3), (2,4.2), and (3.1,1.7), at identical geometry and sampling.

| Mode | Contributing samples | Samples per pixel | Export directory bytes |
|---|---:|---|---:|
| Single | 768 | 1 | 81,632 |
| Single, jitter | 768 | 1 | 81,644 |
| 3×3 grid | 6,912 | 9 | 714,461 |
| 3×3 jitter | 6,912 | 9 | 714,463 |
| Adaptive 1→3 | 3,216 | 1 or 9 | 333,781 |
| Adaptive 2→3, jitter | 4,622 | 4 or 9 | 478,611 |

All six original reconstructions and all eighteen changed-emissivity comparisons
had **max absolute error=0 and relative L2 error=0**, after reproducing float32
casting, weighting and accumulation. Stored-radiance reconstruction was also
exact. This platform's observed equality is stronger than the documented f32
rounding tolerance; it does not promise bit identity across math libraries.
Ordinary PFM bytes were identical with and without export in every mode.

Sample counts, ordering, offsets, screen orientation, per-pixel weight sums,
NaN sentinels, masks and outcomes passed. Both adaptive cases included base and
replacement pixels. Ordinary scenes had no unfinished rays. Separate renders
with max_steps=0 and tol=1e-300 exercised actual max-step and stalled outcomes;
all 768 rays were correctly marked unfinished in each respective case. Unit
tests covered all six outcome codes and absent/nonfinite radius/redshift masks.
Unsupported disk/background/metric settings, nonfinite mass, invalid sampling,
and attempts to overwrite an export were rejected.

The direct visibility example used an assigned 125 microarcsecond horizontal
FOV and 0.6 Jy flux. Across modes, zero-baseline flux differed from 0.6 Jy by at
most 1.11e-16 Jy; conjugate-symmetry error was 0. At q=2.8, g_power=4 on the grid,
V(2e9,-1e9)=(-0.146844665816402+0.164085048341635i) Jy. This is a quadrature example,
not an observer-at-infinity prediction or visibility convergence result.

## Overhead and remaining work

Seven alternating, warmed process-level measurements of the 32×24, 3×3 grid:

* Ordinary median: **53.250 ms**, range 52.292–60.013 ms.
* Export enabled median: **53.978 ms**, range 53.464–57.151 ms.
* Observed incremental median: **0.727 ms (1.37%)**. Timing noise overlaps;
  this is a small-scene measurement, not a general performance guarantee.
* Extra output: **714,461 bytes** (697.7 KiB), versus **9,230 bytes** for the
  unchanged RGB PFM. Records cost 103 bytes per contributing subray, plus headers
  and metadata whose size depends on scene text and paths.

Both timings include process startup, tracing, shading and PFM writes; export
also includes validation, serialization and file writes. Neither requests fsync.
There is no second trace or retained duplicate sample array for export. With the
flag absent, the only new render-path work is checking the optional argument;
the camera helper preserves the existing arithmetic.

Direct subray visibility integration is now possible at fixed disk geometry.
Remaining scientific work is convergence in the actual baseline/noise norm,
finite-observer/angular-scale modeling, and treatment of any unfinished rays.
Changing an emitting edge requires ordered crossings beyond the first opaque
hit, including crossings outside the original annulus; emission inside the ISCO
requires a plunging-flow velocity model. These are intentionally outside schema 1.
