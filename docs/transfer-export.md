# Thin-disk transfer export, schema 1

```sh
# Separate build directory and simulation identity; do not install over a research binary.
cargo build --release -p nullgeo-cli --target-dir target/transfer-v1
target/transfer-v1/release/nullgeo render examples/transfer_disk.toml \
  --transfer-export transfer-data
python examples/transfer_visibility.py transfer-data --q 2.8 --g-power 4
python examples/verify_transfer.py --binary target/transfer-v1/release/nullgeo \
  --output scratch/transfer-verification-new
```

Python examples need Python >=3.11 and NumPy. Use an environment containing NumPy.
Export paths are relative to the process working directory, as are ordinary output
paths. The export directory must not exist; its parent must exist. Every render
still requires an ordinary `[[output]]`. The flag is optional; without it there
is no export validation, additional tracing, sample allocation, or disk I/O.
Both the ordinary image and the export use the same `GeometryBuffer` tracing pass.
No existing files in the export destination are overwritten. If an I/O error
interrupts export, it can leave a partial directory; `metadata.toml` is written
last as the completion marker. Writes are not crash-durable transactions.

## Supported scenes

The `stylized` geometrically thin, infinitely opaque equatorial disk, with a
uniform black sky (and black secondary sky, if supplied), on any CLI spacetime
that already supports circular disk orbits: Kerr and Schwarzschild. The effective inner radius remains
`max(requested_r_in, ISCO)`; the outer radius must exceed it. No plunging-region
emission is introduced. Emissivity parameters and geometry must be finite.

Regular N×N grids, deterministic jitter, and adaptive refinement, with or without
jitter, are all supported. `supersample >= 1` and
`supersample_max >= supersample` are required by the camera. There are no other
sampling modes in the current renderer. Missing disks, blackbody/volume disks,
nonblack/image/diagnostic backgrounds, invalid sampling, and spacetimes without
disk orbits fail clearly. A degenerate first crossing that did not become opaque
also fails export rather than pretending it is an opaque first hit.

Changing `emissivity_index` or `g_power` at fixed geometry requires no retracing.
Adaptive selection for this supported disk depends on ray class and escaped
sky direction, not emissivity. Changing camera, metric, radii, numerical
parameters, or sampling requires a new trace.

## Container and metadata

A directory contains:

* `samples.npy`: uncompressed NPY v1.0, packed structured array, shape `(S,)`,
  C order, 103 bytes per record. Load with `np.load(..., allow_pickle=False)`;
  `mmap_mode='r'` permits out-of-core access. No pickle, compression, lossy
  quantization, platform padding, or new Rust dependencies.
* `metadata.toml`: UTF-8 TOML. `schema = "nullgeo.thin-disk-transfer"`,
  `schema_version = 1`; reject unknown versions. `width`, `height`, `sample_count`,
  `record_bytes`, `array` and `source_scene` describe the arrays/files.
* `scene.toml`: verbatim input scene, including the original output specification.

Metadata includes `resolved_scene` with deserialized defaults, resolved sampling
maximum and stylized emission parameters; `effective_r_in`, `effective_r_out`
(the actual post-ISCO annulus); `camera_energy`, `camera_time`; and the complete
trace controls `rtol`, `atol`, `dl_init`, `dl_min`, `dl_max`, `escape_radius`,
`max_steps`. Top-level trace controls are authoritative for resolved integrator
defaults. `resolved_scene.disk.r_in` is the **requested** radius, not the clamp.
The uniform sky's absent intensity means 1; secondary sky absent means the primary
sky. Neither introduces any external texture dependency in this schema.

Provenance fields: `package_version`, `build_id`, `git_revision`,
`git_state_at_build`, `source_digest`, `rustc`, `target`, `profile`, `parallel`.
`build_id` identifies a build invocation using its timestamp, not a content hash.
Git revision/state are best-effort build-time observations (unknown if Git is
unavailable); dirty source is never represented as a clean revision. A revision
plus `dirty` cannot identify the modifications. `source_digest = "unknown"`
is explicit: this format does not claim a source archive or content digest.
The verification driver separately records the executed binary's SHA-256, build
ID and path. Archive that binary, sources and export for a durable simulation
identity; package version alone is insufficient. Rebuilding is a new simulation
identity even though this change does not bump the package version.

## Sample fields

All multibyte fields are explicitly **little endian**; all fields are scalars.
`<f8` = IEEE binary64, `<f4` = IEEE binary32, `<u8` = unsigned 64-bit,
`|u1` = unsigned byte. Masks use 0/1. Columns appear in the following order:

| Field | dtype | Meaning / units |
|---|---|---|
| `pixel_index` | `<u8` | Flattened top-down pixel index: `j * width + i` |
| `sample_index` | `<u8` | Zero-based contributing-sample index within that pixel |
| `offset_x`, `offset_y` | `<f8` each | Actual fractional pixel offsets; right and down, in (0,1) |
| `screen_u`, `screen_v` | `<f8` each | Actual dimensionless pinhole coordinates; right and up |
| `weight` | `<f8` | Mathematical pixel averaging weight, `1 / N_pixel` rounded to f64 |
| `render_weight` | `<f4` | Renderer weight, `1.0f32 / (N_pixel as f32)` |
| `radius` | `<f8` | First valid annulus intersection's spacetime radius, geometric length units |
| `g` | `<f8` | `nu_observer / nu_emitter`, dimensionless, at that crossing |
| `radiance` | `<f4` | Actual stored gray per-ray disk radiance, arbitrary linear units |
| `has_intersection` | `|u1` | `first_crossing` exists; independent of redshift availability |
| `radius_valid` | `|u1` | Intersection radius exists, finite and >0 |
| `g_available` | `|u1` | Crossing's `Option<g>` is present |
| `g_valid` | `|u1` | Available `g` is finite and >0 |
| `radiance_valid` | `|u1` | Stored radiance is finite |
| `outcome` | `|u1` | Raw ray outcome code, below; not the diagnostic ray class |
| `finished` | `|u1` | Outcome is captured, escaped or opaque saturation; not a convergence claim |
| `steps_accepted`, `steps_rejected` | `<u8` each | Integrator step counts |

Outcome codes: **0** captured, **1** escaped primary, **2** escaped secondary,
**3** saturated at the opaque disk, **4** maximum steps (unfinished),
**5** stalled at minimum step (failed). The supported metrics normally only
escape to the primary side; the secondary code is reserved consistently with
`RayOutcome`. Saturation is normal opaque termination, not numerical failure.
Failures and misses remain in the array with their original weights; never
renormalize over hits or successful rays. Missing radius/g values are NaN;
available nonfinite values are retained losslessly with validity mask 0.
`g_available=0` is shaded black by the renderer even if an intersection exists.
The example rejects available-but-invalid g and nonfinite radiance for analysis.
`finished=1` does not guarantee valid emission quantities; inspect masks too.

Geometric units have G=c=1. Radii and camera positions share the metric's units;
for mass=1 they are conventionally M. For Kerr, radius is the spacetime's
spheroidal radius, not Euclidean Cartesian distance. g uses the actual camera
four-velocity and the renderer's existing circular emitter four-velocity.
The camera remains at finite distance, including the example's radius 85M.

## Coordinates, ordering and weights

Let `W=width`, `H=height`, `a=tan(fov_deg*pi/360)`, and `(dx,dy)` be the
exported offsets. The exact camera formula, shared with ray generation, is:

```
u = (2*(i+dx)/W - 1) * a
v = (1 - 2*(j+dy)/H) * a * H/W
local direction = normalize([1, u, v])  # forward, right, up in camera frame
```

Pixel `(0,0)` is top left; +u goes right, +v goes up. PFM stores rows bottom up,
so flip them when comparing with `pixel_index`. Camera pose, up and observer
velocity define the frame; no astronomical east/north interpretation is built in.

Records are pixel-major, then in the precise order used for float32 shading.
Each grid uses `k=cell_y*n+cell_x`. Without jitter offsets are cell centers.
With jitter the fractions within each cell are radical inverses of `k+1` in
bases 2 and 3 (deterministic, identical pattern at every pixel; no RNG seed).
For refined pixels, the n=`supersample_max` grid **replaces** the base grid.
Discarded base probes are not image contributors and are not exported. Infer
refinement from per-pixel sample counts. Every contributing sample is exported,
including black, captured, failed, and unfinished rays.

Each pixel's `weight` sums to 1 up to f64 rounding; all image weights sum to W*H.
These are flat screen/pixel quadrature weights, not solid angles or total-flux
normalization. `render_weight` sums can differ from 1 at float32 precision.

## Reconstructing and changing emission

For each hit with available g, compute in float64:

```
I = g**g_power * (radius/effective_r_in)**(-emissivity_index)
```

Cast I to f32, multiply by `render_weight` in f32, then accumulate sequentially
in sample order in f32. Missing crossings/g contribute zero. All RGB channels
are identical. `radiance` stores the original f32 emission as an independent
check; radius and g retain their original f64 values. No exposure, gamma, tone
mapping or flux normalization applies to the linear PFM output.

NumPy and Rust pow implementations may round differently. The verifier allows
`8*(N_max+2)*eps_f32*peak_image` absolute error, accounting conservatively for
pow/cast and N sequential accumulation roundoffs for the modest, nonnegative
emissivities tested. It reports actual errors, not just pass/fail. Stored
radiance reconstruction is checked exactly. This tolerance is **not** a universal
bound for arbitrary overflowing/underflowing powers. The driver checks six
sampling variants, independently rerenders three changed emissivity settings,
verifies the ordinary PFM is byte-identical with/without export, validates
layout/masks, exercises actual max-step/stalled rays, and rejects unsupported
scenes and invalid sampling. Rust tests additionally cover all outcome codes
and missing/nonfinite field combinations.

## Direct subray visibilities and limitations

The minimal example assigns the full screen width a chosen angular FOV F
(default 125 microarcseconds, converted to radians). Set
`east = u*F/(2*a)`, `north = v*F/(2*a)`; vertical FOV is F*H/W. This is a
**chosen linear image-plane angular scale for a finite-distance camera**. It
neither puts the 85M observer at infinity nor establishes a mass/distance
calibration. Using physical local sky angles instead would require a different
projection and its solid-angle Jacobian (pinhole solid angle is proportional to
`du*dv/(1+u²+v²)^(3/2)`), as well as appropriate observer modeling.

For baselines `(U,V)` in wavelengths, the flat-sky discrete quadrature is

```
mass_s = I_s * weight_s * (F/W)**2
visibility(U,V) = sum_s mass_s * exp(-2*pi*i*(U*east_s + V*north_s))
```

Masses optionally rescale to a specified compact flux (0.6 Jy in the example).
The zero-baseline value equals the summed flux and real brightness guarantees
`V(-U,-V)=conj(V(U,V))`. These are checked numerically. No pixel sinc factor is
used: each sample is a quadrature node, not a uniform pixel. Real EHT baselines
can be supplied directly in wavelength units; batch many baselines to avoid
allocating an enormous baseline-by-sample phase matrix.

Reproducing PFM pixels establishes **reconstruction fidelity**, not physical
convergence. Direct subray integration removes the separate uniform-pixel
approximation but still requires resolution/subsampling, adaptive selection,
geodesic tolerance and intersection-location convergence studies, measured in
the intended visibility/noise norm. Finite camera distance, angular scale,
opacity and emission prescriptions remain modeling assumptions. Unfinished rays
must be addressed, not silently dropped, for scientific accuracy.

A first-hit table cannot vary the inner edge by masking hits: removing the
opaque first hit can reveal a later crossing. Nor can it recover crossings
excluded by the original annulus. A future geometry trace must continue past
opaque hits and record **ordered multiple equatorial crossings** over the full
candidate radial range, with radius, g or emitter-state data, crossing order,
validity, and final termination/completeness information per ray. Adaptive
sampling can also change when disk edges change. Emission inside the ISCO needs
a physical plunging-flow velocity prescription; this export does not enable it.

Measured errors, timing, sizes and the tested binary identity are recorded in
[the validation report](transfer-verification.md).
