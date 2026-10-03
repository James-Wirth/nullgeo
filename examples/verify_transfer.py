#!/usr/bin/env python3
"""Verify schema v1 with fresh ordinary renders; Python >=3.11 and NumPy only.

Example: python examples/verify_transfer.py --binary target/transfer-v1/release/nullgeo
All generated scenes, renders and exports go in a NEW directory. No installs or
writes to kerr-sbi. Import load(), intensity(), pixels(), visibility() for reuse.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile
import time
import tomllib

import numpy as np


def load(directory):
    directory = Path(directory)
    meta = tomllib.loads((directory / "metadata.toml").read_text())
    assert meta["schema"] == "nullgeo.thin-disk-transfer" and meta["schema_version"] == 1
    samples = np.load(directory / meta["array"], allow_pickle=False, mmap_mode="r")
    assert samples.shape == (meta["sample_count"],)
    assert samples.dtype.itemsize == meta["record_bytes"] == 103
    return samples, meta


def intensity(samples, meta, q=None, g_power=None):
    disk = meta["resolved_scene"]["disk"]
    q = disk["emissivity_index"] if q is None else q
    g_power = disk["g_power"] if g_power is None else g_power
    # Missing redshift is black in the current renderer. Available but invalid
    # physics/numerics must be investigated, not silently treated as a disk miss.
    assert np.all(samples["g_valid"] == samples["g_available"])
    assert np.all(samples["radiance_valid"])
    valid = samples["radius_valid"].astype(bool) & samples["g_valid"].astype(bool)
    result = np.zeros(len(samples), dtype=np.float64)
    result[valid] = samples["g"][valid] ** g_power * (
        samples["radius"][valid] / meta["effective_r_in"]
    ) ** -q
    return result


def pixels(samples, meta, values):
    # Rust casts each emitted value to f32, multiplies by an f32 weight and
    # accumulates sequentially in contributing-sample order. np.add.at preserves it.
    result = np.zeros(meta["width"] * meta["height"], dtype=np.float32)
    weighted = values.astype(np.float32) * samples["render_weight"]
    np.add.at(result, samples["pixel_index"].astype(np.intp), weighted)
    return result.reshape(meta["height"], meta["width"])


def read_pfm(path):
    with open(path, "rb") as f:
        assert f.readline() == b"PF\n"
        w, h = map(int, f.readline().split())
        assert float(f.readline()) == -1.0
        return np.frombuffer(f.read(), dtype="<f4").reshape(h, w, 3)[::-1]


def visibility(samples, meta, values, uv, fov_uas=125.0, flux_jy=0.6):
    """uv in wavelengths; assigned east=screen right, north=screen up.

    This assigns a flat angular image scale to a finite-distance camera. It is
    NOT an observer-at-infinity mapping or an exact solid-angle quadrature.
    """
    fov_rad = fov_uas * np.pi / (180 * 3600 * 1e6)
    screen_width = 2 * np.tan(np.deg2rad(meta["resolved_scene"]["camera"]["fov_deg"]) / 2)
    east = samples["screen_u"] * fov_rad / screen_width
    north = samples["screen_v"] * fov_rad / screen_width
    pixel_area = (fov_rad / meta["width"]) ** 2
    masses = values * samples["weight"] * pixel_area
    if flux_jy is not None:
        masses *= flux_jy / masses.sum()
    uv = np.asarray(uv, dtype=np.float64).reshape(-1, 2)
    # For real EHT datasets batch baselines to bound temporary array memory.
    return np.exp(-2j * np.pi * (uv[:, :1] * east + uv[:, 1:] * north)) @ masses


def radical_inverse(k, base):
    value, denom = 0.0, 1.0
    while k:
        denom *= base
        value += (k % base) / denom
        k //= base
    return value


def check_layout(samples, meta):
    width, height = meta["width"], meta["height"]
    camera = meta["resolved_scene"]["camera"]
    ids = samples["pixel_index"].astype(np.intp)
    counts = np.bincount(ids, minlength=width * height)
    assert len(counts) == width * height
    assert np.all(np.isin(counts, [camera["supersample"] ** 2, camera["supersample_max"] ** 2]))
    np.testing.assert_array_equal(ids, np.repeat(np.arange(width * height), counts))
    np.testing.assert_allclose(np.bincount(ids, weights=samples["weight"]), 1, rtol=0, atol=2e-15)
    np.testing.assert_array_equal(samples["render_weight"], np.float32(1) / counts[ids].astype(np.float32))
    start = 0
    for count in counts:
        chunk = samples[start:start + count]
        np.testing.assert_array_equal(chunk["sample_index"], np.arange(count))
        n = int(np.sqrt(count))
        for k in range(count):
            dx = radical_inverse(k + 1, 2) if camera["jitter"] else 0.5
            dy = radical_inverse(k + 1, 3) if camera["jitter"] else 0.5
            np.testing.assert_allclose([chunk["offset_x"][k], chunk["offset_y"][k]],
                                       [(k % n + dx) / n, (k // n + dy) / n], rtol=0, atol=2e-16)
        start += count
    scale = np.tan(np.deg2rad(camera["fov_deg"]) / 2)
    np.testing.assert_allclose(samples["screen_u"],
        (2 * ((ids % width + samples["offset_x"]) / width) - 1) * scale, rtol=1e-14, atol=1e-16)
    np.testing.assert_allclose(samples["screen_v"],
        (1 - 2 * ((ids // width + samples["offset_y"]) / height)) * scale * height / width,
        rtol=1e-14, atol=1e-16)
    assert np.all(samples["screen_u"][ids % width == 0] < 0)
    assert np.all(samples["screen_v"][ids // width == 0] > 0)
    hit = samples["has_intersection"].astype(bool)
    np.testing.assert_array_equal(samples["radius_valid"], hit)
    assert np.all(np.isnan(samples["radius"][~hit]))
    assert np.all(np.isnan(samples["g"][~samples["g_available"].astype(bool)]))
    assert np.all(samples["g_available"] <= samples["has_intersection"])
    assert np.all(np.isin(samples["outcome"], range(6)))
    np.testing.assert_array_equal(samples["finished"], samples["outcome"] < 4)
    np.testing.assert_array_equal(samples["outcome"] == 3, hit)
    assert np.all(samples["radius"][hit] >= meta["effective_r_in"])
    assert np.all(samples["radius"][hit] <= meta["effective_r_out"])
    return counts


def compare(samples, meta, values, pfm):
    actual = read_pfm(pfm)
    reconstructed = pixels(samples, meta, values)
    np.testing.assert_array_equal(actual[:, :, 0], actual[:, :, 1])
    np.testing.assert_array_equal(actual[:, :, 0], actual[:, :, 2])
    expected = actual[:, :, 0]
    # A conservative f32 rounding bound: O(N*eps) accumulation plus pow/cast.
    # This is a reconstruction tolerance, not a geodesic/visibility error budget.
    n = max(np.bincount(samples["pixel_index"].astype(np.intp)))
    tolerance = 8 * (n + 2) * np.finfo(np.float32).eps * max(float(np.max(np.abs(expected))), 1e-30)
    error = reconstructed.astype(float) - expected
    assert np.max(np.abs(error)) <= tolerance, (np.max(np.abs(error)), tolerance)
    return {"max_abs": float(np.max(np.abs(error))),
            "relative_l2": float(np.linalg.norm(error) / max(np.linalg.norm(expected), 1e-30)),
            "tolerance_abs": tolerance}


def scene_text(output, n, nmax, jitter, q=2.0, power=3.0, max_steps=20000, tol=1e-8):
    return f'''[metric]
kind = "kerr"
mass = 1.0
spin = 0.7
[camera]
position = [-60.10407640085654, 0.0, 60.10407640085654]
fov_deg = 24.0
width = 32
height = 24
supersample = {n}
supersample_max = {nmax}
jitter = {str(jitter).lower()}
[disk]
model = "stylized"
r_in = 8.0
r_out = 18.0
emissivity_index = {q}
g_power = {power}
[sky]
uniform = [0.0, 0.0, 0.0]
[integrator]
tol = {tol}
max_steps = {max_steps}
[[output]]
path = {json.dumps(str(output))}
'''


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--output", type=Path, help="new directory; otherwise use a fresh temporary directory")
    args = parser.parse_args()
    binary = args.binary.resolve()
    if args.output:
        root = args.output.resolve()
        root.mkdir(parents=True, exist_ok=False)
    else:
        root = Path(tempfile.mkdtemp(prefix="nullgeo-transfer-"))
    report = {"binary": str(binary), "binary_sha256": hashlib.sha256(binary.read_bytes()).hexdigest(),
              "run_directory": str(root), "modes": {}}

    def render(config, export=None, success=True):
        command = [str(binary), "render", str(config)]
        if export is not None:
            command += ["--transfer-export", str(export)]
        start = time.perf_counter()
        result = subprocess.run(command, capture_output=True, text=True)
        elapsed = time.perf_counter() - start
        if success:
            assert result.returncode == 0, result.stderr
        else:
            assert result.returncode != 0, command
        return elapsed

    modes = [("single", 1, 1, False), ("single_jitter", 1, 1, True),
             ("grid", 3, 3, False), ("jitter", 3, 3, True),
             ("adaptive", 1, 3, False), ("adaptive_jitter", 2, 3, True)]
    for name, n, nmax, jitter in modes:
        config, pfm, export = root / f"{name}.toml", root / f"{name}.pfm", root / name
        config.write_text(scene_text(pfm, n, nmax, jitter))
        ordinary_time = render(config)
        ordinary_bytes = pfm.read_bytes()
        export_time = render(config, export)
        assert ordinary_bytes == pfm.read_bytes(), "export changed ordinary render"
        samples, meta = load(export)
        counts = check_layout(samples, meta)
        assert np.any(samples["has_intersection"]), "test scene must hit the disk"
        if nmax > n:
            assert len(set(counts)) == 2, "exercise both base and replacement samples"
        values = intensity(samples, meta)
        # The stored reference radiance checks byte order and f32 accumulation exactly.
        np.testing.assert_array_equal(pixels(samples, meta, samples["radiance"]), read_pfm(pfm)[:, :, 0])
        entry = {"samples": len(samples), "sample_counts": sorted(map(int, set(counts))),
                 "original": compare(samples, meta, values, pfm), "variants": [],
                 "ordinary_seconds": ordinary_time, "export_seconds": export_time,
                 "export_bytes": sum(p.stat().st_size for p in export.iterdir()),
                 "pfm_bytes": pfm.stat().st_size,
                 "outcomes": {str(code): int(np.sum(samples["outcome"] == code)) for code in range(6)}}
        for q, power in [(2.8, 3.0), (2.0, 4.2), (3.1, 1.7)]:
            variant_config = root / f"{name}-q{q}-g{power}.toml"
            variant_pfm = root / f"{name}-q{q}-g{power}.pfm"
            variant_config.write_text(scene_text(variant_pfm, n, nmax, jitter, q, power))
            render(variant_config)
            entry["variants"].append({"q": q, "g_power": power,
                **compare(samples, meta, intensity(samples, meta, q, power), variant_pfm)})
        vis = visibility(samples, meta, values, [[0, 0], [2e9, -1e9], [-2e9, 1e9]])
        np.testing.assert_allclose(vis[0], 0.6, rtol=0, atol=2e-15)
        np.testing.assert_allclose(vis[1], vis[2].conjugate(), rtol=0, atol=2e-15)
        entry["visibility_zero_jy"] = [vis[0].real, vis[0].imag]
        entry["visibility_conjugacy_error_jy"] = float(abs(vis[1] - vis[2].conjugate()))
        report["modes"][name] = entry
        report["build_id"] = meta["build_id"]

    for name, steps, tol, expected in [("unfinished", 0, 1e-8, 4), ("stalled", 20000, 1e-300, 5)]:
        config, pfm, export = root / f"{name}.toml", root / f"{name}.pfm", root / name
        config.write_text(scene_text(pfm, 1, 1, False, max_steps=steps, tol=tol))
        render(config, export)
        samples, meta = load(export)
        check_layout(samples, meta)
        assert np.all(samples["outcome"] == expected)
        assert not np.any(samples["finished"])
        assert not np.any(samples["has_intersection"])
        compare(samples, meta, intensity(samples, meta), pfm)

    config = root / "rejected.toml"
    valid = scene_text(root / "rejected.pfm", 1, 1, False)
    for i, text in enumerate([
        valid.replace('model = "stylized"', 'model = "blackbody"').replace('emissivity_index = 2.0\n', '').replace('g_power = 3.0\n', ''),
        valid.replace('uniform = [0.0, 0.0, 0.0]', 'uniform = [1.0, 0.0, 0.0]'),
        valid.replace('supersample = 1', 'supersample = 0'),
        valid.replace('supersample_max = 1', 'supersample_max = 0'),
        valid.replace('[disk]\nmodel = "stylized"\nr_in = 8.0\nr_out = 18.0\nemissivity_index = 2.0\ng_power = 3.0\n', ''),
        valid.replace('kind = "kerr"', 'kind = "reissner-nordstrom"'),
        valid.replace('mass = 1.0', 'mass = nan'),
    ]):
        config.write_text(text)
        export = root / f"rejected-{i}"
        render(config, export, success=False)
        assert not export.exists()

    # Alternating repeated trials reduce warm-up and timing-order effects. Export
    # timings include validation, serialization, and filesystem writes, not fsync.
    config, pfm = root / "benchmark.toml", root / "benchmark.pfm"
    config.write_text(scene_text(pfm, 3, 3, False))
    plain, exported = [], []
    for i in range(7):
        if i % 2:
            exported.append(render(config, root / f"benchmark-{i}"))
            plain.append(render(config))
        else:
            plain.append(render(config))
            exported.append(render(config, root / f"benchmark-{i}"))
    report["benchmark"] = {"ordinary_median_seconds": float(np.median(plain)),
        "export_median_seconds": float(np.median(exported)),
        "overhead_percent": float((np.median(exported) / np.median(plain) - 1) * 100),
        "ordinary_seconds": plain, "export_seconds": exported}
    (root / "verification.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
