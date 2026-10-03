#!/usr/bin/env python3
"""Minimal emissivity reuse and direct visibility example (Python 3.11+, NumPy).

python examples/transfer_visibility.py transfer-data --q 2.8 --g-power 4
See docs/transfer-export.md for the finite-camera angular-scale convention.
"""
import argparse
import numpy as np
from verify_transfer import load, intensity, pixels, visibility

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("export")
parser.add_argument("--q", type=float)
parser.add_argument("--g-power", type=float)
parser.add_argument("--fov-uas", type=float, default=125.0)
args = parser.parse_args()
samples, meta = load(args.export)
values = intensity(samples, meta, args.q, args.g_power)
image = pixels(samples, meta, values)
uv = np.array([[0, 0], [2e9, -1e9], [-2e9, 1e9]])  # wavelengths
vis = visibility(samples, meta, values, uv, fov_uas=args.fov_uas, flux_jy=0.6)
np.testing.assert_allclose(vis[0], 0.6, rtol=0, atol=2e-15)
np.testing.assert_allclose(vis[1], vis[2].conjugate(), rtol=0, atol=2e-15)
print(f"Reconstructed image: {image.shape}; unfinished rays: {(samples['finished'] == 0).sum()}")
for baseline, value in zip(uv, vis):
    print(f"uv={baseline}: V={value} Jy")
