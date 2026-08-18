"""Development helper: report absorbance levels and peak structure."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

DATA = Path(sys.argv[1] if len(sys.argv) > 1 else "hplc-dad-quant/environment/data")


def load(name):
    raw = np.loadtxt(DATA / name, delimiter=",", skiprows=1)
    return raw[:, 0], raw[:, 1:] / 1000.0


manifest = json.loads((DATA / "manifest.json").read_text())
channels = manifest["detector"]["channels_nm"]

print("channels:", channels)
for run in manifest["runs"]:
    t, a = load(run["file"])
    peak_idx = int(np.argmax(a.max(axis=1)))
    print(
        f"{run['file']:32s} n={t.size:5d} dt={np.median(np.diff(t))*60:.3f}s "
        f"maxA={a.max():.3f} at t={t[peak_idx]:.2f} "
        f"per-channel max={np.round(a.max(axis=0), 3).tolist()}"
    )

print()
print("early cluster in sample S1 (1.7-2.9 min), channel 246 and 261:")
t, a = load("runs/sample_S1.csv")
mask = (t > 1.7) & (t < 2.9)
sub_t = t[mask]
c246 = a[mask, 0]
c261 = a[mask, 2]
for i in range(0, sub_t.size, 6):
    print(f"  {sub_t[i]:.3f}  {c246[i]:.4f}  {c261[i]:.4f}")
