"""Development helper: inspect the riboflavin window and the degradant fit."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path("hplc-dad-quant/solution").resolve()))
sys.path.insert(0, str(Path("hplc-dad-quant/tests/generator").resolve()))

import analyze  # noqa: E402
import generate_dataset as gen  # noqa: E402

data = analyze.Dataset(Path("hplc-dad-quant/environment/data"))
stray = analyze.fit_stray_light(data)
print("stray fit:", stray, "true:", gen.STRAY_LIGHT)

spectra, retention = analyze.pure_spectra(data, stray)
print("\nretention times:", {k: round(v, 3) for k, v in retention.items()})
for name, s in spectra.items():
    true = np.array(gen.COMPONENTS[name]["spectrum"])
    print(f"{name:14s} fit={np.round(s,4).tolist()}")
    print(f"{'':14s} tru={np.round(true,4).tolist()}  maxdev={np.abs(s-true).max():.4f}")

u = analyze.degradant_spectrum(data, stray, spectra, retention)
true_u = np.array(gen.COMPONENTS["lumichrome"]["spectrum"])
print("\ndegradant fit:", np.round(u, 4).tolist())
print("degradant tru:", np.round(true_u, 4).tolist())

# What the riboflavin window looks like at 445 and 350 nm for the worst sample.
t, signal = analyze.prepared_signal(data, data.runs_of("sample")[3], stray)
mask = (t > 7.8) & (t < 9.8)
print("\nS4 riboflavin window (t, A350, A445):")
for i in np.where(mask)[0][::4]:
    print(f"  {t[i]:.3f} {signal[i,6]:.4f} {signal[i,7]:.4f}")
