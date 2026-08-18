"""Development helper: locate residual bias in the reference chain."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1] / "hplc-dad-quant"
sys.path.insert(0, str(ROOT / "solution"))
sys.path.insert(0, str(ROOT / "tests" / "generator"))

import analyze  # noqa: E402
import generate_dataset as gen  # noqa: E402

DATA = ROOT / "environment" / "data"

data = analyze.Dataset(DATA)
stray = analyze.fit_stray_light(data)
spectra, retention = analyze.pure_spectra(data, stray)
matrix = analyze.matrix_tail_spectrum(data, stray, spectra, retention)
fitted_u = analyze.degradant_spectrum(data, stray, spectra, retention)
true_u = np.array(gen.COMPONENTS["lumichrome"]["spectrum"])


def theoretical_area(component: str, conc: float, v_inj: float) -> float:
    return gen.COMPONENTS[component]["response"] * conc * v_inj


print("=== measured area vs forward-model area (relative error, %) ===")
print("-- calibration mixes --")
for run in data.runs_of("calibration_mix"):
    areas = analyze.measure_run(data, run, stray, spectra, retention)
    levels = run["concentrations_mg_per_L"]
    row = []
    for name in ["thiamine", "nicotinamide", "pyridoxine", "caffeine", "folic_acid", "riboflavin"]:
        # v_inj is unknown per run, so normalise on caffeine to isolate shape bias
        row.append(areas[name])
    ratios = {
        name: areas[name] / areas["caffeine"] / (theoretical_area(name, levels[name], 1.0) / theoretical_area("caffeine", levels["caffeine"], 1.0))
        for name in ["thiamine", "nicotinamide", "pyridoxine", "folic_acid", "riboflavin"]
    }
    print(run["level"], {k: round((v - 1) * 100, 2) for k, v in ratios.items()})

print("-- samples, ratio to caffeine vs forward model --")
for extras_label, extras in [
    ("fitted extras", {"early": {analyze.MATRIX_COMPONENT: matrix}, "pyridoxine": {analyze.MATRIX_COMPONENT: matrix}, "riboflavin": {analyze.DEGRADANT_COMPONENT: fitted_u}}),
    ("true extras", {"early": {analyze.MATRIX_COMPONENT: np.array(gen.COMPONENTS["ascorbate"]["spectrum"])}, "pyridoxine": {analyze.MATRIX_COMPONENT: np.array(gen.COMPONENTS["ascorbate"]["spectrum"])}, "riboflavin": {analyze.DEGRADANT_COMPONENT: true_u}}),
]:
    print(f"  [{extras_label}]")
    for run in data.runs_of("sample"):
        sample_id = run["sample_id"]
        spec = gen.SAMPLES[sample_id]
        areas = analyze.measure_run(data, run, stray, spectra, retention, extras=extras)
        ratios = {}
        for name in ["thiamine", "nicotinamide", "pyridoxine", "folic_acid", "riboflavin"]:
            model = theoretical_area(name, spec["conc"][name], 1.0) / theoretical_area(
                "caffeine", gen.IS_CONC_MG_L, 1.0
            )
            ratios[name] = areas[name] / areas["caffeine"] / model
        print("   ", sample_id, {k: round((v - 1) * 100, 2) for k, v in ratios.items()})

print()
print("degradant spectrum fitted:", np.round(fitted_u, 4).tolist())
print("degradant spectrum true  :", np.round(true_u, 4).tolist())
print("matrix spectrum fitted   :", np.round(matrix, 4).tolist())
print("matrix spectrum true     :", np.round(np.array(gen.COMPONENTS["ascorbate"]["spectrum"]), 4).tolist())
