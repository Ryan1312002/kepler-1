"""Development harness.

Reports the accuracy of the reference analysis and of a set of deliberately
incomplete analyses, so the tolerance in the verifier can be set where the
correct chain passes comfortably and every shortcut fails.
"""

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
TRUTH = json.loads((ROOT / "tests" / "ground_truth.json").read_text())["samples"]
TOLERANCE = 4.0


def report(label: str, results: dict[str, dict[str, float]]) -> None:
    worst, where, detail = 0.0, "", []
    for sample_id in sorted(TRUTH):
        for analyte, truth in TRUTH[sample_id].items():
            got = results.get(sample_id, {}).get(analyte)
            if got is None or not np.isfinite(got):
                worst, where = float("inf"), f"{sample_id}/{analyte}"
                detail.append("  n/a")
                continue
            error = (got - truth) / truth * 100.0
            detail.append(f"{error:+.1f}")
            if abs(error) > worst:
                worst, where = abs(error), f"{sample_id}/{analyte}"
    verdict = "PASS" if worst <= TOLERANCE else "fail"
    print(f"{label:36s} {verdict}  worst={worst:8.2f}% at {where:18s} [{' '.join(detail)}]")


def variant(
    *,
    stray: float | None = None,
    subtract_before_linearise: bool = False,
    use_internal_standard: bool = True,
    include_matrix: bool = True,
    resolve_degradant: bool = True,
    apply_dilution: bool = True,
    single_channel: bool = False,
    thiamine_main_lobe_only: bool = False,
    ignore_data_rate: bool = False,
    quadratic_calibration: bool = False,
) -> dict[str, dict[str, float]]:
    """Run the reference chain with one link removed."""
    original_corrected = analyze.corrected_signal
    original_trapezoid = np.trapezoid
    if subtract_before_linearise:
        analyze.corrected_signal = lambda observed, blank, f: analyze.linearise(observed - blank, f)
    if ignore_data_rate:
        # Add up data points instead of integrating over time.
        np.trapezoid = lambda y, x=None, axis=0: np.sum(y, axis=axis)
        analyze.np.trapezoid = np.trapezoid
    try:
        data = analyze.Dataset(DATA)
        f = analyze.fit_stray_light(data) if stray is None else stray
        spectra, retention = analyze.pure_spectra(data, f)

        extras: dict[str, dict[str, np.ndarray]] = {}
        if include_matrix:
            hump = analyze.matrix_tail_spectrum(data, f, spectra, retention)
            extras = {
                "early": {analyze.MATRIX_COMPONENT: hump},
                "pyridoxine": {analyze.MATRIX_COMPONENT: hump},
            }
        shape = analyze.standard_profile(data, f, spectra, retention, "riboflavin")

        def measure(run, is_sample):
            if single_channel:
                return single_channel_areas(data, run, f, retention)
            areas = analyze.measure_run(
                data, run, f, spectra, retention, extras=extras if is_sample else None
            )
            if is_sample and resolve_degradant:
                bound = max(3.0 * areas["riboflavin"], 1e-6)
                areas["riboflavin"], _ = analyze.resolve_riboflavin(
                    data, run, f, spectra, retention, shape, bound
                )
            if is_sample and thiamine_main_lobe_only:
                areas["thiamine"] = main_lobe_area(data, run, f, spectra, retention, extras)
            return areas

        points: dict[str, list[tuple[float, float]]] = {a: [] for a in analyze.ANALYTES}
        for run in data.runs_of("calibration_mix"):
            levels = run["concentrations_mg_per_L"]
            areas = measure(run, is_sample=False)
            for analyte in analyze.ANALYTES:
                if use_internal_standard:
                    points[analyte].append(
                        (
                            levels[analyte] / levels[analyze.INTERNAL_STANDARD],
                            areas[analyte] / areas[analyze.INTERNAL_STANDARD],
                        )
                    )
                else:
                    points[analyte].append((levels[analyte], areas[analyte]))
        factors = {}
        curves = {}
        for analyte, pairs in points.items():
            x = np.array([p[0] for p in pairs])
            y = np.array([p[1] for p in pairs])
            factors[analyte] = float(x @ y / (x @ x))
            if quadratic_calibration:
                # y = b1 x + b2 x^2, forced through the origin.
                design = np.column_stack([x, x**2])
                curves[analyte] = np.linalg.lstsq(design, y, rcond=None)[0]

        is_conc = data.manifest["internal_standard"]["concentration_mg_per_L"]
        results: dict[str, dict[str, float]] = {}
        for run in data.runs_of("sample"):
            areas = measure(run, is_sample=True)
            dilution = run.get("dilution_factor", 1.0) if apply_dilution else 1.0
            values = {}
            for analyte in analyze.ANALYTES:
                if quadratic_calibration:
                    ratio = areas[analyte] / areas[analyze.INTERNAL_STANDARD]
                    b1, b2 = curves[analyte]
                    # Invert y = b1 x + b2 x^2 for the positive root.
                    x = (-b1 + np.sqrt(max(b1 * b1 + 4.0 * b2 * ratio, 0.0))) / (2.0 * b2)
                    values[analyte] = float(x) * is_conc * dilution
                elif use_internal_standard:
                    ratio = areas[analyte] / areas[analyze.INTERNAL_STANDARD]
                    values[analyte] = ratio / factors[analyte] * is_conc * dilution
                else:
                    values[analyte] = areas[analyte] / factors[analyte] * dilution
            results[run["sample_id"]] = values
        return results
    finally:
        analyze.corrected_signal = original_corrected
        np.trapezoid = original_trapezoid
        analyze.np.trapezoid = original_trapezoid


def single_channel_areas(data, run, stray, retention) -> dict[str, float]:
    """One channel per analyte, early cluster split at the valley."""
    channel_for = {
        "thiamine": 0,
        "nicotinamide": 2,
        "pyridoxine": 1,
        "caffeine": 3,
        "folic_acid": 3,
        "riboflavin": 7,
    }
    t, signal, _ = analyze.prepared_run(data, run, stray)
    shift = analyze.retention_shift(t, signal, retention)
    layout = analyze.region_layout(retention)
    areas = {}
    for compound, channel in channel_for.items():
        region = "early" if compound in ("thiamine", "nicotinamide") else compound
        lo, hi, _ = layout[region]
        t_win, matrix = analyze.window(t, signal, lo + shift, hi + shift)
        trace = matrix[:, channel]
        if region == "early":
            middle = 0.5 * (retention["thiamine"] + retention["nicotinamide"]) + shift
            near = (t_win > middle - 0.12) & (t_win < middle + 0.12)
            valley = float(t_win[near][int(np.argmin(trace[near]))])
            part = t_win < valley if compound == "thiamine" else t_win >= valley
            areas[compound] = float(np.trapezoid(trace[part], t_win[part]))
        else:
            areas[compound] = float(np.trapezoid(trace, t_win))
    return areas


def main_lobe_area(data, run, stray, spectra, retention, extras) -> float:
    """Thiamine from its tallest lobe only, as a peak picker would report it."""
    t, signal, _ = analyze.prepared_run(data, run, stray)
    shift = analyze.retention_shift(t, signal, retention)
    lo, hi, knowns = analyze.region_layout(retention)["early"]
    t_win, matrix = analyze.window(t, signal, lo + shift, hi + shift)
    stack = [spectra[k] for k in knowns] + list(extras.get("early", {}).values())
    profiles = analyze.solve_profiles(matrix, np.array(stack))
    trace = profiles[:, 0]
    peak = int(np.argmax(trace))
    left, right = peak, peak
    while left > 0 and trace[left - 1] < trace[left]:
        left -= 1
    while right < trace.size - 1 and trace[right + 1] < trace[right]:
        right += 1
    return float(np.trapezoid(trace[left : right + 1], t_win[left : right + 1]))


def main() -> None:
    outcome = analyze.quantify(DATA)
    print(f"stray light: fitted {outcome['stray_light']:.5f}, true {gen.STRAY_LIGHT:.5f}")
    true_u = np.array(gen.COMPONENTS["lumichrome"]["spectrum"])
    for sample_id, spectrum in outcome["degradant_spectra"].items():
        print(f"  degradant spectrum {sample_id} max deviation {np.abs(spectrum - true_u).max():.4f}")
    true_m = np.array(gen.COMPONENTS["ascorbate"]["spectrum"])
    print(f"  matrix spectrum max deviation {np.abs(outcome['matrix_spectrum'] - true_m).max():.4f}")
    print(f"  tolerance used for the PASS/fail column: +-{TOLERANCE:.1f} %")
    print()

    report("reference chain", outcome["samples"])
    print()
    report("no linearisation", variant(stray=0.0))
    report("blank subtracted before invert", variant(subtract_before_linearise=True))
    report("no internal standard", variant(use_internal_standard=False))
    report("matrix hump ignored", variant(include_matrix=False))
    report("degradant ignored", variant(resolve_degradant=False))
    report("dilution factor ignored", variant(apply_dilution=False))
    report("single wavelength + valley drop", variant(single_channel=True))
    report("thiamine main lobe only", variant(thiamine_main_lobe_only=True))
    report("counts instead of time integral", variant(ignore_data_rate=True))
    report("quadratic curve, no linearisation", variant(stray=0.0, quadratic_calibration=True))
    report("quadratic curve on linearised data", variant(quadratic_calibration=True))


if __name__ == "__main__":
    main()
