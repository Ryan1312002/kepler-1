"""Reference analysis for the AQ-2417 vitamin premix data set.

Outline of the processing chain, in the order it has to happen:

  1. characterise the detector from the six-level linearity series. Peak areas
     only become proportional to concentration once the traces are un-mapped
     through the stray-light transmittance relation T_obs = (T_true + f)/(1 + f),
     so f is fitted by demanding that proportionality;
  2. linearise every trace with that f *before* anything is subtracted - the
     photometric distortion acts on the total absorbance reaching the diode
     array, mobile-phase background included;
  3. subtract the linearised gradient blank of the matching sequence, then take
     the small run-to-run drift out with a straight line under each window;
  4. read the pure spectra of the six knowns off the single-compound standards;
  5. the samples contain two components that no standard contains: the tailing
     matrix peak in front of the vitamins, and a riboflavin degradation product
     buried under riboflavin. The first has a selective window of its own; the
     second is recovered by alternating least squares with the riboflavin
     spectrum held fixed. Both go into the fits as extra components;
  6. solve the absorbance matrix channel-wise for elution profiles in every
     window and integrate them over time - dt differs between the two
     sequences - which also collects both lobes of the split thiamine peak;
  7. calibrate on area ratios against the caffeine internal standard, which
     cancels the injection-volume scatter, then apply each sample's dilution
     factor.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize_scalar

ANALYTES = ["thiamine", "nicotinamide", "pyridoxine", "riboflavin", "folic_acid"]
INTERNAL_STANDARD = "caffeine"

MATRIX_COMPONENT = "matrix_tail"
DEGRADANT_COMPONENT = "degradant"

# Half width of the window used to integrate an isolated standard peak and the
# peak-free margin used for the local baseline.
STANDARD_MARGIN = 0.28
BASELINE_MARGIN = 0.26


class Dataset:
    """The runs listed in the manifest, loaded on demand."""

    def __init__(self, data_dir: Path):
        self.data_dir = data_dir
        self.manifest = json.loads((data_dir / "manifest.json").read_text())
        self.channels = np.array(self.manifest["detector"]["channels_nm"], dtype=float)
        self._cache: dict[str, tuple[np.ndarray, np.ndarray]] = {}

    def run(self, relative_path: str) -> tuple[np.ndarray, np.ndarray]:
        if relative_path not in self._cache:
            raw = np.loadtxt(self.data_dir / relative_path, delimiter=",", skiprows=1)
            self._cache[relative_path] = (raw[:, 0], raw[:, 1:] / 1000.0)
        return self._cache[relative_path]

    def runs_of(self, run_type: str) -> list[dict]:
        return [r for r in self.manifest["runs"] if r["run_type"] == run_type]

    def blank(self, sequence: str) -> tuple[np.ndarray, np.ndarray]:
        for run in self.runs_of("blank"):
            if run["sequence"] == sequence:
                return self.run(run["file"])
        raise KeyError(f"no blank run for sequence {sequence!r}")


def linearise(observed: np.ndarray, stray: float) -> np.ndarray:
    """Invert T_obs = (T_true + f)/(1 + f) for absorbance."""
    transmittance = 10.0 ** (-observed) * (1.0 + stray) - stray
    return -np.log10(np.clip(transmittance, 1e-9, None))


def corrected_signal(observed: np.ndarray, blank: np.ndarray, stray: float) -> np.ndarray:
    """Linearise first, subtract the gradient background second."""
    return linearise(observed, stray) - linearise(blank, stray)


def local_baseline(t: np.ndarray, signal: np.ndarray, lo: float, hi: float) -> np.ndarray:
    """Straight line through the peak-free margins on either side of a window."""
    left = (t >= lo - BASELINE_MARGIN) & (t < lo)
    right = (t > hi) & (t <= hi + BASELINE_MARGIN)
    t_left, t_right = t[left].mean(), t[right].mean()
    a_left = np.median(signal[left], axis=0)
    a_right = np.median(signal[right], axis=0)
    slope = (a_right - a_left) / (t_right - t_left)
    return a_left[None, :] + slope[None, :] * (t - t_left)[:, None]


def window(
    t: np.ndarray, signal: np.ndarray, lo: float, hi: float
) -> tuple[np.ndarray, np.ndarray]:
    """Baseline-corrected absorbance matrix inside one window."""
    corrected = signal - local_baseline(t, signal, lo, hi)
    inside = (t >= lo) & (t <= hi)
    return t[inside], corrected[inside]


def apex(t: np.ndarray, signal: np.ndarray, channel: int, guess: float, span: float) -> float:
    near = (t > guess - span) & (t < guess + span)
    trace = signal[near, channel] - np.median(signal[near, channel])
    return float(t[near][int(np.argmax(trace))])


def solve_profiles(matrix: np.ndarray, spectra: np.ndarray) -> np.ndarray:
    """Least-squares elution profiles, shape (time, component)."""
    solution, *_ = np.linalg.lstsq(spectra.T, matrix.T, rcond=None)
    return solution.T


def second_singular_value(matrix: np.ndarray) -> float:
    return float(np.linalg.svd(matrix, compute_uv=False)[1])


def fit_stray_light(data: Dataset) -> float:
    """Fit f by demanding that linearised areas scale with concentration."""
    runs = data.runs_of("linearity_check")
    concentrations = np.array([r["concentration_mg_per_L"] for r in runs])
    _, blank_raw = data.blank("standards")

    def areas_for(stray: float) -> np.ndarray:
        out = []
        for run in runs:
            t, observed = data.run(run["file"])
            signal = corrected_signal(observed, blank_raw, stray)
            t_r = apex(t, signal, 3, 5.44, 0.40)
            t_win, matrix = window(t, signal, t_r - STANDARD_MARGIN, t_r + STANDARD_MARGIN)
            out.append(np.trapezoid(matrix, t_win, axis=0))
        return np.array(out)

    # Only the channels where the standard really absorbs carry information.
    reference = areas_for(0.0)
    useful = np.where(reference.max(axis=0) > 0.02 * reference.max())[0]

    def objective(stray: float) -> float:
        areas = areas_for(stray)[:, useful]
        total = 0.0
        for column in areas.T:
            slope = float(concentrations @ column / (concentrations @ concentrations))
            model = slope * concentrations
            total += float(np.sum((column - model) ** 2) / np.sum(model**2))
        return total

    result = minimize_scalar(
        objective, bounds=(0.0, 0.05), method="bounded", options={"xatol": 1e-7}
    )
    return float(result.x)


def pure_spectra(data: Dataset, stray: float) -> tuple[dict[str, np.ndarray], dict[str, float]]:
    """Peak-normalised spectrum and retention time of every known compound."""
    _, blank_raw = data.blank("standards")

    spectra: dict[str, np.ndarray] = {}
    retention: dict[str, float] = {}
    for run in data.runs_of("single_standard"):
        t, observed = data.run(run["file"])
        signal = corrected_signal(observed, blank_raw, stray)
        rough = signal - np.median(signal, axis=0)
        strongest = int(np.argmax(rough.max(axis=0)))
        t_r = float(t[int(np.argmax(rough[:, strongest]))])
        t_win, matrix = window(t, signal, t_r - STANDARD_MARGIN, t_r + STANDARD_MARGIN)
        area = np.trapezoid(matrix, t_win, axis=0)
        spectra[run["compound"]] = area / area.max()
        retention[run["compound"]] = t_r
    return spectra, retention


def region_layout(retention: dict[str, float]) -> dict[str, tuple[float, float, list[str]]]:
    """Integration windows and the knowns that elute inside each of them.

    Thiamine and nicotinamide are never resolved from one another, so they share
    a window that also has to be wide enough to hold the matrix peak in front of
    them and both lobes of the split thiamine peak.
    """
    return {
        "early": (
            retention["thiamine"] - 1.10,
            retention["nicotinamide"] + 0.55,
            ["thiamine", "nicotinamide"],
        ),
        "pyridoxine": (
            retention["pyridoxine"] - 0.32,
            retention["pyridoxine"] + 0.42,
            ["pyridoxine"],
        ),
        "caffeine": (
            retention["caffeine"] - 0.32,
            retention["caffeine"] + 0.46,
            ["caffeine"],
        ),
        "folic_acid": (
            retention["folic_acid"] - 0.34,
            retention["folic_acid"] + 0.50,
            ["folic_acid"],
        ),
        "riboflavin": (
            retention["riboflavin"] - 0.36,
            retention["riboflavin"] + 1.35,
            ["riboflavin"],
        ),
    }


def prepared_run(data: Dataset, run: dict, stray: float) -> tuple[np.ndarray, np.ndarray, float]:
    """Linearised, background-corrected trace plus the run's retention shift."""
    t, observed = data.run(run["file"])
    blank_t, blank_raw = data.blank(run["sequence"])
    if blank_raw.shape[0] != observed.shape[0]:
        blank_raw = np.column_stack(
            [np.interp(t, blank_t, blank_raw[:, j]) for j in range(blank_raw.shape[1])]
        )
    return t, corrected_signal(observed, blank_raw, stray), 0.0


def retention_shift(
    t: np.ndarray, signal: np.ndarray, retention: dict[str, float]
) -> float:
    """Run-to-run retention drift, measured on the internal standard peak."""
    return apex(t, signal, 3, retention[INTERNAL_STANDARD], 0.30) - retention[INTERNAL_STANDARD]


def measure_run(
    data: Dataset,
    run: dict,
    stray: float,
    spectra: dict[str, np.ndarray],
    retention: dict[str, float],
    extras: dict[str, dict[str, np.ndarray]] | None = None,
) -> dict[str, float]:
    """Integrated profile area per component for one chromatogram."""
    t, signal, _ = prepared_run(data, run, stray)
    shift = retention_shift(t, signal, retention)
    extras = extras or {}

    areas: dict[str, float] = {}
    for name, (lo, hi, knowns) in region_layout(retention).items():
        components = dict(extras.get(name, {}))
        stack = [spectra[k] for k in knowns] + list(components.values())
        labels = list(knowns) + list(components)
        t_win, matrix = window(t, signal, lo + shift, hi + shift)
        profiles = solve_profiles(matrix, np.array(stack))
        for label, area in zip(labels, np.trapezoid(profiles, t_win, axis=0)):
            areas[label] = float(area)
    return areas


def matrix_tail_spectrum(
    data: Dataset, stray: float, spectra: dict[str, np.ndarray], retention: dict[str, float]
) -> np.ndarray:
    """Spectrum of the broad matrix hump, taken off its selective leading edge.

    Nothing else absorbs between the window start and the front of the thiamine
    doublet, so the shape of the spectrum can be read straight off that stretch.
    """
    lo, hi, _ = region_layout(retention)["early"]
    total = np.zeros(len(data.channels))
    for run in data.runs_of("sample"):
        t, signal, _ = prepared_run(data, run, stray)
        shift = retention_shift(t, signal, retention)
        t_win, matrix = window(t, signal, lo + shift, hi + shift)
        selective = (t_win > retention["thiamine"] + shift - 0.62) & (
            t_win < retention["thiamine"] + shift - 0.31
        )
        total += np.trapezoid(matrix[selective], t_win[selective], axis=0)
    return total / total.max()


def standard_profile(
    data: Dataset, stray: float, spectra: dict[str, np.ndarray], retention: dict[str, float], compound: str
) -> tuple[np.ndarray, np.ndarray]:
    """Unit-area elution profile of a compound, taken from its own standard.

    Returned as a function of time relative to the apex, so it can be shifted
    onto a sample run that has drifted and resampled onto its data rate.
    """
    run = next(r for r in data.runs_of("single_standard") if r["compound"] == compound)
    t, signal, _ = prepared_run(data, run, stray)
    lo, hi, _ = region_layout(retention)[compound]
    t_win, matrix = window(t, signal, lo, hi)
    profile = solve_profiles(matrix, np.array([spectra[compound]]))[:, 0]
    return t_win - retention[compound], profile / np.trapezoid(profile, t_win)


def resolve_riboflavin(
    data: Dataset,
    run: dict,
    stray: float,
    spectra: dict[str, np.ndarray],
    retention: dict[str, float],
    shape: tuple[np.ndarray, np.ndarray],
    upper_bound: float,
) -> tuple[float, np.ndarray]:
    """Riboflavin area in a sample by rank annihilation.

    The degradation product sits on riboflavin's tail, so a two-component least
    squares fit is rotationally ambiguous. What is not ambiguous: riboflavin's
    own peak shape is known from its standard, and once the correct amount of
    riboflavin has been subtracted from the window, what is left is one single
    compound, i.e. a rank-one matrix. Scanning the subtracted amount and
    watching the second singular value of the residual therefore reads the area
    off directly, and the leading singular vector of the residual at the optimum
    is the degradation product's spectrum.
    """
    t, signal, _ = prepared_run(data, run, stray)
    drift = retention_shift(t, signal, retention)
    lo, hi, _ = region_layout(retention)["riboflavin"]
    t_win, matrix = window(t, signal, lo + drift, hi + drift)
    relative = t_win - (retention["riboflavin"] + drift)

    shape_t, shape_p = shape
    best: tuple[float, float, np.ndarray] | None = None
    for delta in np.linspace(-0.03, 0.03, 25):
        profile = np.interp(relative - delta, shape_t, shape_p, left=0.0, right=0.0)
        norm = float(np.trapezoid(profile, t_win))
        if norm <= 0.0:
            continue
        profile = profile / norm
        contribution = np.outer(profile, spectra["riboflavin"])

        def objective(area: float) -> float:
            return second_singular_value(matrix - area * contribution)

        result = minimize_scalar(
            objective, bounds=(0.0, upper_bound), method="bounded", options={"xatol": 1e-8}
        )
        if best is None or result.fun < best[0]:
            residual = matrix - float(result.x) * contribution
            _, _, right = np.linalg.svd(residual, full_matrices=False)
            unknown = right[0] * np.sign(right[0].sum())
            best = (float(result.fun), float(result.x), np.clip(unknown / unknown.max(), 0.0, None))

    if best is None:
        raise RuntimeError("rank annihilation failed to find a solution")
    return best[1], best[2]


def response_factors(
    data: Dataset, stray: float, spectra: dict[str, np.ndarray], retention: dict[str, float]
) -> dict[str, float]:
    """Slope of area ratio against concentration ratio over the mixed standards."""
    points: dict[str, list[tuple[float, float]]] = {a: [] for a in ANALYTES}
    for run in data.runs_of("calibration_mix"):
        levels = run["concentrations_mg_per_L"]
        areas = measure_run(data, run, stray, spectra, retention)
        for analyte in ANALYTES:
            points[analyte].append(
                (
                    levels[analyte] / levels[INTERNAL_STANDARD],
                    areas[analyte] / areas[INTERNAL_STANDARD],
                )
            )

    factors: dict[str, float] = {}
    for analyte, pairs in points.items():
        x = np.array([p[0] for p in pairs])
        y = np.array([p[1] for p in pairs])
        factors[analyte] = float(x @ y / (x @ x))
    return factors


def quantify(data_dir: Path) -> dict:
    data = Dataset(data_dir)
    stray = fit_stray_light(data)
    spectra, retention = pure_spectra(data, stray)
    factors = response_factors(data, stray, spectra, retention)
    matrix_spectrum = matrix_tail_spectrum(data, stray, spectra, retention)
    extras = {
        # The hump reaches into the pyridoxine window as well; carrying it there
        # keeps the local baseline error out of the pyridoxine area.
        "early": {MATRIX_COMPONENT: matrix_spectrum},
        "pyridoxine": {MATRIX_COMPONENT: matrix_spectrum},
    }
    shape = standard_profile(data, stray, spectra, retention, "riboflavin")

    is_conc = data.manifest["internal_standard"]["concentration_mg_per_L"]
    results: dict[str, dict[str, float]] = {}
    degradant: dict[str, np.ndarray] = {}
    for run in data.runs_of("sample"):
        areas = measure_run(data, run, stray, spectra, retention, extras=extras)
        bound = max(3.0 * areas["riboflavin"], 1e-6)
        areas["riboflavin"], degradant[run["sample_id"]] = resolve_riboflavin(
            data, run, stray, spectra, retention, shape, bound
        )
        dilution = run.get("dilution_factor", 1.0)
        results[run["sample_id"]] = {
            analyte: round(
                areas[analyte] / areas[INTERNAL_STANDARD] / factors[analyte] * is_conc * dilution,
                4,
            )
            for analyte in ANALYTES
        }
    return {
        "stray_light": stray,
        "samples": results,
        "matrix_spectrum": matrix_spectrum,
        "degradant_spectra": degradant,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", type=Path, default=Path("/app/data"))
    parser.add_argument("--out", type=Path, default=Path("/app/results.json"))
    args = parser.parse_args()

    outcome = quantify(args.data_dir)
    args.out.write_text(json.dumps({"samples": outcome["samples"]}, indent=2) + "\n")
    print(f"fitted stray-light fraction: {outcome['stray_light']:.5f}")
    for sample_id, values in outcome["samples"].items():
        print(f"{sample_id}: " + ", ".join(f"{k}={v:.2f}" for k, v in values.items()))


if __name__ == "__main__":
    main()
