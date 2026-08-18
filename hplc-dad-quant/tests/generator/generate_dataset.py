"""Synthesise the HPLC-DAD data set shipped in environment/data.

This script is the single source of truth for the task. It writes

  * the chromatogram CSV files, the run manifest and the method notes into the
    environment directory (visible to the agent), and
  * ground_truth.json into the tests directory (sealed from the agent).

The forward model, in the order the physical instrument applies it:

  1. every component contributes  area_kj = a_k * s_kj * c_k * v_inj  spread
     over time as a unit-area exponentially modified Gaussian,
  2. the mobile-phase gradient adds its own absorbance background,
  3. the diode array sees the *total* absorbance through a stray-light
     limited photometer:  T_obs = (T_true + f) / (1 + f),
  4. white detector noise is added and the trace is written out in mAU.

Run with  python generate_dataset.py --env-dir <dir> --truth <file>.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.special import erfcx

CHANNELS_NM = [246, 254, 261, 273, 290, 324, 350, 445]

# Peak-normalised UV spectra (unitless, 1.0 at each compound's reference
# channel) followed by the absolute response a_k in AU*min per (mg/L) at that
# reference channel for a nominal 10 uL injection.
COMPONENTS: dict[str, dict] = {
    "thiamine": {
        "spectrum": [1.000, 0.700, 0.880, 0.560, 0.180, 0.030, 0.010, 0.002],
        "response": 0.00280,
        "t_ref": 2.12,
    },
    "nicotinamide": {
        "spectrum": [0.345, 0.700, 1.000, 0.300, 0.048, 0.005, 0.002, 0.001],
        "response": 0.00100,
        "t_ref": 2.38,
    },
    "pyridoxine": {
        "spectrum": [0.590, 1.000, 0.705, 0.245, 0.195, 0.545, 0.095, 0.002],
        "response": 0.00200,
        "t_ref": 3.62,
    },
    "caffeine": {
        "spectrum": [0.300, 0.350, 0.560, 1.000, 0.340, 0.020, 0.005, 0.001],
        "response": 0.00230,
        "t_ref": 5.44,
    },
    "folic_acid": {
        "spectrum": [0.720, 0.620, 0.760, 1.000, 0.815, 0.415, 0.560, 0.012],
        "response": 0.00260,
        "t_ref": 6.86,
    },
    "riboflavin": {
        "spectrum": [0.610, 0.560, 0.740, 1.000, 0.430, 0.290, 0.780, 0.640],
        "response": 0.00340,
        "t_ref": 8.24,
    },
    # Photodegradation product of riboflavin. Present in the stored samples
    # only, never in a freshly prepared standard. Its "concentration" is
    # bookkept in arbitrary units - it is not an analyte and is never reported.
    "lumichrome": {
        "spectrum": [0.640, 0.520, 0.470, 0.560, 0.440, 0.310, 1.000, 0.135],
        "response": 0.00300,
        "t_ref": 8.44,
    },
    # Ascorbic acid from the beverage matrix: several hundred mg/L, barely
    # retained and mass overloaded, so it comes off as a broad hump underneath
    # the first two vitamins. It is what drives the early part of every sample
    # chromatogram deep into the non-linear region of the detector, and it is
    # absent from the standards, which are made up in mobile phase.
    "ascorbate": {
        "spectrum": [1.000, 0.905, 0.640, 0.075, 0.012, 0.003, 0.001, 0.000],
        "response": 0.00168,
        "t_ref": 1.98,
        "sigma": 0.220,
        "tau": 0.140,
    },
}

ANALYTES = ["thiamine", "nicotinamide", "pyridoxine", "riboflavin", "folic_acid"]
INTERNAL_STANDARD = "caffeine"
IS_CONC_MG_L = 40.00

STRAY_LIGHT = 0.0185
NOISE_AU = 3.5e-5
RUN_LENGTH_MIN = 12.0

# Mobile-phase background: a scaled sigmoid ramp per channel.
BG_CHANNEL_SCALE = [1.000, 0.855, 0.720, 0.450, 0.300, 0.165, 0.090, 0.040]
BG_BASE = 0.021
BG_STEP1 = 0.230
BG_STEP2 = 0.105

# Injection-solvent mismatch distorts the least retained analyte in the sample
# sequence into a doublet. Standards are dissolved in mobile phase and stay
# single-peaked.
SPLIT_ANALYTE = "thiamine"
SPLIT_LOBES = [(-0.101, 0.601), (0.158, 0.399)]

SINGLE_STANDARDS = {
    "thiamine": 25.0,
    "nicotinamide": 60.0,
    "pyridoxine": 20.0,
    "caffeine": 40.0,
    "folic_acid": 12.0,
    "riboflavin": 10.0,
}

CAL_MIX_LEVELS = {
    "L1": {"thiamine": 3.20, "nicotinamide": 12.0, "pyridoxine": 3.60, "riboflavin": 2.10, "folic_acid": 1.35},
    "L2": {"thiamine": 8.00, "nicotinamide": 30.0, "pyridoxine": 9.00, "riboflavin": 5.25, "folic_acid": 3.38},
    "L3": {"thiamine": 16.0, "nicotinamide": 62.0, "pyridoxine": 18.0, "riboflavin": 10.5, "folic_acid": 6.75},
    "L4": {"thiamine": 32.0, "nicotinamide": 100.0, "pyridoxine": 30.0, "riboflavin": 18.0, "folic_acid": 9.60},
    "L5": {"thiamine": 44.0, "nicotinamide": 145.0, "pyridoxine": 42.0, "riboflavin": 25.0, "folic_acid": 13.5},
}

LINEARITY_LEVELS = {
    "N1": 5.0,
    "N2": 12.0,
    "N3": 25.0,
    "N4": 45.0,
    "N5": 70.0,
    "N6": 100.0,
}

# Concentration in the solution that was injected. The reported result is this
# value multiplied by the dilution factor.
SAMPLES = {
    "S1": {
        "dilution_factor": 1.0,
        "v_inj": 1.000,
        "conc": {"thiamine": 18.42, "nicotinamide": 96.5, "pyridoxine": 22.74, "riboflavin": 14.18, "folic_acid": 6.85},
        "lumichrome": 3.10,
        "ascorbate": 402.0,
    },
    "S2": {
        "dilution_factor": 1.0,
        "v_inj": 0.9420,
        "conc": {"thiamine": 31.66, "nicotinamide": 128.0, "pyridoxine": 35.91, "riboflavin": 21.57, "folic_acid": 11.28},
        "lumichrome": 5.40,
        "ascorbate": 455.0,
    },
    "S3": {
        "dilution_factor": 4.0,
        "v_inj": 1.0610,
        "conc": {"thiamine": 27.31, "nicotinamide": 116.5, "pyridoxine": 29.63, "riboflavin": 18.92, "folic_acid": 9.41},
        "lumichrome": 6.80,
        "ascorbate": 386.0,
    },
    "S4": {
        "dilution_factor": 1.0,
        "v_inj": 1.0270,
        "conc": {"thiamine": 12.13, "nicotinamide": 91.8, "pyridoxine": 19.44, "riboflavin": 6.72, "folic_acid": 4.23},
        "lumichrome": 9.60,
        "ascorbate": 428.0,
    },
}

CAL_RATE_HZ = 5.0
SAMPLE_RATE_HZ = 4.0


def time_grid(rate_hz: float) -> np.ndarray:
    n = int(round(RUN_LENGTH_MIN * 60.0 * rate_hz)) + 1
    return np.arange(n) / (rate_hz * 60.0)


def peak_shape(
    t: np.ndarray,
    t_r: float,
    sigma: float | None = None,
    tau: float | None = None,
) -> np.ndarray:
    """Unit-area exponentially modified Gaussian centred on t_r."""
    if sigma is None:
        sigma = 0.0260 + 0.00220 * t_r
    if tau is None:
        tau = 0.0160 + 0.00450 * t_r
    u = (t - t_r) / sigma
    k = sigma / tau
    w = (k - u) / np.sqrt(2.0)

    shape = np.zeros_like(t)
    # exp(-u^2/2) * erfcx(w) / (2*tau) is the numerically stable form of the
    # usual EMG expression, except far down the tail where erfcx overflows.
    core = w > -20.0
    shape[core] = np.exp(-0.5 * u[core] ** 2) * erfcx(w[core]) / (2.0 * tau)
    tail = ~core
    if np.any(tail):
        # erfc(w) -> 2 there, leaving the pure exponential tail of the EMG.
        shape[tail] = np.exp(np.clip(0.5 * k**2 - u[tail] * k, -700.0, 700.0)) / tau
    return shape


def background(t: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    ramp = BG_BASE
    ramp = ramp + BG_STEP1 / (1.0 + np.exp(-(t - 4.6) / 1.55))
    ramp = ramp + BG_STEP2 / (1.0 + np.exp(-(t - 9.05) / 0.85))
    bg = ramp[:, None] * np.array(BG_CHANNEL_SCALE)[None, :]
    # Lamp and cell drift make the gradient background reproduce to a few
    # tenths of a percent, not exactly.
    bg = bg * (1.0 + rng.uniform(-0.0025, 0.0025))
    bg = bg + rng.uniform(-0.0010, 0.0010, size=(1, len(CHANNELS_NM)))
    return bg


def apply_stray_light(a_true: np.ndarray) -> np.ndarray:
    transmittance = (10.0 ** (-a_true) + STRAY_LIGHT) / (1.0 + STRAY_LIGHT)
    return -np.log10(transmittance)


def build_run(
    t: np.ndarray,
    contents: dict[str, float],
    *,
    rng: np.random.Generator,
    v_inj: float,
    split: bool,
    drift: float,
) -> np.ndarray:
    absorbance = np.zeros((t.size, len(CHANNELS_NM)))
    for name, conc in contents.items():
        if conc <= 0.0:
            continue
        comp = COMPONENTS[name]
        area = comp["response"] * conc * v_inj
        spectrum = np.array(comp["spectrum"])
        t_r = comp["t_ref"] + drift
        sigma = comp.get("sigma")
        tau = comp.get("tau")
        if split and name == SPLIT_ANALYTE:
            profile = np.zeros_like(t)
            for offset, weight in SPLIT_LOBES:
                profile = profile + weight * peak_shape(t, t_r + offset, sigma, tau)
        else:
            profile = peak_shape(t, t_r, sigma, tau)
        absorbance += area * profile[:, None] * spectrum[None, :]

    total = absorbance + background(t, rng)
    observed = apply_stray_light(total)
    observed = observed + rng.normal(0.0, NOISE_AU, size=observed.shape)
    return observed


def write_csv(path: Path, t: np.ndarray, observed_au: np.ndarray) -> None:
    header = "time_min," + ",".join(f"A{nm}_mAU" for nm in CHANNELS_NM)
    rows = [header]
    milli = observed_au * 1000.0
    for i in range(t.size):
        values = ",".join(f"{v:.3f}" for v in milli[i])
        rows.append(f"{t[i]:.5f},{values}")
    path.write_text("\n".join(rows) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--env-dir", required=True, type=Path)
    parser.add_argument("--truth", required=True, type=Path)
    args = parser.parse_args()

    data_dir = args.env_dir / "data"
    runs_dir = data_dir / "runs"
    runs_dir.mkdir(parents=True, exist_ok=True)
    for stale in runs_dir.glob("*.csv"):
        stale.unlink()

    rng = np.random.default_rng(20240517)
    t_cal = time_grid(CAL_RATE_HZ)
    t_smp = time_grid(SAMPLE_RATE_HZ)

    runs: list[dict] = []

    write_csv(runs_dir / "blank_standards.csv", t_cal, build_run(t_cal, {}, rng=rng, v_inj=1.0, split=False, drift=0.0))
    runs.append(
        {
            "file": "runs/blank_standards.csv",
            "run_type": "blank",
            "sequence": "standards",
            "description": "diluent injection, standards sequence",
        }
    )

    for name, conc in SINGLE_STANDARDS.items():
        drift = float(rng.uniform(-0.035, 0.035))
        observed = build_run(t_cal, {name: conc}, rng=rng, v_inj=1.0, split=False, drift=drift)
        filename = f"single_{name}.csv"
        write_csv(runs_dir / filename, t_cal, observed)
        runs.append(
            {
                "file": f"runs/{filename}",
                "run_type": "single_standard",
                "sequence": "standards",
                "compound": name,
                "concentration_mg_per_L": conc,
            }
        )

    for level, contents in CAL_MIX_LEVELS.items():
        drift = float(rng.uniform(-0.035, 0.035))
        mix = dict(contents)
        mix[INTERNAL_STANDARD] = IS_CONC_MG_L
        observed = build_run(t_cal, mix, rng=rng, v_inj=float(rng.uniform(0.95, 1.05)), split=False, drift=drift)
        filename = f"calmix_{level}.csv"
        write_csv(runs_dir / filename, t_cal, observed)
        runs.append(
            {
                "file": f"runs/{filename}",
                "run_type": "calibration_mix",
                "sequence": "standards",
                "level": level,
                "concentrations_mg_per_L": {**contents, INTERNAL_STANDARD: IS_CONC_MG_L},
            }
        )

    for level, conc in LINEARITY_LEVELS.items():
        drift = float(rng.uniform(-0.035, 0.035))
        observed = build_run(t_cal, {INTERNAL_STANDARD: conc}, rng=rng, v_inj=1.0, split=False, drift=drift)
        filename = f"linearity_{level}.csv"
        write_csv(runs_dir / filename, t_cal, observed)
        runs.append(
            {
                "file": f"runs/{filename}",
                "run_type": "linearity_check",
                "sequence": "standards",
                "compound": INTERNAL_STANDARD,
                "concentration_mg_per_L": conc,
            }
        )

    write_csv(runs_dir / "blank_samples.csv", t_smp, build_run(t_smp, {}, rng=rng, v_inj=1.0, split=False, drift=0.0))
    runs.append(
        {
            "file": "runs/blank_samples.csv",
            "run_type": "blank",
            "sequence": "samples",
            "description": "diluent injection, sample sequence",
        }
    )

    truth: dict[str, dict[str, float]] = {}
    for sample_id, spec in SAMPLES.items():
        contents = dict(spec["conc"])
        contents[INTERNAL_STANDARD] = IS_CONC_MG_L
        contents["lumichrome"] = spec["lumichrome"]
        contents["ascorbate"] = spec["ascorbate"]
        drift = float(rng.uniform(-0.035, 0.035))
        observed = build_run(
            t_smp,
            contents,
            rng=rng,
            v_inj=float(spec["v_inj"]),
            split=True,
            drift=drift,
        )
        filename = f"sample_{sample_id}.csv"
        write_csv(runs_dir / filename, t_smp, observed)
        runs.append(
            {
                "file": f"runs/{filename}",
                "run_type": "sample",
                "sequence": "samples",
                "sample_id": sample_id,
                "dilution_factor": spec["dilution_factor"],
            }
        )
        truth[sample_id] = {
            analyte: round(spec["conc"][analyte] * spec["dilution_factor"], 4) for analyte in ANALYTES
        }

    manifest = {
        "study": "AQ-2417 water soluble vitamin premix, four retained beverage samples",
        "detector": {
            "type": "diode array, 10 mm flow cell",
            "channels_nm": CHANNELS_NM,
            "signal_units": "mAU",
            "acquisition_rate_hz": {"standards": CAL_RATE_HZ, "samples": SAMPLE_RATE_HZ},
            "run_length_min": RUN_LENGTH_MIN,
        },
        "analytes": ANALYTES,
        "internal_standard": {
            "compound": INTERNAL_STANDARD,
            "concentration_mg_per_L": IS_CONC_MG_L,
            "note": (
                "spiked at this concentration into every calibration mix and into every "
                "final sample solution, i.e. after any dilution step"
            ),
        },
        "runs": runs,
    }
    (data_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")

    args.truth.parent.mkdir(parents=True, exist_ok=True)
    args.truth.write_text(
        json.dumps(
            {
                "units": "mg/L in the sample as received, i.e. corrected for the dilution factor",
                "analytes": ANALYTES,
                "samples": truth,
            },
            indent=2,
        )
        + "\n"
    )
    print(f"wrote {len(runs)} runs to {runs_dir}")
    print(f"wrote ground truth to {args.truth}")


if __name__ == "__main__":
    main()
