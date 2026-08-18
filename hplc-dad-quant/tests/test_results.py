"""Verify /app/results.json against the concentrations the solutions were made to.

The grading contract is deliberately blunt: twenty concentrations, each of which
has to land within 4.0 % relative of the value used to prepare the sample. Every
failure mode of an incomplete work-up - leaving the detector non-linearity in,
integrating the split thiamine peak as one lobe, letting the matrix hump or the
riboflavin degradation product into an analyte area, quantifying without the
internal standard, forgetting a dilution factor - lands outside that window on at
least one of the twenty, so there is no way to collect them all except by
resolving the chromatograms properly.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

RESULTS_PATH = Path("/app/results.json")
TRUTH_PATH = Path(__file__).with_name("ground_truth.json")
RELATIVE_TOLERANCE = 0.040

TRUTH = json.loads(TRUTH_PATH.read_text())
ANALYTES = TRUTH["analytes"]
EXPECTED = TRUTH["samples"]
CASES = [(sample_id, analyte) for sample_id in sorted(EXPECTED) for analyte in ANALYTES]


def load_results() -> dict:
    """Parse the agent's output, failing the test rather than erroring out."""
    if not RESULTS_PATH.is_file():
        pytest.fail(f"{RESULTS_PATH} was not written")
    text = RESULTS_PATH.read_text(errors="replace")
    try:
        payload = json.loads(text)
    except json.JSONDecodeError as exc:
        pytest.fail(f"{RESULTS_PATH} is not valid JSON: {exc}")
    if not isinstance(payload, dict):
        pytest.fail(f"{RESULTS_PATH} must hold a JSON object, found {type(payload).__name__}")
    return payload


def load_samples() -> dict:
    payload = load_results()
    samples = payload.get("samples")
    if not isinstance(samples, dict):
        pytest.fail("results.json must contain an object under the key 'samples'")
    return samples


def reported_value(sample_id: str, analyte: str) -> float:
    samples = load_samples()
    if sample_id not in samples:
        pytest.fail(f"no results reported for sample {sample_id}")
    entry = samples[sample_id]
    if not isinstance(entry, dict):
        pytest.fail(f"results for sample {sample_id} must be an object of analyte concentrations")
    if analyte not in entry:
        pytest.fail(f"sample {sample_id} has no value for {analyte}")
    value = entry[analyte]
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        pytest.fail(
            f"{sample_id}/{analyte} must be a JSON number in mg/L, found {value!r}"
        )
    value = float(value)
    if not math.isfinite(value):
        pytest.fail(f"{sample_id}/{analyte} is not a finite number: {value!r}")
    return value


def test_results_file_is_present_and_parses() -> None:
    load_results()


def test_reports_exactly_the_four_samples() -> None:
    samples = load_samples()
    assert set(samples) == set(EXPECTED), (
        f"expected results for {sorted(EXPECTED)}, found {sorted(samples)}"
    )


def test_every_concentration_is_a_positive_number() -> None:
    for sample_id, analyte in CASES:
        value = reported_value(sample_id, analyte)
        assert value > 0.0, f"{sample_id}/{analyte} is not a positive concentration: {value}"


@pytest.mark.parametrize(("sample_id", "analyte"), CASES, ids=[f"{s}-{a}" for s, a in CASES])
def test_concentration_within_tolerance(sample_id: str, analyte: str) -> None:
    expected = float(EXPECTED[sample_id][analyte])
    reported = reported_value(sample_id, analyte)
    deviation = abs(reported - expected) / expected
    assert deviation <= RELATIVE_TOLERANCE, (
        f"{sample_id}/{analyte}: reported {reported:.4g} mg/L, prepared at {expected:.4g} mg/L, "
        f"off by {deviation * 100:.2f} % (limit {RELATIVE_TOLERANCE * 100:.1f} %)"
    )
