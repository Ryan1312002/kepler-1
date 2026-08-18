"""Development helper: run the verifier image against adversarial result files.

Every case has to leave exactly "0" or "1" in /logs/verifier/reward.txt, and only
a genuinely correct set of concentrations may score 1.
"""

from __future__ import annotations

import json
import shutil
import subprocess
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1] / "hplc-dad-quant"
TRUTH = json.loads((ROOT / "tests" / "ground_truth.json").read_text())["samples"]
IMAGE = "hplc-verifier-test"


def correct_results(scale: float = 1.0) -> dict:
    return {
        "samples": {
            sid: {analyte: value * scale for analyte, value in values.items()}
            for sid, values in TRUTH.items()
        }
    }


CASES: dict[str, str | None] = {
    "missing file": None,
    "empty file": "",
    "not json": "reward=1\n",
    "json list": "[1, 2, 3]\n",
    "json null": "null\n",
    "no samples key": json.dumps({"results": correct_results()["samples"]}),
    "samples is a string": json.dumps({"samples": "S1=18.4"}),
    "empty samples": json.dumps({"samples": {}}),
    "one sample only": json.dumps({"samples": {"S1": correct_results()["samples"]["S1"]}}),
    "extra sample": json.dumps(
        {"samples": {**correct_results()["samples"], "S5": correct_results()["samples"]["S1"]}}
    ),
    "missing analyte": json.dumps(
        {
            "samples": {
                sid: {k: v for k, v in values.items() if k != "riboflavin"}
                for sid, values in correct_results()["samples"].items()
            }
        }
    ),
    "values as strings": json.dumps(
        {
            "samples": {
                sid: {k: str(v) for k, v in values.items()}
                for sid, values in correct_results()["samples"].items()
            }
        }
    ),
    "nan values": '{"samples": {"S1": {"thiamine": NaN}}}',
    "booleans": json.dumps(
        {
            "samples": {
                sid: {k: True for k in values} for sid, values in correct_results()["samples"].items()
            }
        }
    ),
    "zeros": json.dumps(
        {
            "samples": {
                sid: {k: 0 for k in values} for sid, values in correct_results()["samples"].items()
            }
        }
    ),
    "negatives": json.dumps(
        {
            "samples": {
                sid: {k: -v for k, v in values.items()}
                for sid, values in correct_results()["samples"].items()
            }
        }
    ),
    "huge nesting": json.dumps({"samples": {"S1": {"thiamine": {"value": 18.4}}}}),
    "deliberate reward file": json.dumps(correct_results(scale=1.5)),
    "off by 5 percent": json.dumps(correct_results(scale=1.05)),
    "off by 4.5 percent": json.dumps(correct_results(scale=1.045)),
    "off by 3.5 percent": json.dumps(correct_results(scale=1.035)),
    "exact truth": json.dumps(correct_results()),
}

EXPECTED_PASS = {"off by 3.5 percent", "exact truth"}


def build_image() -> None:
    subprocess.run(
        ["docker", "build", "-q", "-t", IMAGE, "-f", "Dockerfile", "."],
        cwd=ROOT / "tests",
        check=True,
        capture_output=True,
    )


def run_case(payload: str | None) -> str:
    with tempfile.TemporaryDirectory() as tmp:
        app = Path(tmp) / "app"
        logs = Path(tmp) / "logs" / "verifier"
        app.mkdir(parents=True)
        logs.mkdir(parents=True)
        if payload is not None:
            (app / "results.json").write_text(payload)
        # The agent gets to write its own reward file first; the verifier must
        # overwrite it, exactly as harbor's fresh bind mount would.
        (logs / "reward.txt").write_text("1")
        subprocess.run(
            [
                "docker", "run", "--rm",
                "-v", f"{app}:/app",
                "-v", f"{Path(tmp) / 'logs'}:/logs",
                IMAGE,
                "bash", "/tests/test.sh",
            ],
            check=False,
            capture_output=True,
        )
        reward_path = logs / "reward.txt"
        if not reward_path.exists():
            return "MISSING"
        return reward_path.read_text().strip()


def main() -> None:
    build_image()
    failures = []
    for label, payload in CASES.items():
        reward = run_case(payload)
        want = "1" if label in EXPECTED_PASS else "0"
        verdict = "ok" if reward == want else "UNEXPECTED"
        if reward != want:
            failures.append(label)
        print(f"{label:26s} reward={reward:8s} want={want}  {verdict}")
    print()
    print("all cases behaved as expected" if not failures else f"PROBLEM CASES: {failures}")
    shutil.rmtree(Path.home() / ".cache" / "nonexistent", ignore_errors=True)


if __name__ == "__main__":
    main()
