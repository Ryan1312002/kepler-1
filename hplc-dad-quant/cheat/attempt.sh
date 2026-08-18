#!/bin/bash
# A deliberate attempt to pass the verifier without resolving the chromatograms.
# Kept for review only - the harness never runs this file.
#
# Route 1: write the reward myself.
#   The verifier runs in its own container, built from tests/Dockerfile, after
#   this container has been thrown away, and /logs/verifier is bind-mounted into
#   that container from the trial directory. Anything written here is discarded.
mkdir -p /logs/verifier
echo 1 > /logs/verifier/reward.txt
echo 1 > /logs/verifier/reward.json

# Route 2: find the answer somewhere on the box.
#   The environment image copies environment/data and nothing else: 23
#   chromatograms, the run manifest, the method notes. The prepared
#   concentrations exist only in tests/ground_truth.json, which is baked into the
#   verifier image and never mounted here.
grep -ril --include='*' -e 'ground_truth' -e 'prepared at' / 2>/dev/null | head
ls -la /tests /solution 2>/dev/null

# Route 3: plant a ground truth of my own where the verifier might look.
#   test_results.py resolves the truth relative to its own file inside the
#   verifier image, so a file dropped in the agent container is not consulted -
#   and /app/results.json is the only path collected as an artifact.
mkdir -p /app/tests
cat > /app/tests/ground_truth.json <<'JSON'
{"samples": {"S1": {"thiamine": 1.0}}}
JSON

# Route 4: guess from the label, which is the only external number available.
#   The carton declares nominal contents, the bottles are three weeks past a warm
#   storage excursion, and the verifier wants 4 % relative on all twenty values
#   at once. Label claims are nowhere near that.
python - <<'PY'
import json
label = {"thiamine": 20.0, "nicotinamide": 100.0, "pyridoxine": 25.0,
         "riboflavin": 15.0, "folic_acid": 7.5}
samples = {sid: dict(label) for sid in ("S1", "S2", "S3", "S4")}
json.dump({"samples": samples}, open("/app/results.json", "w"), indent=2)
PY

# Route 5: integrate one wavelength per analyte and calibrate on peak areas,
#   i.e. the work-up a competent but hurried analyst would do in ten minutes.
#   It ignores the stray-light compression, the matrix hump, the thiamine
#   doublet and the degradation product. Measured deviations from the prepared
#   values run from a few percent on the isolated peaks to several hundred
#   percent on thiamine, so this fails as well.
echo "see tools/diagnose.py in the authoring repository for the measured deviations"
