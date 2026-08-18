#!/bin/bash
set -euo pipefail

# The reference work-up. Reads the acquisition data under /app/data and writes
# the twenty concentrations to /app/results.json.
python /solution/analyze.py --data-dir /app/data --out /app/results.json

cat /app/results.json
