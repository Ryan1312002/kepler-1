#!/bin/bash
# Verifier entry point. Every path through this script has to leave a reward
# behind, so the reward file is written up front and only raised to 1 if the
# whole test suite passes.
set -u

mkdir -p /logs/verifier
echo 0 > /logs/verifier/reward.txt

python -m pytest /tests/test_results.py \
    -rA \
    -p no:cacheprovider \
    --ctrf /logs/verifier/ctrf.json
status=$?

if [ "${status}" -eq 0 ]; then
    echo 1 > /logs/verifier/reward.txt
else
    echo 0 > /logs/verifier/reward.txt
fi

exit 0
