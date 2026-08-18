#!/bin/bash
# Build the zip to upload on the Frontier Bench submit tab. The archive holds the
# contents of the task directory, not the directory itself.
set -euo pipefail

root="$(cd "$(dirname "$0")/.." && pwd)"
task_dir="${root}/hplc-dad-quant"
out="${root}/hplc-dad-quant.zip"

rm -f "${out}"
cd "${task_dir}"
zip -qr "${out}" . \
    -x 'jobs/*' \
    -x '*/__pycache__/*' \
    -x '*.pyc'

echo "wrote ${out}"
unzip -l "${out}" | tail -n 5
