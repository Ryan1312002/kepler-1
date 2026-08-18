# afterquery/hplc-dad-quant

Quantify five water soluble vitamins in four beverage samples from a finished
HPLC-DAD acquisition, to 4 % relative on all twenty values at once.

The agent gets 23 chromatograms (eight wavelength channels each), a run manifest
and the analyst's method notes. Getting the numbers out means characterising the
detector's stray-light non-linearity from the linearity series and inverting it on
the raw traces, resolving an early cluster that no single wavelength separates,
recognising that thiamine arrives as a doublet because of the sample diluent,
finding two components that appear in the samples but in none of the standards
(the matrix hump and a riboflavin photodegradation product), calibrating on
internal-standard area ratios, and undoing one 1-in-4 dilution.

## Layout

| path | what it is |
| --- | --- |
| `instruction.md` | the brief handed to the agent |
| `environment/Dockerfile` | python:3.13-slim plus numpy/scipy/pandas/matplotlib |
| `environment/data/` | manifest, method notes and the 23 chromatogram CSVs |
| `solution/analyze.py` | reference work-up; `solution/solve.sh` runs it |
| `tests/test_results.py` | 23 checks: schema plus the twenty concentrations |
| `tests/ground_truth.json` | prepared concentrations, sealed from the agent |
| `tests/generator/generate_dataset.py` | the forward model the data set came from |
| `cheat/attempt.sh` | adversarial attempt, never executed |

## Regenerating the data set

```bash
python tests/generator/generate_dataset.py --env-dir environment --truth tests/ground_truth.json
```

The generator is seeded, so the chromatograms and the ground truth are
reproducible byte for byte.

## Local validation

```bash
harbor run -p . -a oracle -e docker   # scores 1
harbor run -p . -a nop -e docker      # scores 0
```

The authoring repository also carries `tools/diagnose.py`, which runs the
reference chain and nine deliberately incomplete work-ups against the ground
truth. It is what the tolerance in the verifier was set from: the reference is
accurate to 0.46 % worst case, and every shortcut misses by 7 % or more.
