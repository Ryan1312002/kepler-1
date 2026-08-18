# Quantify five vitamins in four beverage samples

You are taking over an unfinished assay. Our lab ran a water soluble vitamin
method on four retained bottles of a premix beverage, together with a full set of
standards, and then the instrument went down for a flow cell swap. The raw
acquisition files are all we have, nobody is going to re-inject anything, and the
release paperwork is waiting on a number for every analyte in every sample.

Everything is under `/app/data`:

- `/app/data/manifest.json` lists all 23 runs, what each one is, the standard
  concentrations, the internal standard level and the dilution factor of each
  sample.
- `/app/data/method_notes.md` is my working notebook for the two sequences:
  chromatography, sample prep, what the detector is doing, and the loose ends I
  never chased down. Read it before you touch the data.
- `/app/data/runs/*.csv` are the chromatograms. First column is retention time
  in minutes, the other eight are the diode array channels in mAU. There is no
  reference wavelength subtracted; what you see is what the detector wrote.

## What I need from you

Report the concentration of thiamine, nicotinamide, pyridoxine, riboflavin and
folic acid in each of the four samples, in mg/L **in the sample as received** -
that is, in the bottle, so anything that was diluted during preparation has to be
put back. Write them to `/app/results.json` in exactly this shape, with JSON
numbers, not strings:

```json
{
  "samples": {
    "S1": {"thiamine": 0.0, "nicotinamide": 0.0, "pyridoxine": 0.0, "riboflavin": 0.0, "folic_acid": 0.0},
    "S2": {"thiamine": 0.0, "nicotinamide": 0.0, "pyridoxine": 0.0, "riboflavin": 0.0, "folic_acid": 0.0},
    "S3": {"thiamine": 0.0, "nicotinamide": 0.0, "pyridoxine": 0.0, "riboflavin": 0.0, "folic_acid": 0.0},
    "S4": {"thiamine": 0.0, "nicotinamide": 0.0, "pyridoxine": 0.0, "riboflavin": 0.0, "folic_acid": 0.0}
  }
}
```

All twenty numbers are checked against the concentrations the four bottles were
actually made up to. A number counts as correct when it is within 4.0 % relative
of that value, and the job is only accepted if all twenty are correct at once. Four percent is a normal release limit for this kind of assay, but it is
tight enough that any effect in these data that you leave uncorrected will push
you outside it, on at least one analyte, and that is enough to fail. There is no
partial credit and no reward for reporting an uncertainty.

## Ground rules

Quantify from the runs in `/app/data/runs` and nothing else. Do not substitute
label claims, typical premix formulations, published concentrations or anything
else from outside for a number you were supposed to measure - the samples are not
a public data set and the answer is not written down anywhere. Reference
absorptivities from the literature will not help you either: the response of this
particular detector on this particular day is what the standards are for.

Anything you install, and anything you write to work the data up, is your own
business. Only `/app/results.json` is read back.

One warning, because it decides whether the numbers come out right. These runs
are not textbook clean. The samples are not the standards, in more than one way,
and the detector is not behaving the way the spec sheet says it should. Work out
what each of those things is doing to the peak areas before you convert an area
into a concentration, rather than after.

You have 14400 seconds to complete this task. Do not cheat by using online solutions or hints specific to this task.
