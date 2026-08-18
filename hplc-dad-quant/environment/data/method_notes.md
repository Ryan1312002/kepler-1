# AQ-2417 - working notes, vitamin premix beverage

Handover notes for the water soluble vitamin assay. Everything below is what I
wrote down while the sequences were running, so it is the whole story as far as
the instrument is concerned.

## Chromatography

Column is a 150 x 4.6 mm, 5 um C18, held at 40 C. Mobile phase A is 20 mM
potassium phosphate at pH 2.8 with 5 mM hexane-1-sulfonate as the ion pairing
reagent, mobile phase B is methanol. Gradient is 5 % B held for 1 min, up to
65 % B at 9.5 min, held to 11 min, then back down; total cycle 12 min with the
re-equilibration cut off the end of the data file. Flow 1.0 mL/min, injection
10 uL.

The phosphate buffer absorbs in the low UV and the methanol ramp shifts the
refractive index, so the baseline climbs through the run on the short
wavelength channels. That rise is real absorbance in the flow cell, not an
electronic offset. There is a diluent injection at the head of each sequence to
record it.

Detector is the diode array with the 10 mm cell. Channels stored are 246, 254,
261, 273, 290, 324, 350 and 445 nm, all in mAU, no reference wavelength
subtraction, no bunching. Standards sequence was acquired at 5 Hz. The sample
sequence had to be re-acquired after the module firmware update and came off at
4 Hz - annoying, but the data are otherwise the same method.

## Standards

Individual stocks of each vitamin plus caffeine were made up in mobile phase A
and injected on their own, one run each, so we have a clean spectrum and
retention time for every compound we expect. Concentrations are in the
manifest.

The five mixed standards are the actual calibration: all five vitamins together
with caffeine as internal standard at 40.00 mg/L in every level. The autosampler
needle seal is worn and the injection volume wanders by a few percent from
injection to injection, which is exactly why the internal standard is in there.

## Detector linearity

The DAD failed its linearity qualification in July. The service engineer's guess
is stray light - the cell window looks fogged and he wants to swap the flow cell
- but the qualification report is not back yet and the part is on order, so the
sequences were run on the detector as it stands. Do not take the 2.0 AU linear
range from the vendor spec sheet at face value.

The six level linearity check with caffeine on its own was run for exactly this
reason. It is the only characterisation of the detector we have.

## Samples

Four retained bottles of the premix, S1 to S4. Ingredient list on the carton
declares thiamine, niacinamide, pyridoxine, riboflavin and folic acid, plus
ascorbic acid at a much higher level than any of the vitamins, plus sugars and
citrate.

The bottles came back from the distributor and sat in the QC bay for about three
weeks before we got to them. It is warm in there and the bench by the window
gets direct sun in the afternoon. Riboflavin is not happy about that.

Sample prep: degas, filter at 0.45 um, spike caffeine to 40.00 mg/L in the final
solution, i.e. after any dilution. S3 is a concentrate and was diluted 1 in 4
before spiking; dilution factors are in the manifest. The premix would not go
fully into mobile phase, so all four sample solutions were made up in 60:40
methanol/water. The standards were not - they are in mobile phase A.

## Things I noticed and did not get to the bottom of

- The first two minutes of the sample runs look nothing like the standards.
- Peak shape early in the sample sequence is not the same as in the standards
  sequence, and it is the same in all four samples.
- The riboflavin region has more going on in the samples than in the standards.
- Nothing was re-injected. There is no more sample left and the instrument is
  down for the flow cell swap, so what is in the data directory is all there is.
