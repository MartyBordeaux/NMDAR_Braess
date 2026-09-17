# Validation

The pipeline was smoke-tested against the supplied frozen **broad-prior** Step19/18/20/21/23 outputs. This is not the final two-prior server run because the complete Step10 all-samples table was not present in the local artifact set.

Development smoke-test gates:

- historical calibration-score replay: PASS;
- maximum absolute historical-score reconstruction error: approximately `6.8e-11`;
- all 50,000 broad candidates processed;
- Step18 primary strong count reproduced as 1,261/50,000 before recalibration;
- output generation completed for score shift, corrected retention, threshold sensitivity, Step20 local clouds, and Step21 expanded basin.

The broad-only smoke test already indicates that the window correction is **material** rather than cosmetic: the retained broad set changes substantially and the number of strong candidates within the recalibrated 5,000-set broad ensemble increases. These development numbers should not be frozen into the manuscript until the full server run uses the original Step10 table containing both broad and reference priors.
