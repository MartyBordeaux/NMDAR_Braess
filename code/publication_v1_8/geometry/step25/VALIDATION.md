# Development validation

This file records code-level validation performed before packaging Step 25. It is not a substitute for the server run on the archived Step-10 file.

The historical Step-10 Latin-hypercube parameter table was deterministically reconstructed from the published Step-10 source ranges, seeds, prior order, and `scipy.stats.qmc.LatinHypercube` construction. The reconstructed broad witness 37318 had inherited calibration forcing approximately `pulse_amplitude_au = 0.2869515067` and `glutamate_tau_ms = 17.49696382`, matching the historical values used downstream.

Checks:

- Step-24 corrected targets reconstructed from all 100,000 stored candidate scores:
  - E/H = 1.0998943768468032
  - L/H = 0.6375504979889217
  - J/H = 257.21994520060724 ms
  - maximum corrected-score replay error = 2.22e-12.
- Historical source-style calibration replay for witness 37318:
  - corrected score = 0.07197988513708545.
- Final Base response replay for witness 37318 at G_peak=0.95, tau_G=10 ms:
  - dt=0.1 ms: r = 1.87451028
  - dt=0.05 ms: r = 1.87446975.
- All 40 Step-24 corrected strong candidates passed the stored-maximum response replay gate.
- In the complete test run, source-style corrected calibration scores for the 40 strong vectors differed from their Step-24 stored values by at most 4.84e-13.
- All 40 remained above the frozen strong threshold at dt=0.05 ms.

The actual server pipeline repeats these gates using `/root/nmda/results_step10.zip`; it does not rely on the reconstructed development table.
