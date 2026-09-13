#!/usr/bin/env python3
"""Fast integrity checks for the manuscript-v1.5 headline results.

This script does not rerun the expensive searches. It verifies the deposited
final tables against the numerical claims used in manuscript v1.5.
"""
from pathlib import Path
import json
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
D = ROOT / "data" / "derived" / "final"
checks = []

def check(name, condition, detail=""):
    checks.append((name, bool(condition), str(detail)))

m = pd.read_csv(D / "step15" / "minus40_sweep_separation.csv")
vc = m["signed_separation_class"].value_counts().to_dict()
check("minus40_classification_5_5_1",
      vc.get("potentiating", 0) == 5 and vc.get("suppressing", 0) == 5 and vc.get("unresolved", 0) == 1,
      vc)

j = json.loads((D / "step20" / "01_calibration_score_reconstruction.json").read_text())
check("step20_exact_calibration_PASS", j.get("status") == "PASS", j.get("status"))
check("step20_zero_acceptance_mismatches",
      j.get("historical_acceptance_classification_mismatches") == 0,
      j.get("historical_acceptance_classification_mismatches"))
check("step20_score_r2_one", abs(float(j.get("score_r2_finite_rows", 0.0)) - 1.0) < 1e-12,
      j.get("score_r2_finite_rows"))

b = pd.read_csv(D / "step21" / "01_expanded_basin_summary.csv")
r = b.loc[(b["factor"] - 1.5).abs().idxmin()]
check("step21_factor1p5_joint_fraction",
      abs(float(r["strong_and_compatible_fixed_fraction"]) - 0.6816666666666666) < 1e-6,
      r["strong_and_compatible_fixed_fraction"])

p = pd.read_csv(D / "step22" / "03_primary_all11_posterior_summary.csv")
row = p[(p["parameter"] == "pi_strong_occupancy") & (p["model"] == "M1_two_state_factor_width")].iloc[0]
check("step22_primary_occupancy_about_0p458", abs(float(row["mean"]) - 0.4583333333) < 0.002, row["mean"])

s = pd.read_csv(D / "step23" / "01_threshold_sensitivity_global.csv")
row = s.iloc[(s["threshold"].astype(float) - 1.6359295439505144).abs().argmin()]
check("step23_frozen_threshold_present", abs(float(row["threshold"]) - 1.6359295439505144) < 1e-6, row["threshold"])
check("step23_global_counts_1261_and_1",
      int(row["precal_strong_anywhere_n"]) == 1261 and int(row["accepted_strong_anywhere_n"]) == 1,
      row.to_dict())

failed = [x for x in checks if not x[1]]
for name, passed, detail in checks:
    print(("PASS" if passed else "FAIL"), name, detail)
if failed:
    raise SystemExit(f"{len(failed)} verification checks failed")
print(f"PUBLICATION V1.5 VERIFICATION PASS ({len(checks)} checks)")
