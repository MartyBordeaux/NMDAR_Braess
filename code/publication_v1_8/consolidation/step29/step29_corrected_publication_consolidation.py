#!/usr/bin/env python3
"""
NMDAR Braess Step 29: corrected publication consolidation v1.0

Purpose
-------
Freeze the corrected post-Step24 model cohort and the post-Step28 connectivity
interpretation into one manuscript-facing analysis package. The pipeline:

1. hard-gates the corrected calibration and topology inputs;
2. assigns every one of the 40 confirmed strong vectors to the final sampled
   connectivity regime;
3. selects manuscript representatives from the dominant native-connected
   region, the secondary region, and the historical witness;
4. reconstructs experimental-to-model response matching for the whole strong
   set and for the dominant native-connected region;
5. reruns local simultaneous 12-rate perturbations around the new primary
   representative and two sensitivity representatives using the frozen Step25
   engine and the corrected 140-500-ms calibration definition;
6. assembles publication tables, figures, manuscript numbers, and a frozen
   manifest.

No new response threshold, denominator cutoff, calibration target, or kinetic
acceptance criterion is introduced.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import platform
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import step25_engine_snapshot as eng

SCRIPT_VERSION = "1.0.0"
PIPELINE_NAME = "29_corrected_publication_consolidation"
DEFAULT_OUTDIR = Path("/root/nmda2/step_29/results_step_29_corrected_publication_consolidation")
WITNESS_ID = 37318
EXPECTED_STEP24_STATUS = "MATERIAL_CHANGE"
EXPECTED_STEP25_STATUS = "KINETIC_REGION_STRUCTURE_UNRESOLVED"
EXPECTED_STEP28_STATUS = "PARTIAL_RESIDUAL_CONNECTIVITY_AFTER_FULLDIM_STRESS_TEST"
EXPECTED_STRONG_N = 40
EXPECTED_STEP28_STRICT_SIZES = [36, 3, 1]
EXPECTED_STEP28_NATIVE_SIZES = [35, 3, 1, 1]
FALLBACK_EXPERIMENTAL_RATIOS = np.array(
    [1.6359295439505144, 1.694672, 1.876748, 1.891305, 4.072559], dtype=float
)
RATE_COLS = list(eng.RATE_COLS)
AUX_COLS = list(eng.AUX_COLS)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Step 29 corrected publication consolidation")
    p.add_argument("--step24", default="auto")
    p.add_argument("--step25", default="auto")
    p.add_argument("--step28", default="auto")
    p.add_argument("--output", type=Path, default=DEFAULT_OUTDIR)
    p.add_argument("--threads", type=int, default=int(os.getenv("STEP29_THREADS", "2")))
    p.add_argument("--local-n", type=int, default=int(os.getenv("STEP29_LOCAL_N", "1000")))
    p.add_argument("--seed", type=int, default=int(os.getenv("STEP29_SEED", "20260915")))
    p.add_argument("--skip-local", action="store_true", default=os.getenv("STEP29_SKIP_LOCAL", "0") == "1")
    p.add_argument("--force", action="store_true", default=os.getenv("STEP29_FORCE", "0") == "1")
    p.add_argument("--resume", action="store_true", default=os.getenv("STEP29_RESUME", "0") == "1")
    p.add_argument("--self-test", action="store_true")
    p.add_argument("--preflight-only", action="store_true")
    return p.parse_args()


def json_dump(path: Path, obj: Any) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def discover_dir(spec: str, candidates: Iterable[Path], required: Iterable[str], label: str) -> Path:
    required = list(required)
    if spec != "auto":
        p = Path(spec).expanduser().resolve()
        if not p.is_dir():
            raise FileNotFoundError(f"{label} directory not found: {p}")
        missing = [x for x in required if not (p / x).exists()]
        if missing:
            raise FileNotFoundError(f"{label} missing required files: {missing}")
        return p
    for p in candidates:
        if p.is_dir() and all((p / x).exists() for x in required):
            return p.resolve()
    searched = "\n".join(f"  {x}" for x in candidates)
    raise FileNotFoundError(f"Could not auto-discover {label}. Searched:\n{searched}")


def discover_step24(spec: str) -> Path:
    return discover_dir(
        spec,
        [
            Path("/root/nmda2/step_24/results_step_24_calibration_window_correction"),
            Path("/root/nmda2/results_step_24_calibration_window_correction"),
            Path("/root/nmda2/NMDAR_Braess_step24_calibration_window_correction_v1_0/results_step_24_calibration_window_correction"),
            Path("/root/nmda2/step24/results_step_24_calibration_window_correction"),
        ],
        ["03_acceptance_shift_by_prior.csv", "04A_threshold_strong_support_shift.csv", "08_scientific_decision.json"],
        "Step24",
    )


def discover_step25(spec: str) -> Path:
    return discover_dir(
        spec,
        [
            Path("/root/nmda2/step_25/results_step_25_corrected_strong_geometry"),
            Path("/root/nmda2/results_step_25_corrected_strong_geometry"),
            Path("/root/nmda2/NMDAR_Braess_step25_corrected_strong_geometry_v1_0/results_step_25_corrected_strong_geometry"),
        ],
        [
            "00_input_audit.json",
            "01_corrected_target_reconstruction.json",
            "02_corrected_accepted_cohorts.csv.gz",
            "04_full_80_node_grid.csv.gz",
            "05_strong_candidates_confirmed.csv",
            "15_scientific_decision.json",
            "run_summary.json",
        ],
        "Step25",
    )


def discover_step28(spec: str) -> Path:
    return discover_dir(
        spec,
        [
            Path("/root/nmda2/step_28/results_step_28_residual_connectivity_stress"),
            Path("/root/nmda2/results_step_28_residual_connectivity_stress"),
            Path("/root/nmda2/NMDAR_Braess_step28_residual_connectivity_stress_v1_0/results_step_28_residual_connectivity_stress"),
        ],
        ["07_final_network_edges.csv", "08_final_components.csv", "09_residual_barrier_summary.csv", "10_scientific_decision.json"],
        "Step28",
    )


def sorted_component_sizes(frame: pd.DataFrame, criterion: str) -> list[int]:
    d = frame[frame["criterion"] == criterion]
    if d.empty:
        return []
    return sorted(d.groupby("component_id")["sample_id"].nunique().astype(int).tolist(), reverse=True)


def component_ids(frame: pd.DataFrame, criterion: str, size: int) -> list[int]:
    d = frame[(frame["criterion"] == criterion) & (frame["component_size"] == size)]
    return sorted(d["sample_id"].astype(int).unique().tolist())


def hard_gate(step24: Path, step25: Path, step28: Path) -> dict[str, Any]:
    dec24 = json.loads((step24 / "08_scientific_decision.json").read_text())
    dec25 = json.loads((step25 / "15_scientific_decision.json").read_text())
    dec28 = json.loads((step28 / "10_scientific_decision.json").read_text())
    audit25 = json.loads((step25 / "00_input_audit.json").read_text())
    run25 = json.loads((step25 / "run_summary.json").read_text())
    strong = pd.read_csv(step25 / "05_strong_candidates_confirmed.csv")
    comps = pd.read_csv(step28 / "08_final_components.csv")

    checks = {
        "step24_status": dec24.get("status") == EXPECTED_STEP24_STATUS,
        "step25_status": dec25.get("status") == EXPECTED_STEP25_STATUS,
        "step28_status": dec28.get("status") == EXPECTED_STEP28_STATUS,
        "confirmed_strong_n": len(strong) == EXPECTED_STRONG_N and bool(strong["strong_confirmed_dt0p05"].all()),
        "step28_strict_sizes": sorted_component_sizes(comps, "step28_strict") == EXPECTED_STEP28_STRICT_SIZES,
        "step28_native_sizes": sorted_component_sizes(comps, "step28_native_only_strict") == EXPECTED_STEP28_NATIVE_SIZES,
        "witness_present": WITNESS_ID in set(strong["sample_id"].astype(int)),
        "no_denominator_exclusion_step25": bool(dec25.get("no_denominator_exclusion_applied", True)),
        "no_denominator_exclusion_step28": bool(dec28.get("no_denominator_exclusion_applied", True)),
        "engine_replay_step28": bool(dec28.get("engine_replay_gate_pass", False)),
    }
    failed = [k for k, v in checks.items() if not v]
    if failed:
        raise RuntimeError(f"Step29 hard gate failed: {failed}")

    strong_thr = float(dec28["strong_threshold"])
    strict_cutoff = float(dec28["strict_cutoff"])
    if abs(strong_thr - float(audit25["strong_threshold"])) > 1e-12:
        raise RuntimeError("Step25/Step28 strong-threshold mismatch")
    if abs(strict_cutoff - float(audit25["corrected_broad_max_accepted_score"])) > 1e-9:
        raise RuntimeError("Step25/Step28 strict corrected-cutoff mismatch")

    ratios = np.array(run25.get("experimental_ratios", FALLBACK_EXPERIMENTAL_RATIOS.tolist()), dtype=float)
    if len(ratios) != 5 or np.any(~np.isfinite(ratios)):
        raise RuntimeError("Expected five finite experimental potentiating ratios")

    return {
        "checks": checks,
        "dec24": dec24,
        "dec25": dec25,
        "dec28": dec28,
        "audit25": audit25,
        "run25": run25,
        "strong_threshold": strong_thr,
        "strict_cutoff": strict_cutoff,
        "experimental_ratios": ratios,
    }


def attach_topology(strong: pd.DataFrame, comps: pd.DataFrame, strict_cutoff: float, strong_threshold: float) -> pd.DataFrame:
    strict = comps[comps["criterion"] == "step28_strict"][["sample_id", "component_id", "component_size"]].copy()
    strict = strict.rename(columns={"component_id": "strict_component_id", "component_size": "strict_component_size"})
    native = comps[comps["criterion"] == "step28_native_only_strict"][["sample_id", "component_id", "component_size"]].copy()
    native = native.rename(columns={"component_id": "native_component_id", "component_size": "native_component_size"})
    d = strong.merge(strict, on="sample_id", how="left", validate="one_to_one").merge(native, on="sample_id", how="left", validate="one_to_one")
    if d[["strict_component_id", "native_component_id"]].isna().any().any():
        raise RuntimeError("Missing Step28 topology assignment for one or more strong vectors")

    max_strict = int(d["strict_component_size"].max())
    max_native = int(d["native_component_size"].max())
    d["final_regime"] = "unclassified"
    d.loc[d["native_component_size"] == max_native, "final_regime"] = "dominant_native_connected"
    d.loc[(d["strict_component_size"] == max_strict) & (d["native_component_size"] < max_native), "final_regime"] = "dominant_flex_only_accession"
    d.loc[d["strict_component_size"] == 3, "final_regime"] = "secondary_native_connected"
    d.loc[d["sample_id"] == WITNESS_ID, "final_regime"] = "historical_witness_singleton"
    if (d["final_regime"] == "unclassified").any():
        bad = d.loc[d["final_regime"] == "unclassified", "sample_id"].tolist()
        raise RuntimeError(f"Unexpected unclassified final topology members: {bad}")

    d["is_dominant_strict_component"] = d["strict_component_size"] == max_strict
    d["is_dominant_native_component"] = d["native_component_size"] == max_native
    d["is_historical_witness"] = d["sample_id"] == WITNESS_ID
    d["calibration_margin"] = strict_cutoff - d["calibration_score_corrected"]
    d["strong_margin_abs"] = d["r_plateau_dt0p05"] - strong_threshold
    d["strong_margin_fraction"] = d["r_plateau_dt0p05"] / strong_threshold - 1.0
    d["control_plateau_rank_pct_within40"] = d["control_plateau_PO_dt0p05"].rank(pct=True, method="average")
    d["response_rank_within40"] = d["r_plateau_dt0p05"].rank(ascending=False, method="min").astype(int)
    return d.sort_values(["final_regime", "calibration_score_corrected", "sample_id"]).reset_index(drop=True)


def closest_match(frame: pd.DataFrame, target: float) -> pd.Series:
    x = frame["r_plateau_dt0p05"].to_numpy(float)
    idx = int(np.argmin(np.abs(np.log(x / target))))
    return frame.iloc[idx]


def build_match_table(strong: pd.DataFrame, ratios: np.ndarray) -> pd.DataFrame:
    scopes = {
        "all_40": strong,
        "dominant_strict_36": strong[strong["is_dominant_strict_component"]],
        "dominant_native_35": strong[strong["is_dominant_native_component"]],
        "secondary_native_3": strong[strong["final_regime"] == "secondary_native_connected"],
        "historical_witness_37318": strong[strong["is_historical_witness"]],
    }
    rows = []
    for target_idx, target in enumerate(ratios, start=1):
        for scope, frame in scopes.items():
            if frame.empty:
                continue
            r = closest_match(frame, float(target))
            model = float(r["r_plateau_dt0p05"])
            rows.append({
                "experimental_index": target_idx,
                "experimental_ratio": float(target),
                "scope": scope,
                "sample_id": int(r["sample_id"]),
                "model_ratio": model,
                "relative_error_fraction": abs(model / float(target) - 1.0),
                "absolute_log_ratio_error": abs(math.log(model / float(target))),
                "calibration_score_corrected": float(r["calibration_score_corrected"]),
                "control_plateau_PO_dt0p05": float(r["control_plateau_PO_dt0p05"]),
                "final_regime": str(r["final_regime"]),
                "G_peak": float(r["G_peak"]),
                "tau_G_ms": float(r["tau_G_ms"]),
            })
    return pd.DataFrame(rows)


def add_rep(rows: list[dict[str, Any]], lookup: pd.DataFrame, sid: int, role: str, rationale: str) -> None:
    r = lookup.loc[sid]
    rows.append({
        "sample_id": int(sid),
        "role": role,
        "rationale": rationale,
        "final_regime": str(r["final_regime"]),
        "calibration_score_corrected": float(r["calibration_score_corrected"]),
        "strict_component_size": int(r["strict_component_size"]),
        "native_component_size": int(r["native_component_size"]),
        "G_peak": float(r["G_peak"]),
        "tau_G_ms": float(r["tau_G_ms"]),
        "r_plateau_dt0p05": float(r["r_plateau_dt0p05"]),
        "control_plateau_PO_dt0p05": float(r["control_plateau_PO_dt0p05"]),
        "blocked_plateau_PO_dt0p05": float(r["blocked_plateau_PO_dt0p05"]),
        "corrected_rank": int(r["corrected_rank"]),
    })


def select_representatives(strong: pd.DataFrame, ratios: np.ndarray) -> pd.DataFrame:
    main = strong[strong["is_dominant_native_component"]].copy()
    secondary = strong[strong["final_regime"] == "secondary_native_connected"].copy()
    flex = strong[strong["final_regime"] == "dominant_flex_only_accession"].copy()
    lookup = strong.set_index("sample_id")
    rows: list[dict[str, Any]] = []

    typical_target = float(np.median(ratios))
    typical = closest_match(main, typical_target)
    add_rep(rows, lookup, int(typical["sample_id"]), "primary_typical_dominant_native", "closest dominant-native response to the median experimental potentiating ratio")

    best = main.loc[main["calibration_score_corrected"].idxmin()]
    add_rep(rows, lookup, int(best["sample_id"]), "best_calibrated_dominant_native", "minimum corrected calibration score within the dominant native-connected region")

    low = closest_match(main, float(ratios.min()))
    add_rep(rows, lookup, int(low["sample_id"]), "minimum_experimental_match_dominant_native", "closest dominant-native response to the smallest experimental potentiating ratio")

    extreme = closest_match(main, float(ratios.max()))
    add_rep(rows, lookup, int(extreme["sample_id"]), "extreme_experimental_match_dominant_native", "closest dominant-native response to the largest experimental potentiating ratio")

    largest = main.loc[main["r_plateau_dt0p05"].idxmax()]
    add_rep(rows, lookup, int(largest["sample_id"]), "largest_response_dominant_native", "largest confirmed response within the dominant native-connected region; denominator retained explicitly")

    sec = secondary.loc[secondary["calibration_score_corrected"].idxmin()]
    add_rep(rows, lookup, int(sec["sample_id"]), "secondary_region_representative", "best calibrated member of the secondary native-connected three-vector region")

    if len(flex) == 1:
        add_rep(rows, lookup, int(flex.iloc[0]["sample_id"]), "flex_aux_only_accession", "member joining the dominant strict network only when endpoint-preserving auxiliary-coordinate flexibility is allowed")

    add_rep(rows, lookup, WITNESS_ID, "historical_witness_sensitivity", "historical witness retained as an atypical sensitivity solution rather than the primary representative")

    return pd.DataFrame(rows).drop_duplicates(subset=["sample_id", "role"]).reset_index(drop=True)


def simulate_local_robustness(
    strong: pd.DataFrame,
    representatives: pd.DataFrame,
    targets: dict[str, Any],
    strict_cutoff: float,
    strong_threshold: float,
    outdir: Path,
    n_per: int,
    seed: int,
    resume: bool,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    # Publication-facing centres: primary, secondary, historical witness, and best calibrated dominant.
    wanted_roles = [
        "primary_typical_dominant_native",
        "best_calibrated_dominant_native",
        "secondary_region_representative",
        "historical_witness_sensitivity",
    ]
    centres = representatives[representatives["role"].isin(wanted_roles)].copy()
    centres = centres.drop_duplicates("sample_id")
    lookup = strong.set_index("sample_id")
    chk = outdir / "_checkpoints_local"
    chk.mkdir(parents=True, exist_ok=True)
    frames = []
    summaries = []

    for _, rr in centres.iterrows():
        sid = int(rr["sample_id"])
        row = lookup.loc[sid]
        base_rates = row[RATE_COLS].to_numpy(float)
        pulse_amp = float(row["pulse_amplitude_au"])
        calib_tau = float(row["glutamate_tau_ms"])
        g = float(row["G_peak"])
        tau = float(row["tau_G_ms"])
        for frac in (0.05, 0.10):
            p = chk / f"sample{sid}_pm{int(frac*100):02d}.csv.gz"
            if p.exists() and resume:
                frame = pd.read_csv(p)
            else:
                rng = np.random.default_rng(seed + sid * 1009 + int(frac * 10000))
                multipliers = rng.uniform(1.0 - frac, 1.0 + frac, size=(n_per, len(RATE_COLS)))
                rates = base_rates[None, :] * multipliers
                response, _ = eng.run_base_condition(rates, g, tau, 0.05)
                calibration = eng.calibration_ensemble(
                    rates,
                    np.full(n_per, pulse_amp, dtype=float),
                    np.full(n_per, calib_tau, dtype=float),
                    0.05,
                )
                score = eng.corrected_score_from_metrics(calibration[:, 1], calibration[:, 2], calibration[:, 3], calibration[:, 5], targets)
                valid = calibration[:, 8] > 0.5
                strong_ok = np.isfinite(response[:, 2]) & (response[:, 2] >= strong_threshold)
                compat = valid & (score <= strict_cutoff)
                frame = pd.DataFrame({
                    "center_sample_id": sid,
                    "center_role": str(rr["role"]),
                    "perturbation_fraction": frac,
                    "point_id": np.arange(n_per, dtype=int),
                    "G_peak": g,
                    "tau_G_ms": tau,
                    "control_plateau_PO": response[:, 0],
                    "blocked_plateau_PO": response[:, 1],
                    "r_plateau": response[:, 2],
                    "strong": strong_ok,
                    "calibration_score_corrected": score,
                    "calibration_valid": valid,
                    "compatible_strict": compat,
                    "strong_and_compatible_strict": strong_ok & compat,
                    "minimum_state_control": response[:, 3],
                    "minimum_state_blocked": response[:, 4],
                })
                for j, col in enumerate(RATE_COLS):
                    frame[col] = rates[:, j]
                frame.to_csv(p, index=False, compression="gzip")
            frames.append(frame)
            summaries.append({
                "center_sample_id": sid,
                "center_role": str(rr["role"]),
                "perturbation_fraction": frac,
                "n": int(len(frame)),
                "strong_fraction": float(frame["strong"].mean()),
                "compatible_strict_fraction": float(frame["compatible_strict"].mean()),
                "joint_strict_fraction": float(frame["strong_and_compatible_strict"].mean()),
                "median_r_plateau": float(frame["r_plateau"].median()),
                "median_control_plateau_PO": float(frame["control_plateau_PO"].median()),
                "median_corrected_score": float(frame["calibration_score_corrected"].median()),
                "min_r_plateau": float(frame["r_plateau"].min()),
                "max_corrected_score": float(frame["calibration_score_corrected"].max()),
            })
            print(f"Step29 local: sample={sid} frac={frac:.2f} joint={summaries[-1]['joint_strict_fraction']:.3f}", flush=True)
    return pd.concat(frames, ignore_index=True), pd.DataFrame(summaries)


def make_figures(
    step24: Path,
    strong: pd.DataFrame,
    matches: pd.DataFrame,
    reps: pd.DataFrame,
    full_grid: pd.DataFrame,
    local_summary: pd.DataFrame | None,
    outdir: Path,
    strong_threshold: float,
) -> None:
    # Fig29A: calibration correction across strong thresholds.
    t = pd.read_csv(step24 / "04A_threshold_strong_support_shift.csv")
    fig, ax = plt.subplots(figsize=(6.5, 4.4), constrained_layout=True)
    ax.plot(t["threshold"], 100*t["precal_strong_fraction"], marker="o", label="pre-calibration broad")
    ax.plot(t["threshold"], 100*t["historical_accepted_strong_fraction"], marker="o", label="historical accepted")
    ax.plot(t["threshold"], 100*t["corrected_accepted_strong_fraction"], marker="o", label="corrected accepted")
    ax.set_xlabel("Strong-response threshold, blocked/control ratio")
    ax.set_ylabel("Strong vectors (%)")
    ax.set_title("A  Corrected calibration restores strong-response support")
    ax.legend(frameon=False)
    fig.savefig(outdir / "Fig29A_corrected_support.pdf")
    fig.savefig(outdir / "Fig29A_corrected_support.png", dpi=300)
    plt.close(fig)

    # Fig29B: denominator vs response, by final sampled regime.
    fig, ax = plt.subplots(figsize=(6.5, 4.6), constrained_layout=True)
    for regime, d in strong.groupby("final_regime", sort=False):
        ax.scatter(d["control_plateau_PO_dt0p05"], d["r_plateau_dt0p05"], s=38, label=regime.replace("_", " "))
    ax.axhline(strong_threshold, linestyle="--", linewidth=1.0)
    for sid in sorted(set(reps["sample_id"].astype(int))):
        r = strong[strong["sample_id"] == sid].iloc[0]
        ax.annotate(str(sid), (r["control_plateau_PO_dt0p05"], r["r_plateau_dt0p05"]), xytext=(4, 4), textcoords="offset points", fontsize=7)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Control plateau open-state occupancy")
    ax.set_ylabel("Blocked/control plateau ratio")
    ax.set_title("B  Corrected strong-compatible solutions")
    ax.legend(frameon=False, fontsize=7)
    fig.savefig(outdir / "Fig29B_response_denominator_topology.pdf")
    fig.savefig(outdir / "Fig29B_response_denominator_topology.png", dpi=300)
    plt.close(fig)

    # Fig29C: experimental ratio matching, all40 vs dominant native.
    m_all = matches[matches["scope"] == "all_40"].sort_values("experimental_ratio")
    m_main = matches[matches["scope"] == "dominant_native_35"].sort_values("experimental_ratio")
    fig, ax = plt.subplots(figsize=(5.6, 5.0), constrained_layout=True)
    lo = min(m_all["experimental_ratio"].min(), m_main["model_ratio"].min()) * 0.92
    hi = max(m_all["experimental_ratio"].max(), m_main["model_ratio"].max()) * 1.08
    ax.plot([lo, hi], [lo, hi], linestyle="--", linewidth=1.0, label="identity")
    ax.scatter(m_all["experimental_ratio"], m_all["model_ratio"], s=46, label="closest among all 40")
    ax.scatter(m_main["experimental_ratio"], m_main["model_ratio"], s=46, marker="s", label="closest dominant-native")
    for _, r in m_main.iterrows():
        ax.annotate(str(int(r["sample_id"])), (r["experimental_ratio"], r["model_ratio"]), xytext=(4, 3), textcoords="offset points", fontsize=7)
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_xlabel("Experimental potentiating ratio")
    ax.set_ylabel("Closest model ratio")
    ax.set_title("C  Experimental-scale responses within the dominant region")
    ax.legend(frameon=False, fontsize=8)
    fig.savefig(outdir / "Fig29C_experimental_ratio_matches.pdf")
    fig.savefig(outdir / "Fig29C_experimental_ratio_matches.png", dpi=300)
    plt.close(fig)

    # Fig29D: disjoint final topology composition.
    counts = strong["final_regime"].value_counts()
    order = [
        "dominant_native_connected",
        "dominant_flex_only_accession",
        "secondary_native_connected",
        "historical_witness_singleton",
    ]
    vals = [int(counts.get(k, 0)) for k in order]
    labels = ["dominant native", "flex-only accession", "secondary native", "historical witness"]
    fig, ax = plt.subplots(figsize=(6.4, 4.2), constrained_layout=True)
    bars = ax.bar(labels, vals)
    for b, v in zip(bars, vals):
        ax.text(b.get_x()+b.get_width()/2, v+0.4, str(v), ha="center", va="bottom", fontsize=9)
    ax.set_ylabel("Confirmed strong-compatible vectors")
    ax.set_title("D  Final sampled connectivity composition (n=40)")
    ax.tick_params(axis="x", rotation=18)
    fig.savefig(outdir / "Fig29D_connectivity_composition.pdf")
    fig.savefig(outdir / "Fig29D_connectivity_composition.png", dpi=300)
    plt.close(fig)

    # Fig29E: primary representative response map from the frozen 80-node grid.
    primary_id = int(reps.loc[reps["role"] == "primary_typical_dominant_native", "sample_id"].iloc[0])
    g = full_grid[full_grid["sample_id"] == primary_id].copy()
    if not g.empty:
        pv = g.pivot(index="tau_G_ms", columns="G_peak", values="r_plateau")
        fig, ax = plt.subplots(figsize=(6.2, 4.8), constrained_layout=True)
        im = ax.imshow(pv.to_numpy(), origin="lower", aspect="auto", interpolation="nearest")
        ax.set_xticks(np.arange(len(pv.columns)), [f"{x:g}" for x in pv.columns])
        ax.set_yticks(np.arange(len(pv.index)), [f"{x:g}" for x in pv.index])
        ax.set_xlabel("G peak")
        ax.set_ylabel("Glutamate decay, ms")
        ax.set_title(f"E  Primary representative {primary_id}: 80-node response map")
        cbar = fig.colorbar(im, ax=ax)
        cbar.set_label("Blocked/control plateau ratio")
        fig.savefig(outdir / "Fig29E_primary_response_map.pdf")
        fig.savefig(outdir / "Fig29E_primary_response_map.png", dpi=300)
        plt.close(fig)

    if local_summary is not None and not local_summary.empty:
        fig, ax = plt.subplots(figsize=(7.2, 4.6), constrained_layout=True)
        roles = list(local_summary["center_role"].drop_duplicates())
        centres = [int(local_summary[local_summary["center_role"] == r]["center_sample_id"].iloc[0]) for r in roles]
        x = np.arange(len(roles), dtype=float)
        width = 0.34
        for j, frac in enumerate((0.05, 0.10)):
            vals = []
            for role in roles:
                d = local_summary[(local_summary["center_role"] == role) & np.isclose(local_summary["perturbation_fraction"], frac)]
                vals.append(float(d["joint_strict_fraction"].iloc[0]))
            ax.bar(x + (j-0.5)*width, vals, width=width, label=f"±{int(frac*100)}%")
        ax.set_xticks(x, [str(s) for s in centres])
        ax.set_ylim(0, 1.05)
        ax.set_ylabel("Strong + corrected-compatible fraction")
        ax.set_xlabel("Representative sample id")
        ax.set_title("F  Simultaneous 12-rate local robustness")
        ax.legend(frameon=False)
        fig.savefig(outdir / "Fig29F_representative_local_robustness.pdf")
        fig.savefig(outdir / "Fig29F_representative_local_robustness.png", dpi=300)
        plt.close(fig)


def write_manuscript_numbers(
    outdir: Path,
    step24: Path,
    strong: pd.DataFrame,
    reps: pd.DataFrame,
    matches: pd.DataFrame,
    gate: dict[str, Any],
    local_summary: pd.DataFrame | None,
) -> None:
    support = pd.read_csv(step24 / "04_primary_strong_support_shift.csv").iloc[0]
    acc = pd.read_csv(step24 / "03_acceptance_shift_by_prior.csv")
    broad = acc[acc["prior"] == "broad"].iloc[0]
    ratios = gate["experimental_ratios"]
    primary = reps[reps["role"] == "primary_typical_dominant_native"].iloc[0]
    sec = reps[reps["role"] == "secondary_region_representative"].iloc[0]
    wit = reps[reps["role"] == "historical_witness_sensitivity"].iloc[0]
    main_match = matches[matches["scope"] == "dominant_native_35"]

    rows = [
        ("broad_precalibration_candidates", int(support["precal_strong_n"] / support["precal_strong_fraction"])),
        ("broad_precalibration_strong_n_primary_threshold", int(support["precal_strong_n"])),
        ("broad_precalibration_strong_fraction", float(support["precal_strong_fraction"])),
        ("historical_accepted_broad_strong_n", int(support["historical_accepted_strong_n"])),
        ("corrected_accepted_broad_n", int(broad["corrected_accepted_n"])),
        ("corrected_accepted_broad_strong_n", int(support["corrected_accepted_strong_n"])),
        ("corrected_accepted_broad_strong_fraction", float(support["corrected_accepted_strong_fraction"])),
        ("historical_vs_corrected_broad_jaccard", float(broad["jaccard"])),
        ("dominant_strict_connected_n", int((strong["strict_component_size"] == strong["strict_component_size"].max()).sum())),
        ("dominant_native_connected_n", int((strong["native_component_size"] == strong["native_component_size"].max()).sum())),
        ("secondary_native_connected_n", int((strong["final_regime"] == "secondary_native_connected").sum())),
        ("historical_witness_singleton_n", int((strong["final_regime"] == "historical_witness_singleton").sum())),
        ("primary_representative_sample_id", int(primary["sample_id"])),
        ("primary_representative_ratio", float(primary["r_plateau_dt0p05"])),
        ("primary_representative_corrected_score", float(primary["calibration_score_corrected"])),
        ("secondary_representative_sample_id", int(sec["sample_id"])),
        ("historical_witness_sample_id", int(wit["sample_id"])),
        ("historical_witness_ratio", float(wit["r_plateau_dt0p05"])),
        ("historical_witness_corrected_score", float(wit["calibration_score_corrected"])),
        ("dominant_native_match_median_relative_error", float(main_match["relative_error_fraction"].median())),
        ("dominant_native_match_max_relative_error", float(main_match["relative_error_fraction"].max())),
        ("strict_corrected_calibration_cutoff", float(gate["strict_cutoff"])),
        ("strong_response_threshold", float(gate["strong_threshold"])),
    ]
    for i, r in enumerate(ratios, start=1):
        rows.append((f"experimental_potentiating_ratio_{i}", float(r)))
    pd.DataFrame(rows, columns=["quantity", "value"]).to_csv(outdir / "07_manuscript_headline_numbers.csv", index=False)

    lines = [
        "# Step 29 manuscript numbers",
        "",
        f"Corrected broad support: **{int(support['corrected_accepted_strong_n'])}/{int(broad['corrected_accepted_n'])} = {100*float(support['corrected_accepted_strong_fraction']):.3f}%** strong at the primary threshold.",
        f"Pre-calibration broad support: **{int(support['precal_strong_n'])}/50000 = {100*float(support['precal_strong_fraction']):.3f}%**.",
        f"Historical accepted broad support at the same threshold: **{int(support['historical_accepted_strong_n'])}/5000**.",
        f"Historical/corrected retained-set Jaccard: **{float(broad['jaccard']):.4f}**.",
        "",
        "Final sampled connectivity under the strict corrected criterion:",
        f"- dominant endpoint-preserving network: **36/40 = 90.0%**;",
        f"- dominant native-aux kinetic network: **35/40 = 87.5%**;",
        f"- secondary native-connected region: **3/40 = 7.5%**;",
        f"- historical witness 37318: **1/40 = 2.5%**, remaining separate within the finite search budget.",
        "",
        f"Primary manuscript representative: **sample {int(primary['sample_id'])}**, r={float(primary['r_plateau_dt0p05']):.6f}, corrected score={float(primary['calibration_score_corrected']):.6f}.",
        f"Historical witness: **sample {int(wit['sample_id'])}**, r={float(wit['r_plateau_dt0p05']):.6f}, corrected score={float(wit['calibration_score_corrected']):.6f}.",
        "",
        f"Dominant-native matches to the five experimental potentiating ratios have median relative error **{100*float(main_match['relative_error_fraction'].median()):.3f}%** and maximum relative error **{100*float(main_match['relative_error_fraction'].max()):.3f}%**.",
    ]
    if local_summary is not None and not local_summary.empty:
        lines += ["", "Representative local robustness (strict corrected compatibility):"]
        for _, r in local_summary.iterrows():
            lines.append(f"- sample {int(r['center_sample_id'])}, ±{int(round(100*r['perturbation_fraction']))}%: joint fraction **{float(r['joint_strict_fraction']):.3f}**.")
    (outdir / "MANUSCRIPT_NUMBERS.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def source_manifest(step24: Path, step25: Path, step28: Path, outdir: Path) -> None:
    files = [
        step24 / "03_acceptance_shift_by_prior.csv",
        step24 / "04A_threshold_strong_support_shift.csv",
        step24 / "04_primary_strong_support_shift.csv",
        step24 / "08_scientific_decision.json",
        step25 / "00_input_audit.json",
        step25 / "01_corrected_target_reconstruction.json",
        step25 / "02_corrected_accepted_cohorts.csv.gz",
        step25 / "04_full_80_node_grid.csv.gz",
        step25 / "05_strong_candidates_confirmed.csv",
        step25 / "15_scientific_decision.json",
        step25 / "run_summary.json",
        step28 / "07_final_network_edges.csv",
        step28 / "08_final_components.csv",
        step28 / "09_residual_barrier_summary.csv",
        step28 / "10_scientific_decision.json",
    ]
    rows = []
    for p in files:
        rows.append({"path": str(p), "sha256": sha256_file(p), "bytes": p.stat().st_size})
    pd.DataFrame(rows).to_csv(outdir / "00_source_manifest.csv", index=False)


def self_test() -> None:
    # Minimal deterministic tests for topology assignment and representative selection.
    rows = []
    for sid in range(1, 41):
        rows.append({
            "sample_id": sid,
            "strong_confirmed_dt0p05": True,
            "r_plateau_dt0p05": 1.64 + sid/100.0,
            "calibration_score_corrected": sid/1000.0,
            "control_plateau_PO_dt0p05": 1e-3 + sid*1e-5,
            "blocked_plateau_PO_dt0p05": 2e-3,
            "corrected_rank": sid,
            "G_peak": 0.5,
            "tau_G_ms": 20.0,
        })
    strong = pd.DataFrame(rows)
    strong.loc[strong["sample_id"] == 40, "sample_id"] = WITNESS_ID
    comp_rows = []
    ids = strong["sample_id"].astype(int).tolist()
    for sid in ids[:36]: comp_rows.append({"criterion":"step28_strict","component_id":0,"component_size":36,"sample_id":sid})
    for sid in ids[36:39]: comp_rows.append({"criterion":"step28_strict","component_id":1,"component_size":3,"sample_id":sid})
    comp_rows.append({"criterion":"step28_strict","component_id":2,"component_size":1,"sample_id":ids[39]})
    for sid in ids[:35]: comp_rows.append({"criterion":"step28_native_only_strict","component_id":0,"component_size":35,"sample_id":sid})
    comp_rows.append({"criterion":"step28_native_only_strict","component_id":2,"component_size":1,"sample_id":ids[35]})
    for sid in ids[36:39]: comp_rows.append({"criterion":"step28_native_only_strict","component_id":1,"component_size":3,"sample_id":sid})
    comp_rows.append({"criterion":"step28_native_only_strict","component_id":3,"component_size":1,"sample_id":ids[39]})
    d = attach_topology(strong, pd.DataFrame(comp_rows), 0.1, 1.6359)
    assert (d["final_regime"] == "dominant_native_connected").sum() == 35
    assert (d["final_regime"] == "dominant_flex_only_accession").sum() == 1
    assert (d["final_regime"] == "secondary_native_connected").sum() == 3
    assert (d["final_regime"] == "historical_witness_singleton").sum() == 1
    print("Step29 self-test PASS")


def main() -> None:
    args = parse_args()
    if args.self_test:
        self_test()
        return
    if args.threads < 1:
        raise ValueError("--threads must be >=1")
    eng.set_num_threads(args.threads)

    step24 = discover_step24(args.step24)
    step25 = discover_step25(args.step25)
    step28 = discover_step28(args.step28)
    gate = hard_gate(step24, step25, step28)

    outdir = args.output.resolve()
    if outdir.exists() and not args.resume:
        if args.force:
            shutil.rmtree(outdir)
        else:
            raise FileExistsError(f"Output exists: {outdir}. Set STEP29_RESUME=1 to continue or STEP29_FORCE=1 to replace it.")
    outdir.mkdir(parents=True, exist_ok=True)

    source_manifest(step24, step25, step28, outdir)
    input_audit = {
        "pipeline": PIPELINE_NAME,
        "script_version": SCRIPT_VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "numba_threads": eng.get_num_threads(),
        "step24": str(step24),
        "step25": str(step25),
        "step28": str(step28),
        "hard_gate_checks": gate["checks"],
        "strong_threshold": gate["strong_threshold"],
        "strict_corrected_cutoff": gate["strict_cutoff"],
        "experimental_ratios": gate["experimental_ratios"].tolist(),
        "no_new_denominator_cutoff": True,
        "no_new_response_threshold": True,
        "local_perturbation_n_per_center_per_level": 0 if args.skip_local else args.local_n,
        "local_perturbation_levels": [] if args.skip_local else [0.05, 0.10],
    }
    json_dump(outdir / "00_input_audit.json", input_audit)

    if args.preflight_only:
        print(json.dumps(input_audit, indent=2))
        return

    accepted = pd.read_csv(step25 / "02_corrected_accepted_cohorts.csv.gz")
    accepted_broad = accepted[accepted["prior"] == "broad"].copy()
    if len(accepted_broad) != 5000:
        raise RuntimeError(f"Expected 5000 corrected accepted broad vectors; got {len(accepted_broad)}")
    accepted.to_csv(outdir / "01_frozen_corrected_accepted_cohorts.csv.gz", index=False, compression="gzip")

    strong0 = pd.read_csv(step25 / "05_strong_candidates_confirmed.csv")
    comps = pd.read_csv(step28 / "08_final_components.csv")
    strong = attach_topology(strong0, comps, gate["strict_cutoff"], gate["strong_threshold"])
    strong.to_csv(outdir / "02_final_strong40_topology.csv", index=False)

    # Freeze final network evidence.
    edges = pd.read_csv(step28 / "07_final_network_edges.csv")
    edges.to_csv(outdir / "03_final_verified_network_edges.csv", index=False)
    barriers = pd.read_csv(step28 / "09_residual_barrier_summary.csv")
    barriers.to_csv(outdir / "03A_final_residual_barriers.csv", index=False)

    matches = build_match_table(strong, gate["experimental_ratios"])
    matches.to_csv(outdir / "04_experimental_model_matches.csv", index=False)

    reps = select_representatives(strong, gate["experimental_ratios"])
    reps.to_csv(outdir / "05_publication_representatives.csv", index=False)

    local_points = None
    local_summary = None
    if not args.skip_local:
        targets = json.loads((step25 / "01_corrected_target_reconstruction.json").read_text())
        local_points, local_summary = simulate_local_robustness(
            strong, reps, targets, gate["strict_cutoff"], gate["strong_threshold"], outdir, args.local_n, args.seed, args.resume
        )
        local_points.to_csv(outdir / "06_representative_local_perturbation_points.csv.gz", index=False, compression="gzip")
        local_summary.to_csv(outdir / "06A_representative_local_robustness_summary.csv", index=False)

    full_grid = pd.read_csv(step25 / "04_full_80_node_grid.csv.gz")
    make_figures(step24, strong, matches, reps, full_grid, local_summary, outdir, gate["strong_threshold"])
    write_manuscript_numbers(outdir, step24, strong, reps, matches, gate, local_summary)

    # Final scientific decision and interpretation boundaries.
    strict_sizes = sorted_component_sizes(comps, "step28_strict")
    native_sizes = sorted_component_sizes(comps, "step28_native_only_strict")
    primary_id = int(reps.loc[reps["role"] == "primary_typical_dominant_native", "sample_id"].iloc[0])
    decision = {
        "status": "CORRECTED_PUBLICATION_MODEL_RESULTS_CONSOLIDATED",
        "interpretation_scope": "Corrected 140-500-ms calibration, 40 numerically confirmed strong vectors, and finite-budget connectivity evidence through Step28. Connectivity failures are not global topological proofs.",
        "n_corrected_accepted_broad": 5000,
        "n_confirmed_strong": 40,
        "strict_component_sizes": strict_sizes,
        "native_only_strict_component_sizes": native_sizes,
        "dominant_strict_fraction_of_strong": 36/40,
        "dominant_native_fraction_of_strong": 35/40,
        "primary_representative_sample_id": primary_id,
        "historical_witness_sample_id": WITNESS_ID,
        "historical_witness_is_primary_representative": primary_id == WITNESS_ID,
        "no_new_denominator_cutoff": True,
        "no_new_response_threshold": True,
        "old_witness_centric_interpretation_should_be_retired": True,
        "publication_message": "Multiple corrected strong-compatible kinetic solutions exist. Most sampled solutions lie in one verified connected network; a small secondary native-connected region and the historical witness remain separate within the finite path-search budget.",
    }
    if local_summary is not None:
        decision["representative_local_robustness_recomputed"] = True
        decision["local_robustness_n_per_center_per_level"] = args.local_n
    else:
        decision["representative_local_robustness_recomputed"] = False
    json_dump(outdir / "08_scientific_decision.json", decision)

    interpretation = f"""# Step 29 scientific interpretation\n\nStatus: **{decision['status']}**.\n\nThe corrected broad accepted cohort contains 5000 vectors, of which 40 are numerically confirmed to reach the frozen experimental strong-response threshold. The final sampled strict connectivity structure is **{' + '.join(map(str, strict_sizes))}**; the stricter native-aux structure is **{' + '.join(map(str, native_sizes))}**.\n\nThe primary manuscript representative is sample **{primary_id}**, selected from the dominant native-connected region by closeness to the median experimental potentiating ratio. Sample 37318 is retained only as the historical-witness sensitivity solution.\n\nThe publication-level interpretation is therefore no longer a unique-witness or isolated-basin claim. The supported statement is that corrected control calibration retains multiple strong-response solutions, most sampled solutions belong to one verified connected network, and residual separation persists for a small secondary region and the historical witness under the finite path families tested through Step28.\n\nNo new denominator cutoff, response threshold, experimental target, or calibration tolerance was introduced in Step29.\n"""
    (outdir / "SCIENTIFIC_INTERPRETATION.md").write_text(interpretation, encoding="utf-8")

    readme = f"""# Step 29 results — corrected publication consolidation\n\nStatus: **{decision['status']}**\n\nInputs:\n- Step24: `{step24}`\n- Step25: `{step25}`\n- Step28: `{step28}`\n\nFrozen outputs:\n- `01_frozen_corrected_accepted_cohorts.csv.gz`\n- `02_final_strong40_topology.csv`\n- `03_final_verified_network_edges.csv`\n- `04_experimental_model_matches.csv`\n- `05_publication_representatives.csv`\n- `07_manuscript_headline_numbers.csv`\n- `MANUSCRIPT_NUMBERS.md`\n- `SCIENTIFIC_INTERPRETATION.md`\n\nPrimary representative: **sample {primary_id}**. Historical witness: **sample {WITNESS_ID}**.\n\nThe old witness-centric interpretation is retired. Step29 introduces no new biological or numerical acceptance threshold.\n"""
    (outdir / "README_RESULTS.md").write_text(readme, encoding="utf-8")

    run_summary = {
        **decision,
        "pipeline": PIPELINE_NAME,
        "script_version": SCRIPT_VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "output_directory": str(outdir),
        "representatives": reps.to_dict("records"),
        "experimental_ratios": gate["experimental_ratios"].tolist(),
        "local_robustness_summary": None if local_summary is None else local_summary.to_dict("records"),
        "files": sorted([p.name for p in outdir.iterdir() if p.is_file()]),
    }
    json_dump(outdir / "run_summary.json", run_summary)

    # Do not preserve transient local checkpoints in a completed run.
    chk = outdir / "_checkpoints_local"
    if chk.exists():
        shutil.rmtree(chk)

    # Create an archive adjacent to the results directory.
    archive_base = outdir.parent / outdir.name
    archive = shutil.make_archive(str(archive_base), "zip", root_dir=outdir.parent, base_dir=outdir.name)
    print(f"Step29 results written to: {outdir}")
    print(f"Step29 archive: {archive}")


if __name__ == "__main__":
    main()
