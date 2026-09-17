#!/usr/bin/env python3
"""Step 25: corrected strong-ensemble validation and kinetic geometry.

This pipeline starts from the Step-24 correction of the control-calibration
integral (140--500 ms) and asks whether the corrected strong compatible Base
solutions form one sampled kinetic region or several distinct regimes.

The analysis deliberately keeps three calculations separate:
1. historical Step-10 control calibration, replayed with the original source
   integrator and the corrected J(140--500) score;
2. final Base-model strong-response testing on the frozen 80-node forcing grid;
3. geometry/local perturbation diagnostics in the 12 kinetic-rate coordinates.

No new denominator threshold is introduced. Small control plateaus are reported
as a diagnostic, not used as an exclusion rule.
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
import zipfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from numba import get_num_threads, njit, prange, set_num_threads
from scipy.spatial.distance import cdist
from scipy.sparse.csgraph import minimum_spanning_tree
from sklearn.cluster import OPTICS
from sklearn.metrics import adjusted_rand_score


SCRIPT_VERSION = "1.0.0"
STRONG_DEFAULT = 1.6359295439505144
EXPERIMENTAL_RATIOS = np.array(
    [1.6359295439505144, 1.694672, 1.876748, 1.891305, 4.072559],
    dtype=float,
)
G_PEAK_GRID = np.array([0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50, 0.65, 0.80, 0.95], dtype=float)
TAU_GRID = np.array([3.0, 5.0, 7.5, 10.0, 15.0, 20.0, 30.0, 40.0], dtype=float)
PULSE_COUNT = 25
PULSE_INTERVAL_MS = 5.1
TRAIN_END_MS = (PULSE_COUNT - 1) * PULSE_INTERVAL_MS

RATE_COLS = [
    "kon_A_per_ms",
    "kon_B_per_ms",
    "koff_A_per_ms",
    "koff_B_per_ms",
    "kf_plus_per_ms",
    "ks_plus_per_ms",
    "kf_minus_per_ms",
    "ks_minus_per_ms",
    "kd1_plus_per_ms",
    "kd1_minus_per_ms",
    "kd2_plus_per_ms",
    "kd2_minus_per_ms",
]
AUX_COLS = ["pulse_amplitude_au", "glutamate_tau_ms"]

# Numeric indices for the 12-rate vectors.
KON_A = 0
KON_B = 1
KOFF_A = 2
KOFF_B = 3
KF_PLUS = 4
KS_PLUS = 5
KF_MINUS = 6
KS_MINUS = 7
KD1_PLUS = 8
KD1_MINUS = 9
KD2_PLUS = 10
KD2_MINUS = 11


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Step 25 corrected strong-ensemble validation and geometry")
    p.add_argument("--step24", default="auto", help="Step-24 results directory or 'auto'")
    p.add_argument("--step10", default="auto", help="Step-10 full candidate CSV/ZIP or 'auto'")
    p.add_argument(
        "--output",
        type=Path,
        default=Path("/root/nmda2/step_25/results_step_25_corrected_strong_geometry"),
    )
    p.add_argument("--threads", type=int, default=2)
    p.add_argument("--local-n", type=int, default=1000)
    p.add_argument("--seed", type=int, default=20260913)
    p.add_argument("--max-centres", type=int, default=8)
    p.add_argument("--grid-dt-ms", type=float, default=0.1)
    p.add_argument("--confirm-dt-ms", type=float, default=0.05)
    p.add_argument("--calibration-dt-ms", type=float, default=0.05)
    p.add_argument("--screen-replay-rel-tol", type=float, default=0.01)
    p.add_argument("--screen-replay-abs-tol", type=float, default=0.01)
    p.add_argument("--allow-replay-mismatch", action="store_true")
    p.add_argument("--force", action="store_true")
    p.add_argument("--preflight-only", action="store_true")
    p.add_argument("--self-test", action="store_true")
    return p.parse_args()


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def json_dump(path: Path, data: Any) -> None:
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def condition_id(g: float, tau: float) -> str:
    return f"G{int(round(g * 1000)):03d}_T{int(round(tau * 10)):03d}"


def discover_step24(spec: str) -> Path:
    if spec != "auto":
        p = Path(spec).expanduser().resolve()
        if not p.is_dir():
            raise FileNotFoundError(p)
        return p
    candidates = [
        Path("/root/nmda2/step_24/results_step_24_calibration_window_correction"),
        Path("/root/nmda2/results_step_24_calibration_window_correction"),
        Path("/root/nmda2/NMDAR_Braess_step24_calibration_window_correction_v1_0/results_step_24_calibration_window_correction"),
    ]
    required = {"02_corrected_candidate_scores.csv.gz", "08_scientific_decision.json"}
    for p in candidates:
        if p.is_dir() and required.issubset({x.name for x in p.iterdir()}):
            return p.resolve()
    root = Path("/root/nmda2")
    if root.exists():
        for marker in root.rglob("08_scientific_decision.json"):
            p = marker.parent
            if (p / "02_corrected_candidate_scores.csv.gz").exists() and (p / "04B_corrected_accepted_top_response_candidates.csv").exists():
                return p.resolve()
    raise FileNotFoundError("Could not auto-discover Step-24 results")


def discover_step10(spec: str, step24: Path) -> tuple[Path, str | None]:
    if spec != "auto":
        p = Path(spec).expanduser().resolve()
        if not p.exists():
            raise FileNotFoundError(p)
        return p, None
    audit = json.loads((step24 / "00_input_audit.json").read_text(encoding="utf-8"))
    source = str(audit.get("step10_source", ""))
    if "::" in source:
        left, member = source.split("::", 1)
        p = Path(left)
        if p.exists():
            return p.resolve(), member
    candidates = [Path("/root/nmda/results_step10.zip"), Path("/root/nmda/results_step10_tables.zip")]
    for p in candidates:
        if p.exists():
            return p.resolve(), None
    raise FileNotFoundError("Could not auto-discover full Step-10 candidate table")


def read_step10(path: Path, preferred_member: str | None) -> tuple[pd.DataFrame, str]:
    usecols = ["prior", "sample_id", *RATE_COLS, *AUX_COLS]
    if path.suffix.lower() == ".zip":
        with zipfile.ZipFile(path) as zf:
            names = zf.namelist()
            member = preferred_member if preferred_member in names else None
            if member is None:
                matches = [n for n in names if n.endswith("matched_parameter_all_samples.csv")]
                if not matches:
                    matches = [n for n in names if n.endswith(".csv") and "all_samples" in n]
                if not matches:
                    raise FileNotFoundError(f"No full candidate CSV inside {path}")
                member = sorted(matches, key=len)[0]
            with zf.open(member) as f:
                df = pd.read_csv(f, usecols=usecols)
        return df, f"{path}::{member}"
    if path.is_file():
        return pd.read_csv(path, usecols=usecols), str(path)
    matches = list(path.rglob("matched_parameter_all_samples.csv"))
    if not matches:
        raise FileNotFoundError(f"No matched_parameter_all_samples.csv below {path}")
    selected = sorted(matches, key=lambda x: len(str(x)))[0]
    return pd.read_csv(selected, usecols=usecols), str(selected)


def load_step24(step24: Path) -> dict[str, Any]:
    files = {
        "audit": "00_input_audit.json",
        "scores": "02_corrected_candidate_scores.csv.gz",
        "shift": "03_acceptance_shift_by_prior.csv",
        "support": "04_primary_strong_support_shift.csv",
        "top": "04B_corrected_accepted_top_response_candidates.csv",
        "decision": "08_scientific_decision.json",
    }
    for name in files.values():
        if not (step24 / name).exists():
            raise FileNotFoundError(step24 / name)
    audit = json.loads((step24 / files["audit"]).read_text(encoding="utf-8"))
    decision = json.loads((step24 / files["decision"]).read_text(encoding="utf-8"))
    scores = pd.read_csv(step24 / files["scores"])
    shift = pd.read_csv(step24 / files["shift"])
    support = pd.read_csv(step24 / files["support"])
    top = pd.read_csv(step24 / files["top"])
    threshold = float(audit.get("strong_threshold", support.iloc[0]["threshold"]))
    expected_strong_n = int(round(float(support.iloc[0]["corrected_accepted_strong_n"])))
    strong = top[top["r_max_confirmed_or_screen"] >= threshold].copy()
    if len(strong) != expected_strong_n:
        raise RuntimeError(
            f"Step24 strong count mismatch: top table gives {len(strong)}, support table gives {expected_strong_n}"
        )
    if decision.get("status") != "MATERIAL_CHANGE":
        raise RuntimeError(f"Step24 status is {decision.get('status')!r}, expected MATERIAL_CHANGE")
    return {
        "audit": audit,
        "decision": decision,
        "scores": scores,
        "shift": shift,
        "support": support,
        "top": top,
        "strong": strong,
        "threshold": threshold,
        "expected_strong_n": expected_strong_n,
    }


def recover_corrected_targets(scores: pd.DataFrame) -> dict[str, Any]:
    h = scores["control_plateau_PO_140_250"].to_numpy(float)
    e = scores["control_early_PO_140_200"].to_numpy(float)
    l = scores["control_late_PO_200_500"].to_numpy(float)
    j = scores["control_integrated_open_state_140_500_ms"].to_numpy(float)
    s = scores["calibration_score_corrected"].to_numpy(float)
    mask = (
        np.isfinite(h) & np.isfinite(e) & np.isfinite(l) & np.isfinite(j) & np.isfinite(s)
        & (h > 0) & (e > 0) & (l > 0) & (j > 0)
    )
    x = np.column_stack([np.log(e[mask] / h[mask]), np.log(l[mask] / h[mask]), np.log(j[mask] / h[mask])])
    weights = np.array([1.0 / 0.35**2, 1.0 / 0.35**2, 0.5 / 0.50**2], dtype=float)
    y = s[mask] - np.sum(weights * x * x, axis=1)
    design = np.column_stack([x, np.ones(len(x))])
    beta = np.linalg.lstsq(design, y, rcond=None)[0]
    mu = -beta[:3] / (2.0 * weights)
    targets = np.exp(mu)
    replay = np.sum(weights * (x - mu) ** 2, axis=1)
    max_err = float(np.max(np.abs(replay - s[mask])))
    intercept_expected = float(np.sum(weights * mu * mu))
    if max_err > 1e-8:
        raise RuntimeError(f"Corrected-score target reconstruction failed: max error {max_err:g}")
    return {
        "early_over_primary": float(targets[0]),
        "late_over_primary": float(targets[1]),
        "integral_over_primary_ms": float(targets[2]),
        "early_log_tolerance": 0.35,
        "late_log_tolerance": 0.35,
        "integral_log_tolerance": 0.50,
        "integral_penalty_weight": 0.5,
        "n_finite": int(mask.sum()),
        "max_abs_score_replay_error": max_err,
        "fitted_intercept": float(beta[3]),
        "expected_intercept": intercept_expected,
    }


def corrected_score_from_metrics(h: np.ndarray, e: np.ndarray, l: np.ndarray, j: np.ndarray, targets: dict[str, Any]) -> np.ndarray:
    result = np.full_like(h, np.inf, dtype=float)
    valid = np.isfinite(h) & np.isfinite(e) & np.isfinite(l) & np.isfinite(j) & (h > 0) & (e > 0) & (l > 0) & (j > 0)
    if np.any(valid):
        er = e[valid] / h[valid]
        lr = l[valid] / h[valid]
        jr = j[valid] / h[valid]
        result[valid] = (
            (np.log(er / targets["early_over_primary"]) / 0.35) ** 2
            + (np.log(lr / targets["late_over_primary"]) / 0.35) ** 2
            + 0.5 * (np.log(jr / targets["integral_over_primary_ms"]) / 0.50) ** 2
        )
    return result


# -----------------------------------------------------------------------------
# Historical Step-10 control replay: explicit midpoint, exactly as source logic.
# -----------------------------------------------------------------------------

@njit(cache=True)
def _calib_rhs(state: np.ndarray, g: float, rates: np.ndarray) -> np.ndarray:
    p0, p1a, p1b, p2, pd1, pd2, pc1, pc2, po, pro = state
    kon_a, kon_b, koff_a, koff_b = rates[0], rates[1], rates[2], rates[3]
    kfp, ksp, kfm, ksm = rates[4], rates[5], rates[6], rates[7]
    kd1p, kd1m, kd2p, kd2m = rates[8], rates[9], rates[10], rates[11]

    f0a = kon_a * g * p0
    f0b = kon_b * g * p0
    r_a = koff_a * p1a
    r_b = koff_b * p1b
    f_ab = kon_b * g * p1a
    r_2a = koff_b * p2
    f_ba = kon_a * g * p1b
    r_2b = koff_a * p2
    f_d1 = kd1p * p2
    r_d1 = kd1m * pd1
    f_d2 = kd2p * p2
    r_d2 = kd2m * pd2
    f_c1 = kfp * p2
    r_c1 = kfm * pc1
    f_c2 = ksp * p2
    r_c2 = ksm * pc2
    c1_o = ksp * pc1
    o_c1 = ksm * po
    c2_o = kfp * pc2
    o_c2 = kfm * po

    out = np.empty(10, dtype=np.float64)
    out[0] = -f0a - f0b + r_a + r_b
    out[1] = f0a - r_a - f_ab + r_2a
    out[2] = f0b - r_b - f_ba + r_2b
    out[3] = f_ab + f_ba - r_2a - r_2b - f_d1 - f_d2 - f_c1 - f_c2 + r_d1 + r_d2 + r_c1 + r_c2
    out[4] = f_d1 - r_d1
    out[5] = f_d2 - r_d2
    out[6] = f_c1 - r_c1 - c1_o + o_c1
    out[7] = f_c2 - r_c2 - c2_o + o_c2
    out[8] = c1_o + c2_o - o_c1 - o_c2
    out[9] = 0.0
    return out


@njit(cache=True)
def _calib_one(rates: np.ndarray, pulse_amp: float, tau_ms: float, dt_ms: float) -> np.ndarray:
    duration_ms = 600.0
    steps = int(round(duration_ms / dt_ms))
    state = np.zeros(10, dtype=np.float64)
    state[0] = 1.0
    g = 0.0
    decay = math.exp(-dt_ms / tau_ms)
    next_pulse = 0
    peak = 0.0
    h_sum = 0.0
    h_n = 0
    e_sum = 0.0
    e_n = 0
    l_sum = 0.0
    l_n = 0
    j_full = 0.0
    mass_error = 0.0
    minimum_state = 0.0

    for step in range(steps + 1):
        t = step * dt_ms
        while next_pulse < PULSE_COUNT and t + 0.5 * dt_ms >= next_pulse * PULSE_INTERVAL_MS:
            g += pulse_amp
            next_pulse += 1
        g_eff = g / (1.0 + g)
        po = state[8]
        if po > peak:
            peak = po
        j_full += po * dt_ms
        if 140.0 <= t < 250.0:
            h_sum += po
            h_n += 1
        if 140.0 <= t < 200.0:
            e_sum += po
            e_n += 1
        if 200.0 <= t < 500.0:
            l_sum += po
            l_n += 1

        mass = 0.0
        for k in range(10):
            mass += state[k]
            if state[k] < minimum_state:
                minimum_state = state[k]
        err = abs(mass - 1.0)
        if err > mass_error:
            mass_error = err
        if step == steps:
            break
        d1 = _calib_rhs(state, g_eff, rates)
        midpoint = state + 0.5 * dt_ms * d1
        g_next = g * decay
        g_next_eff = g_next / (1.0 + g_next)
        dmid = _calib_rhs(midpoint, 0.5 * (g_eff + g_next_eff), rates)
        state = state + dt_ms * dmid
        g = g_next

    out = np.empty(9, dtype=np.float64)
    h = h_sum / max(1, h_n)
    e = e_sum / max(1, e_n)
    l = l_sum / max(1, l_n)
    out[0] = peak
    out[1] = h
    out[2] = e
    out[3] = l
    out[4] = j_full
    out[5] = 60.0 * e + 300.0 * l
    out[6] = mass_error
    out[7] = minimum_state
    out[8] = 1.0 if mass_error < 1e-6 and minimum_state > -1e-7 else 0.0
    return out


@njit(parallel=True, cache=True)
def calibration_ensemble(rates: np.ndarray, pulse_amp: np.ndarray, tau_ms: np.ndarray, dt_ms: float) -> np.ndarray:
    n = rates.shape[0]
    out = np.empty((n, 9), dtype=np.float64)
    for i in prange(n):
        out[i, :] = _calib_one(rates[i], pulse_amp[i], tau_ms[i], dt_ms)
    return out


# -----------------------------------------------------------------------------
# Final Base response model: 8 explicit occupancies, P2 by conservation, RK4.
# -----------------------------------------------------------------------------

@njit(cache=True)
def _base_rhs(y: np.ndarray, g: float, rates: np.ndarray, lambda_b: float) -> np.ndarray:
    # Explicit state order: P0, P1A, P1B, PD1, PD2, PC1, PC2, PO.
    p0, p1a, p1b, pd1, pd2, pc1, pc2, po = y
    p2 = 1.0 - (p0 + p1a + p1b + pd1 + pd2 + pc1 + pc2 + po)
    kon_a, kon_b, koff_a, koff_b = rates[0], rates[1], rates[2], rates[3]
    kfp, ksp, kfm, ksm = rates[4], rates[5], rates[6], rates[7]
    kd1p, kd1m, kd2p, kd2m = rates[8], rates[9], rates[10], rates[11]

    out = np.empty(8, dtype=np.float64)
    out[0] = -g * kon_a * p0 - lambda_b * g * kon_b * p0 + koff_a * p1a + lambda_b * koff_b * p1b
    out[1] = g * kon_a * p0 - koff_a * p1a - g * kon_b * p1a + koff_b * p2
    out[2] = lambda_b * (g * kon_b * p0 - koff_b * p1b - g * kon_a * p1b + koff_a * p2)
    out[3] = kd1p * p2 - kd1m * pd1
    out[4] = kd2p * p2 - kd2m * pd2
    out[5] = kfp * p2 - kfm * pc1 - ksp * pc1 + ksm * po
    out[6] = ksp * p2 - ksm * pc2 - kfp * pc2 + kfm * po
    out[7] = ksp * pc1 + kfp * pc2 - (ksm + kfm) * po
    return out


@njit(cache=True)
def _median_prefix(values: np.ndarray, n: int) -> float:
    tmp = np.empty(n, dtype=np.float64)
    for i in range(n):
        tmp[i] = values[i]
    tmp.sort()
    if n % 2 == 1:
        return tmp[n // 2]
    return 0.5 * (tmp[n // 2 - 1] + tmp[n // 2])


@njit(cache=True)
def _base_plateau_one(rates: np.ndarray, g0: np.ndarray, gm: np.ndarray, g1: np.ndarray, dt_ms: float, lambda_b: float) -> tuple[float, float]:
    y = np.zeros(8, dtype=np.float64)
    y[0] = 1.0
    # Each retained interval is 3.5 ms wide, so 96 slots safely covers dt >= 0.04 ms.
    samples = np.empty((19, 96), dtype=np.float64)
    counts = np.zeros(19, dtype=np.int64)
    minimum_state = 0.0
    nsteps = len(g0)

    for step in range(nsteps):
        h = dt_ms
        k1 = _base_rhs(y, g0[step], rates, lambda_b)
        k2 = _base_rhs(y + 0.5 * h * k1, gm[step], rates, lambda_b)
        k3 = _base_rhs(y + 0.5 * h * k2, gm[step], rates, lambda_b)
        k4 = _base_rhs(y + h * k3, g1[step], rates, lambda_b)
        y = y + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        p2 = 1.0
        for q in range(8):
            p2 -= y[q]
            if y[q] < minimum_state:
                minimum_state = y[q]
        if p2 < minimum_state:
            minimum_state = p2

        t1 = (step + 1) * dt_ms
        interval = int(math.floor((t1 + 1e-10) / PULSE_INTERVAL_MS))
        if 5 <= interval <= 23:
            lo = interval * PULSE_INTERVAL_MS + 0.8
            hi = (interval + 1) * PULSE_INTERVAL_MS - 0.8
            if t1 >= lo - 1e-10 and t1 <= hi + 1e-10:
                idx = interval - 5
                c = counts[idx]
                if c < 96:
                    samples[idx, c] = y[7]
                    counts[idx] = c + 1

    interval_medians = np.empty(19, dtype=np.float64)
    for idx in range(19):
        if counts[idx] <= 0:
            return np.nan, minimum_state
        interval_medians[idx] = _median_prefix(samples[idx], counts[idx])
    plateau = _median_prefix(interval_medians, 19)
    return plateau, minimum_state


@njit(parallel=True, cache=True)
def base_ratio_ensemble(rates: np.ndarray, g0: np.ndarray, gm: np.ndarray, g1: np.ndarray, dt_ms: float) -> np.ndarray:
    n = rates.shape[0]
    out = np.empty((n, 5), dtype=np.float64)
    for i in prange(n):
        control, min_c = _base_plateau_one(rates[i], g0, gm, g1, dt_ms, 1.0)
        blocked, min_b = _base_plateau_one(rates[i], g0, gm, g1, dt_ms, 0.0)
        ratio = np.nan
        if np.isfinite(control) and np.isfinite(blocked) and control > 1e-14:
            ratio = blocked / control
        out[i, 0] = control
        out[i, 1] = blocked
        out[i, 2] = ratio
        out[i, 3] = min_c
        out[i, 4] = min_b
    return out


def forcing_arrays(g_peak: float, tau_ms: float, dt_ms: float) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    if not (0.0 < g_peak < 1.0):
        raise ValueError("G_peak must be in (0,1)")
    pulse_times = np.arange(PULSE_COUNT, dtype=float) * PULSE_INTERVAL_MS
    raw_peak = g_peak / (1.0 - g_peak)
    denom = np.sum(np.exp(-(TRAIN_END_MS - pulse_times) / tau_ms))
    amp = raw_peak / denom
    nsteps = int(round(TRAIN_END_MS / dt_ms))
    t0 = np.arange(nsteps, dtype=float) * dt_ms
    tm = t0 + 0.5 * dt_ms
    t1 = t0 + dt_ms

    def evaluate(times: np.ndarray) -> np.ndarray:
        raw = np.zeros_like(times)
        before = times < TRAIN_END_MS - 1e-12
        after = ~before
        for pulse in pulse_times:
            mask = before & (times >= pulse - 1e-12)
            raw[mask] += amp * np.exp(-(times[mask] - pulse) / tau_ms)
        raw[after] = raw_peak * np.exp(-(times[after] - TRAIN_END_MS) / tau_ms)
        return raw / (1.0 + raw)

    return evaluate(t0), evaluate(tm), evaluate(t1), float(amp)


def run_base_condition(rates: np.ndarray, g: float, tau: float, dt: float) -> tuple[np.ndarray, float]:
    g0, gm, g1, amp = forcing_arrays(g, tau, dt)
    return base_ratio_ensemble(rates, g0, gm, g1, dt), amp


def full_grid_screen(strong: pd.DataFrame, outdir: Path, dt: float, force: bool) -> pd.DataFrame:
    chk = outdir / "_checkpoints" / "full_grid"
    if force and chk.exists():
        shutil.rmtree(chk)
    chk.mkdir(parents=True, exist_ok=True)
    rates = strong[RATE_COLS].to_numpy(float)
    sample_ids = strong["sample_id"].to_numpy(int)
    rows: list[pd.DataFrame] = []
    total = len(G_PEAK_GRID) * len(TAU_GRID)
    done = 0
    for g in G_PEAK_GRID:
        for tau in TAU_GRID:
            cid = condition_id(float(g), float(tau))
            p = chk / f"{cid}.csv.gz"
            if not p.exists():
                metrics, forcing_amp = run_base_condition(rates, float(g), float(tau), dt)
                frame = pd.DataFrame({
                    "sample_id": sample_ids,
                    "condition_id": cid,
                    "G_peak": float(g),
                    "tau_G_ms": float(tau),
                    "forcing_pulse_increment": forcing_amp,
                    "control_plateau_PO": metrics[:, 0],
                    "blocked_plateau_PO": metrics[:, 1],
                    "r_plateau": metrics[:, 2],
                    "minimum_state_control": metrics[:, 3],
                    "minimum_state_blocked": metrics[:, 4],
                    "dt_ms": dt,
                })
                frame.to_csv(p, index=False, compression="gzip")
            rows.append(pd.read_csv(p))
            done += 1
            print(f"Step25 full grid: {done}/{total} {cid}", flush=True)
    return pd.concat(rows, ignore_index=True)


def confirm_max_nodes(strong: pd.DataFrame, maxima: pd.DataFrame, dt: float) -> pd.DataFrame:
    outputs = []
    for (g, tau), local in maxima.groupby(["G_peak", "tau_G_ms"], sort=True):
        ids = local["sample_id"].astype(int).tolist()
        sub = strong.set_index("sample_id").loc[ids].reset_index()
        metrics, amp = run_base_condition(sub[RATE_COLS].to_numpy(float), float(g), float(tau), dt)
        frame = pd.DataFrame({
            "sample_id": sub["sample_id"].to_numpy(int),
            "G_peak": float(g),
            "tau_G_ms": float(tau),
            "forcing_pulse_increment": amp,
            "control_plateau_PO_dt0p05": metrics[:, 0],
            "blocked_plateau_PO_dt0p05": metrics[:, 1],
            "r_plateau_dt0p05": metrics[:, 2],
            "minimum_state_control_dt0p05": metrics[:, 3],
            "minimum_state_blocked_dt0p05": metrics[:, 4],
        })
        outputs.append(frame)
    return pd.concat(outputs, ignore_index=True)


def build_geometry(accepted_broad: pd.DataFrame, strong_ids: set[int], outdir: Path, seed: int) -> dict[str, Any]:
    logs = np.log(accepted_broad[RATE_COLS].to_numpy(float))
    mean = logs.mean(axis=0)
    sd = logs.std(axis=0, ddof=1)
    if np.any(sd <= 0):
        raise RuntimeError("Zero variance in accepted broad log-rate coordinates")
    z = (logs - mean) / sd
    u, s, vt = np.linalg.svd(z, full_matrices=False)
    eigen = (s * s) / (len(z) - 1)
    explained = eigen / eigen.sum()
    coords = z @ vt[:6].T
    geo = accepted_broad[["sample_id", "calibration_score_corrected"]].copy()
    geo["is_strong_step24"] = geo["sample_id"].astype(int).isin(strong_ids)
    for i in range(6):
        geo[f"PC{i+1}"] = coords[:, i]
    geo.to_csv(outdir / "08_geometry_coordinates.csv.gz", index=False, compression="gzip")

    strong_mask = geo["is_strong_step24"].to_numpy(bool)
    strong_geo = geo.loc[strong_mask].reset_index(drop=True)
    strong_z = z[strong_mask]
    d = cdist(strong_z, strong_z)
    pair_rows = []
    for i in range(len(strong_geo)):
        for j in range(i + 1, len(strong_geo)):
            pair_rows.append({
                "sample_id_a": int(strong_geo.loc[i, "sample_id"]),
                "sample_id_b": int(strong_geo.loc[j, "sample_id"]),
                "standardized_log_rate_distance": float(d[i, j]),
            })
    pd.DataFrame(pair_rows).to_csv(outdir / "09_geometry_pairwise_distances.csv", index=False)

    d_nn = d.copy()
    np.fill_diagonal(d_nn, np.inf)
    nearest_idx = np.argmin(d_nn, axis=1)
    nn = pd.DataFrame({
        "sample_id": strong_geo["sample_id"].astype(int),
        "nearest_sample_id": strong_geo.loc[nearest_idx, "sample_id"].astype(int).to_numpy(),
        "nearest_distance": d_nn[np.arange(len(d_nn)), nearest_idx],
    })
    nn.to_csv(outdir / "10_geometry_nearest_neighbors.csv", index=False)

    # OPTICS sensitivity: no prespecified number of clusters.
    assignments: dict[str, np.ndarray] = {}
    sensitivity_rows = []
    n = len(strong_z)
    for ms in (3, 4, 5):
        for xi in (0.03, 0.05, 0.10):
            min_cluster_size = max(ms, int(math.ceil(0.15 * n)))
            model = OPTICS(min_samples=ms, xi=xi, min_cluster_size=min_cluster_size, metric="euclidean", cluster_method="xi")
            labels = model.fit_predict(strong_z)
            key = f"ms{ms}_xi{xi:.2f}"
            assignments[key] = labels.copy()
            clusters = sorted(set(labels) - {-1})
            sensitivity_rows.append({
                "configuration": key,
                "min_samples": ms,
                "xi": xi,
                "min_cluster_size": min_cluster_size,
                "cluster_count": len(clusters),
                "noise_fraction": float(np.mean(labels < 0)),
            })
            if ms == 4 and math.isclose(xi, 0.05):
                canonical_model = model
                canonical_labels = labels.copy()
    sens = pd.DataFrame(sensitivity_rows)
    sens.to_csv(outdir / "11_optics_sensitivity.csv", index=False)
    canonical = strong_geo[["sample_id"]].copy()
    canonical["optics_label"] = canonical_labels
    canonical["optics_reachability"] = canonical_model.reachability_
    canonical["optics_core_distance"] = canonical_model.core_distances_
    order_rank = np.empty(n, dtype=int)
    order_rank[canonical_model.ordering_] = np.arange(n)
    canonical["optics_order_rank"] = order_rank
    canonical.to_csv(outdir / "11A_optics_assignments.csv", index=False)

    keys = list(assignments)
    ari_values = []
    ari_rows = []
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            ari = float(adjusted_rand_score(assignments[keys[i]], assignments[keys[j]]))
            ari_values.append(ari)
            ari_rows.append({"config_a": keys[i], "config_b": keys[j], "adjusted_rand_index": ari})
    pd.DataFrame(ari_rows).to_csv(outdir / "11B_optics_pairwise_ari.csv", index=False)

    # MST over the 40 strong points.
    mst = minimum_spanning_tree(d).tocoo()
    mst_rows = []
    for a, b, val in zip(mst.row, mst.col, mst.data):
        mst_rows.append({
            "sample_id_a": int(strong_geo.loc[a, "sample_id"]),
            "sample_id_b": int(strong_geo.loc[b, "sample_id"]),
            "distance": float(val),
        })
    mst_df = pd.DataFrame(mst_rows)
    mst_df.to_csv(outdir / "11C_strong_mst_edges.csv", index=False)

    # Compare within-strong nearest-neighbour compactness with random 40-point subsets.
    rng = np.random.default_rng(seed)
    strong_median_nn = float(np.median(nn["nearest_distance"]))
    random_medians = np.empty(500, dtype=float)
    for rep in range(len(random_medians)):
        idx = rng.choice(len(z), size=n, replace=False)
        local_d = cdist(z[idx], z[idx])
        np.fill_diagonal(local_d, np.inf)
        random_medians[rep] = np.median(np.min(local_d, axis=1))
    p_compact = float((1 + np.sum(random_medians <= strong_median_nn)) / (len(random_medians) + 1))
    pd.DataFrame({"random_subset_median_nearest_distance": random_medians}).to_csv(
        outdir / "11D_random_subset_compactness.csv", index=False
    )

    return {
        "log_rate_mean": dict(zip(RATE_COLS, map(float, mean))),
        "log_rate_sd": dict(zip(RATE_COLS, map(float, sd))),
        "pca_explained_variance_fraction": [float(x) for x in explained[:6]],
        "canonical_optics_cluster_count": int(len(set(canonical_labels) - {-1})),
        "canonical_optics_noise_fraction": float(np.mean(canonical_labels < 0)),
        "optics_pairwise_ari_median": float(np.median(ari_values)) if ari_values else float("nan"),
        "strong_median_nearest_distance": strong_median_nn,
        "random_subset_median_nearest_distance_median": float(np.median(random_medians)),
        "compactness_one_sided_empirical_p": p_compact,
        "mst_edge_median": float(mst_df["distance"].median()) if len(mst_df) else float("nan"),
        "mst_edge_max": float(mst_df["distance"].max()) if len(mst_df) else float("nan"),
        "canonical_labels": {int(sid): int(lbl) for sid, lbl in zip(strong_geo["sample_id"], canonical_labels)},
        "strong_z": strong_z,
        "strong_ids_order": strong_geo["sample_id"].astype(int).to_numpy(),
        "pairwise_distance_matrix": d,
    }


def select_representatives(confirmed: pd.DataFrame, geometry: dict[str, Any], max_centres: int) -> pd.DataFrame:
    candidates: list[tuple[int, str]] = []
    c = confirmed.copy()
    c = c[c["strong_confirmed_dt0p05"]].copy()
    if c.empty:
        raise RuntimeError("No numerically confirmed strong candidates")
    candidates.append((int(c.loc[c["calibration_score_corrected"].idxmin(), "sample_id"]), "best_calibrated_strong"))
    exp_typical = float(np.median(EXPERIMENTAL_RATIOS))
    idx = np.argmin(np.abs(np.log(c["r_plateau_dt0p05"].to_numpy(float) / exp_typical)))
    candidates.append((int(c.iloc[idx]["sample_id"]), "experimental_typical_match"))
    idx = np.argmin(np.abs(np.log(c["r_plateau_dt0p05"].to_numpy(float) / EXPERIMENTAL_RATIOS[-1])))
    candidates.append((int(c.iloc[idx]["sample_id"]), "extreme_experimental_match"))
    candidates.append((int(c.loc[c["r_plateau_dt0p05"].idxmax(), "sample_id"]), "largest_confirmed_response"))
    idx = np.argmin(np.abs(np.log(c["r_plateau_dt0p05"].to_numpy(float) / EXPERIMENTAL_RATIOS[0])))
    candidates.append((int(c.iloc[idx]["sample_id"]), "minimum_experimental_match"))

    labels = geometry["canonical_labels"]
    ids = geometry["strong_ids_order"]
    d = geometry["pairwise_distance_matrix"]
    for label in sorted(set(labels.values()) - {-1}):
        member_positions = [i for i, sid in enumerate(ids) if labels[int(sid)] == label]
        if not member_positions:
            continue
        local = d[np.ix_(member_positions, member_positions)]
        medoid_pos = member_positions[int(np.argmin(local.sum(axis=1)))]
        candidates.append((int(ids[medoid_pos]), f"optics_cluster_{label}_medoid"))

    roles: dict[int, list[str]] = {}
    for sid, role in candidates:
        roles.setdefault(sid, []).append(role)
    ordered = list(roles.items())[:max_centres]
    rows = []
    lookup = c.set_index("sample_id")
    for sid, role_list in ordered:
        if sid not in lookup.index:
            continue
        r = lookup.loc[sid]
        rows.append({
            "sample_id": sid,
            "roles": ";".join(role_list),
            "optics_label": labels.get(sid, -1),
            "calibration_score_corrected": float(r["calibration_score_corrected"]),
            "G_peak": float(r["G_peak"]),
            "tau_G_ms": float(r["tau_G_ms"]),
            "r_plateau_dt0p05": float(r["r_plateau_dt0p05"]),
            "control_plateau_PO_dt0p05": float(r["control_plateau_PO_dt0p05"]),
        })
    return pd.DataFrame(rows)


def local_perturbations(
    reps: pd.DataFrame,
    strong_full: pd.DataFrame,
    targets: dict[str, Any],
    cutoff_mid: float,
    cutoff_strict: float,
    strong_threshold: float,
    accepted_broad: pd.DataFrame,
    full_broad: pd.DataFrame,
    outdir: Path,
    n_per: int,
    seed: int,
    calibration_dt: float,
    response_dt: float,
    force: bool,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    chk = outdir / "_checkpoints" / "local"
    if force and chk.exists():
        shutil.rmtree(chk)
    chk.mkdir(parents=True, exist_ok=True)
    center_lookup = strong_full.set_index("sample_id")
    accepted_min = accepted_broad[RATE_COLS].min().to_numpy(float)
    accepted_max = accepted_broad[RATE_COLS].max().to_numpy(float)
    full_min = full_broad[RATE_COLS].min().to_numpy(float)
    full_max = full_broad[RATE_COLS].max().to_numpy(float)
    files = []
    summary_rows = []

    for _, rep in reps.iterrows():
        sid = int(rep["sample_id"])
        center = center_lookup.loc[sid]
        base_rates = center[RATE_COLS].to_numpy(float)
        pulse_amp = float(center["pulse_amplitude_au"])
        calib_tau = float(center["glutamate_tau_ms"])
        g = float(rep["G_peak"])
        tau = float(rep["tau_G_ms"])
        for frac in (0.05, 0.10):
            tag = f"sample{sid}_pm{int(round(frac*100)):02d}"
            path = chk / f"{tag}.csv.gz"
            if not path.exists():
                rng = np.random.default_rng(seed + sid * 101 + int(round(frac * 10000)))
                multipliers = rng.uniform(1.0 - frac, 1.0 + frac, size=(n_per, len(RATE_COLS)))
                rates = base_rates[None, :] * multipliers
                response, _ = run_base_condition(rates, g, tau, response_dt)
                calibration = calibration_ensemble(
                    rates,
                    np.full(n_per, pulse_amp, dtype=float),
                    np.full(n_per, calib_tau, dtype=float),
                    calibration_dt,
                )
                score = corrected_score_from_metrics(calibration[:, 1], calibration[:, 2], calibration[:, 3], calibration[:, 5], targets)
                valid_calib = calibration[:, 8] > 0.5
                compatible_mid = valid_calib & (score <= cutoff_mid)
                compatible_strict = valid_calib & (score <= cutoff_strict)
                strong = np.isfinite(response[:, 2]) & (response[:, 2] >= strong_threshold)
                within_acc = np.all((rates >= accepted_min) & (rates <= accepted_max), axis=1)
                within_full = np.all((rates >= full_min) & (rates <= full_max), axis=1)
                frame = pd.DataFrame({
                    "center_sample_id": sid,
                    "roles": str(rep["roles"]),
                    "optics_label": int(rep["optics_label"]),
                    "perturbation_fraction": frac,
                    "point_id": np.arange(n_per, dtype=int),
                    "G_peak": g,
                    "tau_G_ms": tau,
                    "control_plateau_PO": response[:, 0],
                    "blocked_plateau_PO": response[:, 1],
                    "r_plateau": response[:, 2],
                    "strong": strong,
                    "calibration_H": calibration[:, 1],
                    "calibration_E": calibration[:, 2],
                    "calibration_L": calibration[:, 3],
                    "calibration_J_full_0_600_ms": calibration[:, 4],
                    "calibration_J_corrected_140_500_ms": calibration[:, 5],
                    "calibration_score_corrected": score,
                    "calibration_valid": valid_calib,
                    "compatible_midpoint": compatible_mid,
                    "compatible_strict_max_accepted": compatible_strict,
                    "strong_and_compatible": strong & compatible_mid,
                    "within_corrected_accepted_rate_bounds": within_acc,
                    "within_full_broad_sample_rate_bounds": within_full,
                    "minimum_state_response_control": response[:, 3],
                    "minimum_state_response_blocked": response[:, 4],
                    "minimum_state_calibration": calibration[:, 7],
                })
                for k, col in enumerate(RATE_COLS):
                    frame[col] = rates[:, k]
                frame.to_csv(path, index=False, compression="gzip")
            frame = pd.read_csv(path)
            files.append(frame)
            summary_rows.append({
                "center_sample_id": sid,
                "roles": str(rep["roles"]),
                "optics_label": int(rep["optics_label"]),
                "perturbation_fraction": frac,
                "n": len(frame),
                "strong_fraction": float(frame["strong"].mean()),
                "compatible_fraction": float(frame["compatible_midpoint"].mean()),
                "joint_fraction": float(frame["strong_and_compatible"].mean()),
                "joint_fraction_within_full_broad_bounds": (
                    float(frame.loc[frame["within_full_broad_sample_rate_bounds"], "strong_and_compatible"].mean())
                    if frame["within_full_broad_sample_rate_bounds"].any() else float("nan")
                ),
                "median_r_plateau": float(frame["r_plateau"].median()),
                "median_control_plateau_PO": float(frame["control_plateau_PO"].median()),
                "median_corrected_score": float(frame["calibration_score_corrected"].median()),
            })
            print(f"Step25 local: center={sid} perturb={frac:.0%} n={len(frame)}", flush=True)
    return pd.concat(files, ignore_index=True), pd.DataFrame(summary_rows)


def experiment_matches(confirmed: pd.DataFrame) -> pd.DataFrame:
    rows = []
    vals = confirmed["r_plateau_dt0p05"].to_numpy(float)
    for target in EXPERIMENTAL_RATIOS:
        idx = int(np.nanargmin(np.abs(np.log(vals / target))))
        row = confirmed.iloc[idx]
        rows.append({
            "experimental_ratio": float(target),
            "sample_id": int(row["sample_id"]),
            "model_ratio_dt0p05": float(row["r_plateau_dt0p05"]),
            "relative_error": float(abs(row["r_plateau_dt0p05"] / target - 1.0)),
            "absolute_log_ratio_distance": float(abs(math.log(row["r_plateau_dt0p05"] / target))),
            "G_peak": float(row["G_peak"]),
            "tau_G_ms": float(row["tau_G_ms"]),
            "calibration_score_corrected": float(row["calibration_score_corrected"]),
        })
    return pd.DataFrame(rows)


def make_figures(confirmed: pd.DataFrame, geometry_coords: pd.DataFrame, local_summary: pd.DataFrame, outdir: Path) -> None:
    # A: response ratio vs raw control denominator.
    fig, ax = plt.subplots(figsize=(6.4, 4.2), constrained_layout=True)
    ax.scatter(confirmed["control_plateau_PO_dt0p05"], confirmed["r_plateau_dt0p05"], s=28)
    ax.axhline(STRONG_DEFAULT, linestyle="--", linewidth=0.9)
    for sid in (10052, 27084, 37318, 1900):
        row = confirmed[confirmed["sample_id"] == sid]
        if not row.empty:
            r = row.iloc[0]
            ax.annotate(str(sid), (r["control_plateau_PO_dt0p05"], r["r_plateau_dt0p05"]), xytext=(4, 4), textcoords="offset points", fontsize=7)
    ax.set_xlabel("Control inter-pulse plateau, $P_O$")
    ax.set_ylabel("Blocked/control plateau ratio")
    ax.set_yscale("log")
    ax.set_title("A  Strong candidates: magnitude vs denominator")
    fig.savefig(outdir / "Fig25A_response_denominator.pdf")
    fig.savefig(outdir / "Fig25A_response_denominator.png", dpi=300)
    plt.close(fig)

    # B: PCA geometry.
    fig, ax = plt.subplots(figsize=(6.4, 4.4), constrained_layout=True)
    bg = geometry_coords[~geometry_coords["is_strong_step24"]]
    st = geometry_coords[geometry_coords["is_strong_step24"]]
    ax.scatter(bg["PC1"], bg["PC2"], s=5, alpha=0.12, label="corrected accepted broad")
    ax.scatter(st["PC1"], st["PC2"], s=30, label="strong")
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_title("B  Corrected accepted rate geometry")
    ax.legend(frameon=False, fontsize=8)
    fig.savefig(outdir / "Fig25B_pca_geometry.pdf")
    fig.savefig(outdir / "Fig25B_pca_geometry.png", dpi=300)
    plt.close(fig)

    # C: local robustness by representative centre.
    if not local_summary.empty:
        fig, ax = plt.subplots(figsize=(7.2, 4.6), constrained_layout=True)
        reps = list(dict.fromkeys(local_summary["center_sample_id"].astype(int).tolist()))
        x = np.arange(len(reps), dtype=float)
        width = 0.34
        for offset, frac in ((-width / 2, 0.05), (width / 2, 0.10)):
            local = local_summary[np.isclose(local_summary["perturbation_fraction"], frac)].set_index("center_sample_id")
            vals = [local.loc[sid, "joint_fraction"] if sid in local.index else np.nan for sid in reps]
            ax.bar(x + offset, vals, width=width, label=f"±{int(frac*100)}%")
        ax.set_xticks(x, [str(v) for v in reps], rotation=45, ha="right")
        ax.set_ylim(0, 1)
        ax.set_ylabel("Strong + corrected-compatible fraction")
        ax.set_xlabel("Representative centre sample ID")
        ax.set_title("C  Local joint robustness across strong centres")
        ax.legend(frameon=False)
        fig.savefig(outdir / "Fig25C_local_robustness.pdf")
        fig.savefig(outdir / "Fig25C_local_robustness.png", dpi=300)
        plt.close(fig)

    # D: calibration score vs response ratio.
    fig, ax = plt.subplots(figsize=(6.4, 4.2), constrained_layout=True)
    ax.scatter(confirmed["calibration_score_corrected"], confirmed["r_plateau_dt0p05"], s=28)
    ax.set_xlabel("Corrected control-calibration score")
    ax.set_ylabel("Confirmed maximum response ratio")
    ax.set_yscale("log")
    ax.set_title("D  Calibration quality and amplification")
    fig.savefig(outdir / "Fig25D_score_response.pdf")
    fig.savefig(outdir / "Fig25D_score_response.png", dpi=300)
    plt.close(fig)


def self_test() -> None:
    g0, gm, g1, amp = forcing_arrays(0.95, 10.0, 0.1)
    if abs(g1[-1] - 0.95) > 1e-12:
        raise AssertionError(f"Forcing endpoint mismatch: {g1[-1]}")
    rates = np.array([[0.8, 0.35, 0.30, 0.10, 2.5, 0.45, 0.45, 0.45, 0.25, 0.09, 0.55, 0.01]], dtype=float)
    out = base_ratio_ensemble(rates, g0, gm, g1, 0.1)
    if not np.isfinite(out[0, 2]) or out[0, 0] <= 0:
        raise AssertionError("Base response self-test failed")
    cal = calibration_ensemble(rates, np.array([0.7]), np.array([10.0]), 0.05)
    if cal[0, 8] < 0.5 or not np.all(np.isfinite(cal[0, :6])):
        raise AssertionError("Historical calibration replay self-test failed")
    print("SELF-TEST PASSED")


def main() -> None:
    args = parse_args()
    if args.self_test:
        if args.threads > 0:
            set_num_threads(args.threads)
        self_test()
        return
    if args.local_n < 100:
        raise ValueError("--local-n must be at least 100")
    if args.grid_dt_ms <= 0 or args.confirm_dt_ms <= 0 or args.calibration_dt_ms <= 0:
        raise ValueError("All dt values must be positive")
    if args.confirm_dt_ms > args.grid_dt_ms:
        raise ValueError("Confirmation dt must be <= grid dt")
    if args.threads > 0:
        set_num_threads(args.threads)

    step24 = discover_step24(args.step24)
    s24 = load_step24(step24)
    targets = recover_corrected_targets(s24["scores"])
    broad_shift = s24["shift"].loc[s24["shift"]["prior"] == "broad"].iloc[0]
    cutoff_mid = float(broad_shift["corrected_cutoff_midpoint"])
    cutoff_strict = float(broad_shift["corrected_max_accepted_score"])
    threshold = float(s24["threshold"])

    preflight = {
        "step24": str(step24),
        "step24_status": s24["decision"].get("status"),
        "corrected_strong_n": int(len(s24["strong"])),
        "strong_threshold": threshold,
        "corrected_broad_cutoff_midpoint": cutoff_mid,
        "corrected_broad_max_accepted_score": cutoff_strict,
        "recovered_targets": targets,
    }
    if args.preflight_only:
        print(json.dumps(preflight, indent=2))
        return

    step10_path, preferred_member = discover_step10(args.step10, step24)
    step10, step10_source = read_step10(step10_path, preferred_member)
    if set(step10["prior"].unique()) != {"broad", "reference"}:
        raise RuntimeError(f"Expected broad+reference in Step10, found {sorted(step10['prior'].unique())}")
    if len(step10[step10["prior"] == "broad"]) != 50000 or len(step10[step10["prior"] == "reference"]) != 50000:
        raise RuntimeError("Expected 50,000 Step10 candidates per prior")

    outdir = args.output.expanduser().resolve()
    if args.force and outdir.exists():
        shutil.rmtree(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    input_audit = {
        **preflight,
        "script_version": SCRIPT_VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "step10_source": step10_source,
        "step10_archive_sha256": sha256_file(step10_path) if step10_path.is_file() else None,
        "step10_rows": int(len(step10)),
        "step10_rows_by_prior": step10.groupby("prior").size().astype(int).to_dict(),
        "frozen_forcing_grid_G_peak": G_PEAK_GRID.tolist(),
        "frozen_forcing_grid_tau_G_ms": TAU_GRID.tolist(),
        "pulse_count": PULSE_COUNT,
        "pulse_interval_ms": PULSE_INTERVAL_MS,
        "train_end_ms": TRAIN_END_MS,
        "response_endpoint": "nested median of 19 retained inter-pulse PO windows after discarding first five intervals",
        "calibration_integral": "J_140_500_ms = 60*E_140_200 + 300*L_200_500",
        "no_new_denominator_cutoff": True,
        "threads": get_num_threads(),
    }
    json_dump(outdir / "00_input_audit.json", input_audit)
    json_dump(outdir / "01_corrected_target_reconstruction.json", targets)

    scores = s24["scores"].copy()
    # Step24 already carries the inherited calibration forcing columns. Merge only
    # the 12 kinetic rates from Step10 to avoid duplicate pulse/tau suffixes.
    merge_cols = ["prior", "sample_id", *RATE_COLS]
    accepted = scores[scores["accepted_corrected"]].merge(step10[merge_cols], on=["prior", "sample_id"], how="left", validate="one_to_one")
    if accepted[RATE_COLS].isna().any().any():
        raise RuntimeError("Failed to merge corrected accepted cohorts to Step10 rates")
    accepted.to_csv(outdir / "02_corrected_accepted_cohorts.csv.gz", index=False, compression="gzip")
    accepted_broad = accepted[accepted["prior"] == "broad"].copy().reset_index(drop=True)
    full_broad = step10[step10["prior"] == "broad"].copy().reset_index(drop=True)

    strong = s24["strong"].merge(
        accepted_broad[["sample_id", "calibration_score_corrected", "corrected_rank", "control_plateau_PO_140_250", "control_early_PO_140_200", "control_late_PO_200_500", "control_charge_PO_ms", "control_integrated_open_state_140_500_ms", *RATE_COLS, *AUX_COLS]],
        on=["sample_id", "calibration_score_corrected"], how="left", validate="one_to_one"
    )
    if strong[RATE_COLS].isna().any().any():
        raise RuntimeError("Failed to merge Step24 strong candidates to rate vectors")
    strong.to_csv(outdir / "03_strong_candidates_step24.csv", index=False)

    # Exact source-style calibration replay on the 40 original vectors.
    cal = calibration_ensemble(
        strong[RATE_COLS].to_numpy(float),
        strong["pulse_amplitude_au"].to_numpy(float),
        strong["glutamate_tau_ms"].to_numpy(float),
        args.calibration_dt_ms,
    )
    cal_score = corrected_score_from_metrics(cal[:, 1], cal[:, 2], cal[:, 3], cal[:, 5], targets)
    calib_audit = strong[["sample_id", "calibration_score_corrected"]].copy()
    calib_audit["score_replayed_from_source_integrator"] = cal_score
    calib_audit["abs_score_difference"] = np.abs(cal_score - calib_audit["calibration_score_corrected"].to_numpy(float))
    calib_audit["H_replayed"] = cal[:, 1]
    calib_audit["E_replayed"] = cal[:, 2]
    calib_audit["L_replayed"] = cal[:, 3]
    calib_audit["J_full_0_600_replayed"] = cal[:, 4]
    calib_audit["J_corrected_140_500_replayed"] = cal[:, 5]
    calib_audit["valid_replayed"] = cal[:, 8] > 0.5
    calib_audit.to_csv(outdir / "03A_strong_calibration_replay_audit.csv", index=False)
    max_calib_diff = float(calib_audit["abs_score_difference"].max())
    if max_calib_diff > 1e-6:
        raise RuntimeError(f"Source-style corrected calibration replay mismatch: max score difference {max_calib_diff:g}")

    grid = full_grid_screen(strong, outdir, args.grid_dt_ms, args.force)
    grid.to_csv(outdir / "04_full_80_node_grid.csv.gz", index=False, compression="gzip")
    max_idx = grid.groupby("sample_id")["r_plateau"].idxmax()
    maxima = grid.loc[max_idx].copy().sort_values("sample_id").reset_index(drop=True)
    maxima = maxima.rename(columns={
        "control_plateau_PO": "control_plateau_PO_dt0p1",
        "blocked_plateau_PO": "blocked_plateau_PO_dt0p1",
        "r_plateau": "r_plateau_dt0p1",
        "minimum_state_control": "minimum_state_control_dt0p1",
        "minimum_state_blocked": "minimum_state_blocked_dt0p1",
    })

    prelim = strong[["sample_id", "r_max_confirmed_or_screen", "G_peak_at_max", "tau_G_ms_at_max"]].merge(
        maxima[["sample_id", "G_peak", "tau_G_ms", "r_plateau_dt0p1"]], on="sample_id", how="left", validate="one_to_one"
    )
    prelim["r_abs_difference"] = np.abs(prelim["r_plateau_dt0p1"] - prelim["r_max_confirmed_or_screen"])
    prelim["r_rel_difference"] = prelim["r_abs_difference"] / np.maximum(np.abs(prelim["r_max_confirmed_or_screen"]), 1e-12)
    prelim["max_node_same"] = np.isclose(prelim["G_peak"], prelim["G_peak_at_max"]) & np.isclose(prelim["tau_G_ms"], prelim["tau_G_ms_at_max"])
    prelim["replay_pass"] = (prelim["r_abs_difference"] <= args.screen_replay_abs_tol) | (prelim["r_rel_difference"] <= args.screen_replay_rel_tol)
    prelim.to_csv(outdir / "04A_step24_response_replay_audit.csv", index=False)
    replay_fraction = float(prelim["replay_pass"].mean())
    witness_ok = bool(prelim.loc[prelim["sample_id"] == 37318, "replay_pass"].all())
    if (replay_fraction < 0.95 or not witness_ok) and not args.allow_replay_mismatch:
        json_dump(outdir / "04B_REPLAY_GATE_FAILED.json", {
            "replay_pass_fraction": replay_fraction,
            "witness_37318_pass": witness_ok,
            "instruction": "Do not interpret downstream Step25 outputs; inspect final Base replay implementation or rerun with explicit override only after audit.",
        })
        raise RuntimeError(f"Final Base response replay gate failed: pass fraction={replay_fraction:.3f}, witness_pass={witness_ok}")

    confirm = confirm_max_nodes(strong, maxima[["sample_id", "G_peak", "tau_G_ms"]], args.confirm_dt_ms)
    confirmed = maxima.merge(confirm, on=["sample_id", "G_peak", "tau_G_ms"], how="left", validate="one_to_one")
    confirmed = confirmed.merge(
        strong[["sample_id", "calibration_score_corrected", "corrected_rank", "control_plateau_PO_140_250", "control_early_PO_140_200", "control_late_PO_200_500", "control_charge_PO_ms", "control_integrated_open_state_140_500_ms", *RATE_COLS, *AUX_COLS]],
        on="sample_id", how="left", validate="one_to_one"
    )
    confirmed["r_abs_dt_difference"] = np.abs(confirmed["r_plateau_dt0p1"] - confirmed["r_plateau_dt0p05"])
    confirmed["r_rel_dt_difference"] = confirmed["r_abs_dt_difference"] / np.maximum(np.abs(confirmed["r_plateau_dt0p05"]), 1e-12)
    confirmed["strong_dt0p1"] = confirmed["r_plateau_dt0p1"] >= threshold
    confirmed["strong_confirmed_dt0p05"] = confirmed["r_plateau_dt0p05"] >= threshold
    confirmed["strong_classification_same"] = confirmed["strong_dt0p1"] == confirmed["strong_confirmed_dt0p05"]
    confirmed = confirmed.sort_values("r_plateau_dt0p05", ascending=False).reset_index(drop=True)
    confirmed.to_csv(outdir / "05_strong_candidates_confirmed.csv", index=False)

    matches = experiment_matches(confirmed)
    matches.to_csv(outdir / "06_experimental_ratio_matches.csv", index=False)

    denom = confirmed[["sample_id", "r_plateau_dt0p05", "control_plateau_PO_dt0p05", "blocked_plateau_PO_dt0p05", "calibration_score_corrected", "G_peak", "tau_G_ms"]].copy()
    denom["control_plateau_percentile_within_40"] = denom["control_plateau_PO_dt0p05"].rank(pct=True, method="average")
    denom["response_rank_within_40"] = denom["r_plateau_dt0p05"].rank(ascending=False, method="min").astype(int)
    denom.to_csv(outdir / "07_denominator_diagnostics.csv", index=False)

    geometry = build_geometry(accepted_broad, set(strong["sample_id"].astype(int)), outdir, args.seed)
    reps = select_representatives(confirmed, geometry, args.max_centres)
    reps.to_csv(outdir / "12_representative_centres.csv", index=False)

    local_points, local_summary = local_perturbations(
        reps,
        confirmed,
        targets,
        cutoff_mid,
        cutoff_strict,
        threshold,
        accepted_broad,
        full_broad,
        outdir,
        args.local_n,
        args.seed,
        args.calibration_dt_ms,
        args.grid_dt_ms,
        args.force,
    )
    local_points.to_csv(outdir / "13_local_perturbation_points.csv.gz", index=False, compression="gzip")
    local_summary.to_csv(outdir / "14_local_perturbation_summary.csv", index=False)

    # Decision: finite-sample geometry can indicate regions but cannot prove global connectivity.
    canonical_clusters = int(geometry["canonical_optics_cluster_count"])
    ari_med = float(geometry["optics_pairwise_ari_median"])
    robust_cluster_labels = set()
    for _, row in local_summary[np.isclose(local_summary["perturbation_fraction"], 0.05)].iterrows():
        if row["joint_fraction"] >= 0.25 and int(row["optics_label"]) >= 0:
            robust_cluster_labels.add(int(row["optics_label"]))
    if canonical_clusters >= 2 and len(robust_cluster_labels) >= 2 and ari_med >= 0.50:
        region_status = "MULTIPLE_ROBUST_KINETIC_REGIONS_INDICATED"
    elif canonical_clusters <= 1 and geometry["canonical_optics_noise_fraction"] < 0.50:
        region_status = "SINGLE_DOMINANT_SAMPLED_REGION_COMPATIBLE"
    else:
        region_status = "KINETIC_REGION_STRUCTURE_UNRESOLVED"

    numerical_pass = bool(confirmed["strong_classification_same"].all())
    decision = {
        "status": region_status,
        "interpretation_scope": "Finite corrected strong point-cloud geometry plus local perturbation support; not proof of global basin connectivity or physiological population structure.",
        "step24_strong_n": int(len(strong)),
        "confirmed_strong_n_dt0p05": int(confirmed["strong_confirmed_dt0p05"].sum()),
        "dt_strong_classification_agreement": float(confirmed["strong_classification_same"].mean()),
        "max_relative_dt_difference": float(confirmed["r_rel_dt_difference"].max()),
        "response_replay_pass_fraction": replay_fraction,
        "canonical_optics_cluster_count": canonical_clusters,
        "canonical_optics_noise_fraction": float(geometry["canonical_optics_noise_fraction"]),
        "optics_sensitivity_pairwise_ari_median": ari_med,
        "robust_optics_clusters_with_pm5_joint_fraction_ge_0p25": sorted(robust_cluster_labels),
        "strong_compactness_empirical_p": float(geometry["compactness_one_sided_empirical_p"]),
        "best_calibrated_strong_sample_id": int(confirmed.loc[confirmed["calibration_score_corrected"].idxmin(), "sample_id"]),
        "largest_confirmed_response_sample_id": int(confirmed.loc[confirmed["r_plateau_dt0p05"].idxmax(), "sample_id"]),
        "numerical_gate_pass": numerical_pass,
        "no_denominator_exclusion_applied": True,
    }
    json_dump(outdir / "15_scientific_decision.json", decision)

    geometry_coords = pd.read_csv(outdir / "08_geometry_coordinates.csv.gz")
    make_figures(confirmed, geometry_coords, local_summary, outdir)

    report_lines = [
        "# Step 25 corrected strong-ensemble validation and geometry",
        "",
        f"Step-24 corrected strong candidates: **{len(strong)}**.",
        f"Numerically confirmed at dt={args.confirm_dt_ms:g} ms: **{int(confirmed['strong_confirmed_dt0p05'].sum())}/{len(strong)}**.",
        f"Final Base replay gate pass fraction against Step-24 stored maxima: **{replay_fraction:.3f}**.",
        f"Canonical OPTICS cluster count: **{canonical_clusters}**; noise fraction: **{geometry['canonical_optics_noise_fraction']:.3f}**.",
        f"Decision: **{region_status}**.",
        "",
        "No new control-denominator cutoff was used. Raw control plateaus are stored in `07_denominator_diagnostics.csv`.",
        "Local perturbations modify all 12 kinetic rates simultaneously while retaining each centre's inherited historical calibration forcing.",
        "",
        "The one-versus-several-region conclusion is intentionally limited to sampled geometry and local perturbation evidence; it is not a global topological proof.",
    ]
    (outdir / "README_RESULTS.md").write_text("\n".join(report_lines) + "\n", encoding="utf-8")

    run_summary = {
        "pipeline_step": "25_corrected_strong_ensemble_validation_and_geometry",
        "script_version": SCRIPT_VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "numba_threads": get_num_threads(),
        "step24": str(step24),
        "step10_source": step10_source,
        "strong_threshold": threshold,
        "grid_dt_ms": args.grid_dt_ms,
        "confirmation_dt_ms": args.confirm_dt_ms,
        "calibration_dt_ms": args.calibration_dt_ms,
        "local_n_per_centre_per_level": args.local_n,
        "local_perturbation_levels": [0.05, 0.10],
        "representative_centres": reps.to_dict("records"),
        "experimental_ratios": EXPERIMENTAL_RATIOS.tolist(),
        "decision": decision,
    }
    json_dump(outdir / "run_summary.json", run_summary)

    if (outdir / "_checkpoints").exists():
        shutil.rmtree(outdir / "_checkpoints")
    print(f"Step25 results written to: {outdir}", flush=True)


if __name__ == "__main__":
    main()
