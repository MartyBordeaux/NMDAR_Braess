#!/usr/bin/env python3
"""NMDAR Braess Step 28 v1.0 - residual connectivity stress test.

Scientific purpose
------------------
Step 27 connected 36/40 corrected strong/control-compatible vectors, leaving
three residual blocks: {27498, 29959}, {6113}, and {37318}. Step 28 performs a
stronger finite-budget connectivity stress test using endpoint-preserving
multi-waypoint paths in the full 12-dimensional standardized log-rate space.

Primary path family (native_aux)
    - 12 kinetic rates are allowed to move through 2 or 3 freely optimized
      internal waypoints in full standardized log-rate space;
    - pulse_amplitude_au and glutamate_tau_ms follow their endpoint-preserving
      geometric interpolation.

Secondary path family (flex_aux)
    - the same kinetic path is used;
    - auxiliary calibration coordinates also receive freely optimized internal
      waypoints, with both endpoints preserved exactly.

All internal waypoint coordinates are constrained to the original broad-prior
box recovered from the full Step-10 candidate table. Search uses the 40 known
strong vectors and the 5000 corrected accepted broad vectors only as seeds;
these do not alter the frozen acceptance criteria.

A candidate bridge is accepted only after 61-point verification. At every path
point, corrected calibration is computed at dt=0.05 ms and the full frozen
80-node response grid is scanned at dt=0.1 ms, with the best node reconfirmed
at dt=0.05 ms and a full dt=0.05-ms fallback if needed. The strict corrected
cutoff is the primary acceptance criterion. No denominator exclusion is used.

Failure to connect residual components within this finite search budget is not
proof of global topological disconnection.
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
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from numba import get_num_threads, set_num_threads

import step25_engine_snapshot as eng


SCRIPT_VERSION = "1.0.0"
EXPECTED_STEP25_STATUS = "KINETIC_REGION_STRUCTURE_UNRESOLVED"
EXPECTED_STEP27_STATUS = "PARTIAL_CURVED_CONNECTIVITY_ONLY"
EXPECTED_STRICT_COMPONENTS = [
    frozenset({
        232, 1091, 1900, 3687, 6863, 7273, 8859, 10052, 12080, 13577,
        13715, 15395, 16950, 17407, 18319, 26689, 27084, 27495, 27872,
        28451, 28634, 29419, 30902, 31498, 32755, 34287, 36402, 38506,
        39172, 41698, 43651, 44207, 46580, 46829, 47327, 49029,
    }),
    frozenset({27498, 29959}),
    frozenset({6113}),
    frozenset({37318}),
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Step 28 residual strong-compatible connectivity stress test")
    p.add_argument("--step25", default="auto")
    p.add_argument("--step27", default="auto")
    p.add_argument(
        "--output",
        type=Path,
        default=Path("/root/nmda2/step_28/results_step_28_residual_connectivity_stress"),
    )
    p.add_argument("--threads", type=int, default=2)
    p.add_argument("--search-points", type=int, default=17)
    p.add_argument("--verify-points", type=int, default=61)
    p.add_argument("--population", type=int, default=32)
    p.add_argument("--generations", type=int, default=20)
    p.add_argument("--restarts", type=int, default=2)
    p.add_argument("--response-top", type=int, default=10)
    p.add_argument("--pair-attempts", type=int, default=6)
    p.add_argument("--deep-pair-attempts", type=int, default=2)
    p.add_argument("--max-component-pairs", type=int, default=6)
    p.add_argument("--max-waypoints", type=int, default=3, choices=[2, 3, 4])
    p.add_argument("--initial-sigma", type=float, default=0.9)
    p.add_argument("--verify-trigger", type=float, default=1.12)
    p.add_argument("--flex-trigger", type=float, default=1.50)
    p.add_argument("--seed", type=int, default=20260914)
    p.add_argument("--grid-dt-ms", type=float, default=0.1)
    p.add_argument("--confirm-dt-ms", type=float, default=0.05)
    p.add_argument("--calibration-dt-ms", type=float, default=0.05)
    p.add_argument("--disable-flex-aux", action="store_true")
    p.add_argument("--allow-accepted-bounds-fallback", action="store_true")
    p.add_argument("--force", action="store_true")
    p.add_argument("--resume", action="store_true")
    p.add_argument("--preflight-only", action="store_true")
    p.add_argument("--self-test", action="store_true")
    return p.parse_args()


def json_dump(path: Path, obj: Any) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def discover_step25(spec: str) -> Path:
    if spec != "auto":
        p = Path(spec).expanduser().resolve()
        if not p.is_dir():
            raise FileNotFoundError(p)
        return p
    candidates = [
        Path("/root/nmda2/step_25/results_step_25_corrected_strong_geometry"),
        Path("/root/nmda2/results_step_25_corrected_strong_geometry"),
        Path("/root/nmda2/NMDAR_Braess_step25_corrected_strong_geometry_v1_0/results_step_25_corrected_strong_geometry"),
    ]
    required = {"05_strong_candidates_confirmed.csv", "15_scientific_decision.json", "02_corrected_accepted_cohorts.csv.gz"}
    for p in candidates:
        if p.is_dir() and required.issubset({x.name for x in p.iterdir()}):
            return p.resolve()
    root = Path("/root/nmda2")
    if root.exists():
        for marker in root.rglob("15_scientific_decision.json"):
            q = marker.parent
            if required.issubset({x.name for x in q.iterdir()}):
                return q.resolve()
    raise FileNotFoundError("Could not auto-discover Step-25 results")


def discover_step27(spec: str) -> Path:
    if spec != "auto":
        p = Path(spec).expanduser().resolve()
        if not p.is_dir():
            raise FileNotFoundError(p)
        return p
    candidates = [
        Path("/root/nmda2/step_27/results_step_27_curved_bridge_optimization"),
        Path("/root/nmda2/results_step_27_curved_bridge_optimization"),
        Path("/root/nmda2/NMDAR_Braess_step27_curved_bridge_optimization_v1_0/results_step_27_curved_bridge_optimization"),
    ]
    required = {"04_bridge_attempts.csv", "07_verified_network_edges.csv", "08_final_components.csv", "10_scientific_decision.json"}
    for p in candidates:
        if p.is_dir() and required.issubset({x.name for x in p.iterdir()}):
            return p.resolve()
    root = Path("/root/nmda2")
    if root.exists():
        for marker in root.rglob("10_scientific_decision.json"):
            q = marker.parent
            if required.issubset({x.name for x in q.iterdir()}):
                return q.resolve()
    raise FileNotFoundError("Could not auto-discover Step-27 results")


@dataclass
class Inputs:
    step25: Path
    step27: Path
    audit25: dict[str, Any]
    decision25: dict[str, Any]
    decision27: dict[str, Any]
    targets: dict[str, Any]
    confirmed: pd.DataFrame
    accepted: pd.DataFrame
    geometry: pd.DataFrame
    pairwise: pd.DataFrame
    attempts27: pd.DataFrame
    network27: pd.DataFrame
    components27: pd.DataFrame
    threshold: float
    cutoff_mid: float
    cutoff_strict: float


def load_inputs(step25: Path, step27: Path) -> Inputs:
    req25 = [
        "00_input_audit.json",
        "01_corrected_target_reconstruction.json",
        "02_corrected_accepted_cohorts.csv.gz",
        "05_strong_candidates_confirmed.csv",
        "08_geometry_coordinates.csv.gz",
        "09_geometry_pairwise_distances.csv",
        "15_scientific_decision.json",
    ]
    req27 = [
        "04_bridge_attempts.csv",
        "07_verified_network_edges.csv",
        "08_final_components.csv",
        "10_scientific_decision.json",
    ]
    for name in req25:
        if not (step25 / name).exists():
            raise FileNotFoundError(step25 / name)
    for name in req27:
        if not (step27 / name).exists():
            raise FileNotFoundError(step27 / name)

    audit25 = json.loads((step25 / "00_input_audit.json").read_text(encoding="utf-8"))
    decision25 = json.loads((step25 / "15_scientific_decision.json").read_text(encoding="utf-8"))
    decision27 = json.loads((step27 / "10_scientific_decision.json").read_text(encoding="utf-8"))
    targets = json.loads((step25 / "01_corrected_target_reconstruction.json").read_text(encoding="utf-8"))
    confirmed = pd.read_csv(step25 / "05_strong_candidates_confirmed.csv")
    accepted = pd.read_csv(step25 / "02_corrected_accepted_cohorts.csv.gz")
    accepted = accepted[accepted["prior"].astype(str).str.lower().eq("broad")].copy()
    geometry = pd.read_csv(step25 / "08_geometry_coordinates.csv.gz")
    pairwise = pd.read_csv(step25 / "09_geometry_pairwise_distances.csv")
    attempts27 = pd.read_csv(step27 / "04_bridge_attempts.csv")
    network27 = pd.read_csv(step27 / "07_verified_network_edges.csv")
    components27 = pd.read_csv(step27 / "08_final_components.csv")

    if decision25.get("status") != EXPECTED_STEP25_STATUS:
        raise RuntimeError(f"Unexpected Step25 status: {decision25.get('status')!r}")
    if decision27.get("status") != EXPECTED_STEP27_STATUS:
        raise RuntimeError(f"Unexpected Step27 status: {decision27.get('status')!r}")
    if len(confirmed) != 40 or not confirmed["strong_confirmed_dt0p05"].astype(bool).all():
        raise RuntimeError("Step25 confirmed strong table is not the expected 40/40 cohort")
    if len(accepted) != 5000:
        raise RuntimeError(f"Expected 5000 corrected accepted broad rows; got {len(accepted)}")
    if not network27["strict_cutoff_pass"].astype(bool).all():
        raise RuntimeError("Step27 verified network contains an edge failing the strict cutoff")

    strict = components27[components27["criterion"].astype(str).eq("verified_strict")].copy()
    groups = [frozenset(g["sample_id"].astype(int).tolist()) for _, g in strict.groupby("component_id")]
    groups = sorted(groups, key=lambda x: (-len(x), min(x)))
    expected = sorted(EXPECTED_STRICT_COMPONENTS, key=lambda x: (-len(x), min(x)))
    if groups != expected:
        raise RuntimeError(f"Step27 strict component freeze mismatch: got {[sorted(x) for x in groups]}")

    threshold = float(audit25["strong_threshold"])
    cutoff_mid = float(audit25["corrected_broad_cutoff_midpoint"])
    cutoff_strict = float(audit25["corrected_broad_max_accepted_score"])
    if abs(float(decision27["strong_threshold"]) - threshold) > 1e-12:
        raise RuntimeError("Step25/27 strong-threshold mismatch")
    if abs(float(decision27["midpoint_cutoff"]) - cutoff_mid) > 1e-12:
        raise RuntimeError("Step25/27 midpoint-cutoff mismatch")
    if abs(float(decision27["strict_cutoff"]) - cutoff_strict) > 1e-12:
        raise RuntimeError("Step25/27 strict-cutoff mismatch")

    return Inputs(
        step25=step25,
        step27=step27,
        audit25=audit25,
        decision25=decision25,
        decision27=decision27,
        targets=targets,
        confirmed=confirmed,
        accepted=accepted,
        geometry=geometry,
        pairwise=pairwise,
        attempts27=attempts27,
        network27=network27,
        components27=components27,
        threshold=threshold,
        cutoff_mid=cutoff_mid,
        cutoff_strict=cutoff_strict,
    )


def replay_gate(inp: Inputs, calib_dt: float, confirm_dt: float) -> dict[str, Any]:
    ids = [1900, 6113, 10052, 27498, 29959, 37318]
    sub = inp.confirmed.set_index("sample_id").loc[ids].reset_index()
    rates = sub[eng.RATE_COLS].to_numpy(float)
    pulse = sub["pulse_amplitude_au"].to_numpy(float)
    tau = sub["glutamate_tau_ms"].to_numpy(float)
    met = eng.calibration_ensemble(rates, pulse, tau, calib_dt)
    score = eng.corrected_score_from_metrics(met[:, 1], met[:, 2], met[:, 3], met[:, 5], inp.targets)
    score_err = np.abs(score - sub["calibration_score_corrected"].to_numpy(float))
    r_err = []
    for _, row in sub.iterrows():
        out, _ = eng.run_base_condition(
            row[eng.RATE_COLS].to_numpy(float)[None, :],
            float(row["G_peak"]),
            float(row["tau_G_ms"]),
            confirm_dt,
        )
        r_err.append(abs(float(out[0, 2]) - float(row["r_plateau_dt0p05"])))
    max_score_err = float(np.max(score_err))
    max_r_err = float(np.max(r_err))
    if max_score_err > 1e-8 or max_r_err > 1e-8:
        raise RuntimeError(f"Engine replay gate failed: score={max_score_err:g}, response={max_r_err:g}")
    return {
        "engine_replay_gate_pass": True,
        "replayed_sample_ids": ids,
        "max_abs_corrected_score_error": max_score_err,
        "max_abs_response_ratio_error": max_r_err,
    }


def load_full_broad_bounds(inp: Inputs, allow_fallback: bool) -> tuple[dict[str, np.ndarray], str]:
    candidates = [Path("/root/nmda/results_step10.zip"), Path("/root/nmda/results_step10_tables.zip")]
    for p in candidates:
        if not p.exists():
            continue
        try:
            df, src = eng.read_step10(p, None)
            broad = df[df["prior"].astype(str).str.lower().eq("broad")].copy()
            if len(broad) >= 40000:
                rate_log = np.log(broad[eng.RATE_COLS].to_numpy(float))
                aux_log = np.log(broad[eng.AUX_COLS].to_numpy(float))
                return {
                    "rate_log_min": rate_log.min(axis=0),
                    "rate_log_max": rate_log.max(axis=0),
                    "aux_log_min": aux_log.min(axis=0),
                    "aux_log_max": aux_log.max(axis=0),
                }, src
        except Exception:
            pass
    if not allow_fallback:
        raise FileNotFoundError(
            "Full broad Step10 candidate table was not found. Expected /root/nmda/results_step10.zip. "
            "Use --allow-accepted-bounds-fallback only for diagnostic reuse."
        )
    broad = inp.accepted
    rate_log = np.log(broad[eng.RATE_COLS].to_numpy(float))
    aux_log = np.log(broad[eng.AUX_COLS].to_numpy(float))
    return {
        "rate_log_min": rate_log.min(axis=0),
        "rate_log_max": rate_log.max(axis=0),
        "aux_log_min": aux_log.min(axis=0),
        "aux_log_max": aux_log.max(axis=0),
    }, "accepted-broad empirical fallback"


def geometry_stats(inp: Inputs, bounds: dict[str, np.ndarray]) -> dict[str, np.ndarray]:
    rate_logs = np.log(inp.accepted[eng.RATE_COLS].to_numpy(float))
    rmean = rate_logs.mean(axis=0)
    rsd = rate_logs.std(axis=0, ddof=1)
    aux_logs = np.log(inp.accepted[eng.AUX_COLS].to_numpy(float))
    amean = aux_logs.mean(axis=0)
    asd = aux_logs.std(axis=0, ddof=1)
    if np.any(rsd <= 0) or np.any(asd <= 0):
        raise RuntimeError("Zero variance in accepted-broad standardized coordinates")
    return {
        "rate_log_mean": rmean,
        "rate_log_sd": rsd,
        "aux_log_mean": amean,
        "aux_log_sd": asd,
        "rate_z_min": (bounds["rate_log_min"] - rmean) / rsd,
        "rate_z_max": (bounds["rate_log_max"] - rmean) / rsd,
        "aux_z_min": (bounds["aux_log_min"] - amean) / asd,
        "aux_z_max": (bounds["aux_log_max"] - amean) / asd,
    }


class UnionFind:
    def __init__(self, nodes: list[int]):
        self.parent = {int(x): int(x) for x in nodes}
        self.rank = {int(x): 0 for x in nodes}

    def find(self, x: int) -> int:
        x = int(x)
        if self.parent[x] != x:
            self.parent[x] = self.find(self.parent[x])
        return self.parent[x]

    def union(self, a: int, b: int) -> bool:
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return False
        if self.rank[ra] < self.rank[rb]:
            ra, rb = rb, ra
        self.parent[rb] = ra
        if self.rank[ra] == self.rank[rb]:
            self.rank[ra] += 1
        return True

    def groups(self) -> list[list[int]]:
        d: dict[int, list[int]] = {}
        for x in self.parent:
            d.setdefault(self.find(x), []).append(x)
        return sorted((sorted(v) for v in d.values()), key=lambda v: (-len(v), v[0]))


def strict_component_map(inp: Inputs) -> dict[int, int]:
    s = inp.components27[inp.components27["criterion"].astype(str).eq("verified_strict")].copy()
    return {int(r.sample_id): int(r.component_id) for r in s.itertuples(index=False)}


def target_forcing_nodes(inp: Inputs) -> list[tuple[float, float]]:
    nodes = list(dict.fromkeys(
        (float(g), float(t))
        for g, t in inp.confirmed[["G_peak", "tau_G_ms"]].itertuples(index=False, name=None)
    ))
    return nodes


def full_grid_nodes() -> list[tuple[float, float]]:
    return [(float(g), float(t)) for g in eng.G_PEAK_GRID for t in eng.TAU_GRID]


def calibration_eval(rates: np.ndarray, pulse: np.ndarray, tau: np.ndarray, inp: Inputs, dt: float) -> dict[str, np.ndarray]:
    met = eng.calibration_ensemble(rates, pulse.astype(float), tau.astype(float), dt)
    score = eng.corrected_score_from_metrics(met[:, 1], met[:, 2], met[:, 3], met[:, 5], inp.targets)
    valid = met[:, 8] > 0.5
    score[~valid] = np.inf
    return {
        "score": score,
        "H": met[:, 1],
        "E": met[:, 2],
        "L": met[:, 3],
        "J_full": met[:, 4],
        "J_corr": met[:, 5],
        "valid": valid,
        "minimum_state": met[:, 7],
    }


def response_scan_nodes(rates: np.ndarray, nodes: list[tuple[float, float]], dt: float) -> dict[str, np.ndarray]:
    n = len(rates)
    best_r = np.full(n, -np.inf, dtype=float)
    best_control = np.full(n, np.nan, dtype=float)
    best_blocked = np.full(n, np.nan, dtype=float)
    best_g = np.full(n, np.nan, dtype=float)
    best_tau = np.full(n, np.nan, dtype=float)
    min_state_c = np.full(n, np.nan, dtype=float)
    min_state_b = np.full(n, np.nan, dtype=float)
    for g, tau in nodes:
        out, _ = eng.run_base_condition(rates, float(g), float(tau), dt)
        r = out[:, 2]
        better = np.isfinite(r) & (r > best_r)
        if np.any(better):
            best_r[better] = r[better]
            best_control[better] = out[better, 0]
            best_blocked[better] = out[better, 1]
            best_g[better] = float(g)
            best_tau[better] = float(tau)
            min_state_c[better] = out[better, 3]
            min_state_b[better] = out[better, 4]
    best_r[~np.isfinite(best_r)] = np.nan
    return {
        "r": best_r,
        "control": best_control,
        "blocked": best_blocked,
        "G_peak": best_g,
        "tau_G_ms": best_tau,
        "min_state_control": min_state_c,
        "min_state_blocked": min_state_b,
    }


def confirm_best_nodes(rates: np.ndarray, coarse: dict[str, np.ndarray], dt: float, threshold: float) -> dict[str, np.ndarray]:
    n = len(rates)
    result = {k: np.array(v, copy=True) for k, v in coarse.items()}
    r_confirm = np.full(n, np.nan, dtype=float)
    c_confirm = np.full(n, np.nan, dtype=float)
    b_confirm = np.full(n, np.nan, dtype=float)
    minc = np.full(n, np.nan, dtype=float)
    minb = np.full(n, np.nan, dtype=float)
    unique = sorted({
        (float(coarse["G_peak"][i]), float(coarse["tau_G_ms"][i]))
        for i in range(n) if np.isfinite(coarse["G_peak"][i])
    })
    for g, tau in unique:
        idx = np.where((np.abs(coarse["G_peak"] - g) < 1e-12) & (np.abs(coarse["tau_G_ms"] - tau) < 1e-12))[0]
        out, _ = eng.run_base_condition(rates[idx], g, tau, dt)
        r_confirm[idx] = out[:, 2]
        c_confirm[idx] = out[:, 0]
        b_confirm[idx] = out[:, 1]
        minc[idx] = out[:, 3]
        minb[idx] = out[:, 4]
    fail = ~np.isfinite(r_confirm) | (r_confirm < threshold)
    if np.any(fail):
        exact = response_scan_nodes(rates[fail], full_grid_nodes(), dt)
        r_confirm[fail] = exact["r"]
        c_confirm[fail] = exact["control"]
        b_confirm[fail] = exact["blocked"]
        minc[fail] = exact["min_state_control"]
        minb[fail] = exact["min_state_blocked"]
        result["G_peak"][fail] = exact["G_peak"]
        result["tau_G_ms"][fail] = exact["tau_G_ms"]
    result["r"] = r_confirm
    result["control"] = c_confirm
    result["blocked"] = b_confirm
    result["min_state_control"] = minc
    result["min_state_blocked"] = minb
    return result


def row_rate_z(row: pd.Series, geom: dict[str, np.ndarray]) -> np.ndarray:
    return (np.log(row[eng.RATE_COLS].to_numpy(float)) - geom["rate_log_mean"]) / geom["rate_log_sd"]


def row_aux_z(row: pd.Series, geom: dict[str, np.ndarray]) -> np.ndarray:
    return (np.log(row[eng.AUX_COLS].to_numpy(float)) - geom["aux_log_mean"]) / geom["aux_log_sd"]


def variable_bounds(n_waypoints: int, mode: str, geom: dict[str, np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    lo = np.tile(geom["rate_z_min"], n_waypoints)
    hi = np.tile(geom["rate_z_max"], n_waypoints)
    if mode == "flex_aux":
        lo = np.concatenate([lo, np.tile(geom["aux_z_min"], n_waypoints)])
        hi = np.concatenate([hi, np.tile(geom["aux_z_max"], n_waypoints)])
    elif mode != "native_aux":
        raise ValueError(mode)
    return lo.astype(float), hi.astype(float)


def linear_internal_variables(row_a: pd.Series, row_b: pd.Series, n_waypoints: int, mode: str, geom: dict[str, np.ndarray]) -> np.ndarray:
    za, zb = row_rate_z(row_a, geom), row_rate_z(row_b, geom)
    frac = np.arange(1, n_waypoints + 1, dtype=float) / (n_waypoints + 1.0)
    internal = np.vstack([(1.0 - f) * za + f * zb for f in frac]).reshape(-1)
    if mode == "flex_aux":
        aa, ab = row_aux_z(row_a, geom), row_aux_z(row_b, geom)
        aux = np.vstack([(1.0 - f) * aa + f * ab for f in frac]).reshape(-1)
        internal = np.concatenate([internal, aux])
    return internal.astype(float)


def path_from_variables(
    row_a: pd.Series,
    row_b: pd.Series,
    variables: np.ndarray,
    lam: np.ndarray,
    geom: dict[str, np.ndarray],
    n_waypoints: int,
    mode: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    variables = np.asarray(variables, dtype=float)
    nr = n_waypoints * len(eng.RATE_COLS)
    expected = nr + (n_waypoints * 2 if mode == "flex_aux" else 0)
    if variables.size != expected:
        raise ValueError(f"Expected {expected} waypoint variables, got {variables.size}")

    za, zb = row_rate_z(row_a, geom), row_rate_z(row_b, geom)
    z_internal = variables[:nr].reshape(n_waypoints, len(eng.RATE_COLS))
    knots = np.linspace(0.0, 1.0, n_waypoints + 2)
    z_knots = np.vstack([za, z_internal, zb])
    z = np.column_stack([np.interp(lam, knots, z_knots[:, j]) for j in range(z_knots.shape[1])])
    rate_log = geom["rate_log_mean"][None, :] + z * geom["rate_log_sd"][None, :]
    rates = np.exp(rate_log)

    aux_a_log = np.log(row_a[eng.AUX_COLS].to_numpy(float))
    aux_b_log = np.log(row_b[eng.AUX_COLS].to_numpy(float))
    if mode == "native_aux":
        aux_log = (1.0 - lam[:, None]) * aux_a_log[None, :] + lam[:, None] * aux_b_log[None, :]
        aux_z_path = (aux_log - geom["aux_log_mean"][None, :]) / geom["aux_log_sd"][None, :]
    elif mode == "flex_aux":
        aux_internal = variables[nr:].reshape(n_waypoints, 2)
        aa, ab = row_aux_z(row_a, geom), row_aux_z(row_b, geom)
        aux_knots = np.vstack([aa, aux_internal, ab])
        aux_z_path = np.column_stack([np.interp(lam, knots, aux_knots[:, j]) for j in range(2)])
        aux_log = geom["aux_log_mean"][None, :] + aux_z_path * geom["aux_log_sd"][None, :]
    else:
        raise ValueError(mode)
    aux = np.exp(aux_log)
    return rates, aux[:, 0], aux[:, 1], z, aux_z_path


def path_stretch(z: np.ndarray) -> float:
    direct = float(np.linalg.norm(z[-1] - z[0]))
    length = float(np.sum(np.linalg.norm(np.diff(z, axis=0), axis=1)))
    if direct <= 1e-12:
        return 1.0 if length <= 1e-12 else 1e6
    return length / direct


def verify_path(
    row_a: pd.Series,
    row_b: pd.Series,
    variables: np.ndarray,
    n_waypoints: int,
    mode: str,
    inp: Inputs,
    geom: dict[str, np.ndarray],
    n_points: int,
    grid_dt: float,
    confirm_dt: float,
    calib_dt: float,
) -> tuple[dict[str, Any], pd.DataFrame]:
    lam = np.linspace(0.0, 1.0, n_points)
    rates, pulse, tau, z, aux_z = path_from_variables(row_a, row_b, variables, lam, geom, n_waypoints, mode)
    calib = calibration_eval(rates, pulse, tau, inp, calib_dt)
    coarse = response_scan_nodes(rates, full_grid_nodes(), grid_dt)
    confirmed = confirm_best_nodes(rates, coarse, confirm_dt, inp.threshold)

    strong = np.isfinite(confirmed["r"]) & (confirmed["r"] >= inp.threshold)
    compat_mid = np.isfinite(calib["score"]) & (calib["score"] <= inp.cutoff_mid)
    compat_strict = np.isfinite(calib["score"]) & (calib["score"] <= inp.cutoff_strict)
    pass_mid = bool(strong.all() and compat_mid.all())
    pass_strict = bool(strong.all() and compat_strict.all())

    sid_a, sid_b = int(row_a["sample_id"]), int(row_b["sample_id"])
    summary = {
        "sample_id_a": sid_a,
        "sample_id_b": sid_b,
        "mode": mode,
        "n_waypoints": int(n_waypoints),
        "n_points": int(n_points),
        "max_calibration_score": float(np.nanmax(calib["score"])),
        "min_r_dt0p05": float(np.nanmin(confirmed["r"])),
        "lambda_at_max_score": float(lam[int(np.nanargmax(calib["score"]))]),
        "lambda_at_min_r": float(lam[int(np.nanargmin(confirmed["r"]))]),
        "min_control_plateau_at_support_node_dt0p05": float(np.nanmin(confirmed["control"])),
        "midpoint_cutoff_pass": pass_mid,
        "strict_cutoff_pass": pass_strict,
        "strong_pass": bool(strong.all()),
        "path_stretch_standardized_log_rate": float(path_stretch(z)),
    }
    frame = pd.DataFrame({
        "sample_id_a": sid_a,
        "sample_id_b": sid_b,
        "mode": mode,
        "n_waypoints": int(n_waypoints),
        "lambda": lam,
        "calibration_score": calib["score"],
        "calibration_H": calib["H"],
        "calibration_E": calib["E"],
        "calibration_L": calib["L"],
        "calibration_J_corrected_140_500_ms": calib["J_corr"],
        "calibration_J_full_0_600_ms": calib["J_full"],
        "r_plateau_dt0p05": confirmed["r"],
        "control_plateau_PO_dt0p05": confirmed["control"],
        "blocked_plateau_PO_dt0p05": confirmed["blocked"],
        "G_peak": confirmed["G_peak"],
        "tau_G_ms": confirmed["tau_G_ms"],
        "minimum_state_control": confirmed["min_state_control"],
        "minimum_state_blocked": confirmed["min_state_blocked"],
        "compatible_midpoint": compat_mid,
        "compatible_strict": compat_strict,
        "strong": strong,
    })
    for k, col in enumerate(eng.RATE_COLS):
        frame[col] = rates[:, k]
    frame["pulse_amplitude_au"] = pulse
    frame["glutamate_tau_ms"] = tau
    return summary, frame


def component_priority(inp: Inputs, comp_map: dict[int, int]) -> pd.DataFrame:
    pair_lookup = {}
    for r in inp.pairwise.itertuples(index=False):
        a, b = int(r.sample_id_a), int(r.sample_id_b)
        key = tuple(sorted((a, b)))
        pair_lookup[key] = float(r.standardized_log_rate_distance)

    prior_barrier_by_pair: dict[tuple[int, int], float] = {}
    for r in inp.attempts27.itertuples(index=False):
        a, b = int(r.sample_id_a), int(r.sample_id_b)
        vals = []
        sj = getattr(r, "search_joint_barrier", np.nan)
        if np.isfinite(sj):
            vals.append(float(sj))
        ms = getattr(r, "max_calibration_score", np.nan)
        mr = getattr(r, "min_r_dt0p05", np.nan)
        if np.isfinite(ms) and np.isfinite(mr) and mr > 0:
            vals.append(max(float(ms) / inp.cutoff_strict, inp.threshold / float(mr)))
        if vals:
            key = tuple(sorted((a, b)))
            prior_barrier_by_pair[key] = min(prior_barrier_by_pair.get(key, np.inf), min(vals))

    rows = []
    comps = sorted(set(comp_map.values()))
    for i, ca in enumerate(comps):
        for cb in comps[i + 1:]:
            members_a = [sid for sid, c in comp_map.items() if c == ca]
            members_b = [sid for sid, c in comp_map.items() if c == cb]
            best_dist = np.inf
            best_prior = np.inf
            best_pair = None
            for a in members_a:
                for b in members_b:
                    key = tuple(sorted((a, b)))
                    d = pair_lookup.get(key, np.inf)
                    pb = prior_barrier_by_pair.get(key, np.inf)
                    if d < best_dist:
                        best_dist = d
                        best_pair = key
                    best_prior = min(best_prior, pb)
            rows.append({
                "component_a": int(ca),
                "component_b": int(cb),
                "size_a": len(members_a),
                "size_b": len(members_b),
                "best_step27_joint_barrier": float(best_prior) if np.isfinite(best_prior) else np.nan,
                "min_standardized_log_rate_distance": float(best_dist),
                "nearest_sample_id_a": int(best_pair[0]) if best_pair else np.nan,
                "nearest_sample_id_b": int(best_pair[1]) if best_pair else np.nan,
            })
    out = pd.DataFrame(rows)
    out["priority_barrier"] = out["best_step27_joint_barrier"].fillna(np.inf)
    out = out.sort_values(["priority_barrier", "min_standardized_log_rate_distance", "component_a", "component_b"]).drop(columns="priority_barrier")
    return out.reset_index(drop=True)


def endpoint_pairs_for_components(
    inp: Inputs,
    comp_map: dict[int, int],
    ca: int,
    cb: int,
    n: int,
) -> pd.DataFrame:
    pw = inp.pairwise.copy()
    dlookup = {
        tuple(sorted((int(r.sample_id_a), int(r.sample_id_b)))): float(r.standardized_log_rate_distance)
        for r in pw.itertuples(index=False)
    }
    prior = {}
    for r in inp.attempts27.itertuples(index=False):
        a, b = int(r.sample_id_a), int(r.sample_id_b)
        if {comp_map[a], comp_map[b]} != {int(ca), int(cb)}:
            continue
        vals = []
        sj = getattr(r, "search_joint_barrier", np.nan)
        if np.isfinite(sj):
            vals.append(float(sj))
        ms = getattr(r, "max_calibration_score", np.nan)
        mr = getattr(r, "min_r_dt0p05", np.nan)
        if np.isfinite(ms) and np.isfinite(mr) and mr > 0:
            vals.append(max(float(ms) / inp.cutoff_strict, inp.threshold / float(mr)))
        if vals:
            key = tuple(sorted((a, b)))
            prior[key] = min(prior.get(key, np.inf), min(vals))

    members_a = sorted([sid for sid, c in comp_map.items() if c == ca])
    members_b = sorted([sid for sid, c in comp_map.items() if c == cb])
    candidates = []
    for a in members_a:
        for b in members_b:
            key = tuple(sorted((a, b)))
            candidates.append({
                "sample_id_a": int(a),
                "sample_id_b": int(b),
                "distance": float(dlookup.get(key, np.inf)),
                "step27_best_joint_barrier": float(prior.get(key, np.nan)),
            })
    df = pd.DataFrame(candidates)
    chosen = []
    seen = set()
    # First preserve the best Step27 near-feasible attempts.
    prior_df = df[np.isfinite(df["step27_best_joint_barrier"])].sort_values(["step27_best_joint_barrier", "distance"])
    for r in prior_df.itertuples(index=False):
        key = tuple(sorted((int(r.sample_id_a), int(r.sample_id_b))))
        if key not in seen:
            chosen.append(r._asdict())
            seen.add(key)
        if len(chosen) >= max(1, n // 2):
            break
    # Then nearest pairs.
    for r in df.sort_values(["distance", "sample_id_a", "sample_id_b"]).itertuples(index=False):
        key = tuple(sorted((int(r.sample_id_a), int(r.sample_id_b))))
        if key not in seen:
            chosen.append(r._asdict())
            seen.add(key)
        if len(chosen) >= n:
            break
    return pd.DataFrame(chosen).head(n).reset_index(drop=True)


def seed_variables(
    row_a: pd.Series,
    row_b: pd.Series,
    inp: Inputs,
    geom: dict[str, np.ndarray],
    n_waypoints: int,
    mode: str,
    max_seeds: int,
) -> list[np.ndarray]:
    seeds: list[np.ndarray] = []
    linear = linear_internal_variables(row_a, row_b, n_waypoints, mode, geom)
    seeds.append(linear)

    za, zb = row_rate_z(row_a, geom), row_rate_z(row_b, geom)
    frac = np.arange(1, n_waypoints + 1, dtype=float) / (n_waypoints + 1.0)
    base = np.vstack([(1 - f) * za + f * zb for f in frac])

    strong = inp.confirmed.copy()
    strong_z = (np.log(strong[eng.RATE_COLS].to_numpy(float)) - geom["rate_log_mean"]) / geom["rate_log_sd"]
    strong_aux = (np.log(strong[eng.AUX_COLS].to_numpy(float)) - geom["aux_log_mean"]) / geom["aux_log_sd"]
    for rank_shift in range(0, min(6, max_seeds)):
        zlist = []
        alist = []
        used: set[int] = set()
        for j in range(n_waypoints):
            d = np.linalg.norm(strong_z - base[j][None, :], axis=1)
            order = np.argsort(d)
            pick = None
            for idx in order[rank_shift:]:
                sid = int(strong.iloc[int(idx)]["sample_id"])
                if sid not in used and sid not in {int(row_a["sample_id"]), int(row_b["sample_id"])}:
                    pick = int(idx)
                    break
            if pick is None:
                pick = int(order[0])
            used.add(int(strong.iloc[pick]["sample_id"]))
            zlist.append(strong_z[pick])
            alist.append(strong_aux[pick])
        var = np.asarray(zlist, float).reshape(-1)
        if mode == "flex_aux":
            var = np.concatenate([var, np.asarray(alist, float).reshape(-1)])
        seeds.append(var)
        if len(seeds) >= max_seeds:
            break

    # Calibration-manifold seed: nearest corrected accepted vector to each direct internal knot.
    acc_z = (np.log(inp.accepted[eng.RATE_COLS].to_numpy(float)) - geom["rate_log_mean"]) / geom["rate_log_sd"]
    acc_aux = (np.log(inp.accepted[eng.AUX_COLS].to_numpy(float)) - geom["aux_log_mean"]) / geom["aux_log_sd"]
    zlist = []
    alist = []
    for j in range(n_waypoints):
        idx = int(np.argmin(np.linalg.norm(acc_z - base[j][None, :], axis=1)))
        zlist.append(acc_z[idx])
        alist.append(acc_aux[idx])
    var = np.asarray(zlist, float).reshape(-1)
    if mode == "flex_aux":
        var = np.concatenate([var, np.asarray(alist, float).reshape(-1)])
    seeds.append(var)
    return seeds[:max_seeds]


def evaluate_population_calibration(
    pop: np.ndarray,
    row_a: pd.Series,
    row_b: pd.Series,
    lam: np.ndarray,
    inp: Inputs,
    geom: dict[str, np.ndarray],
    n_waypoints: int,
    mode: str,
    dt: float,
) -> tuple[np.ndarray, list[dict[str, Any]]]:
    all_rates = []
    all_pulse = []
    all_tau = []
    cache: list[dict[str, Any]] = []
    for var in pop:
        rates, pulse, tau, z, _ = path_from_variables(row_a, row_b, var, lam, geom, n_waypoints, mode)
        all_rates.append(rates)
        all_pulse.append(pulse)
        all_tau.append(tau)
        cache.append({"rates": rates, "pulse": pulse, "tau": tau, "z": z, "stretch": path_stretch(z)})
    nr = len(lam)
    met = calibration_eval(
        np.vstack(all_rates),
        np.concatenate(all_pulse),
        np.concatenate(all_tau),
        inp,
        dt,
    )
    scores = met["score"].reshape(len(pop), nr)
    max_score = np.nanmax(scores, axis=1)
    for i, c in enumerate(cache):
        c["max_score"] = float(max_score[i])
    return max_score, cache


def evaluate_selected_response(
    cache: list[dict[str, Any]],
    indices: list[int],
    nodes: list[tuple[float, float]],
    dt: float,
) -> dict[int, float]:
    if not indices:
        return {}
    npath = len(cache[indices[0]]["rates"])
    rates = np.vstack([cache[i]["rates"] for i in indices])
    scan = response_scan_nodes(rates, nodes, dt)
    rr = scan["r"].reshape(len(indices), npath)
    return {idx: float(np.nanmin(rr[j])) for j, idx in enumerate(indices)}


def search_waypoint_path(
    row_a: pd.Series,
    row_b: pd.Series,
    inp: Inputs,
    geom: dict[str, np.ndarray],
    args: argparse.Namespace,
    n_waypoints: int,
    mode: str,
    seed_offset: int,
) -> dict[str, Any]:
    lo, hi = variable_bounds(n_waypoints, mode, geom)
    base = linear_internal_variables(row_a, row_b, n_waypoints, mode, geom)
    base = np.clip(base, lo, hi)
    seed_pool = [np.clip(x, lo, hi) for x in seed_variables(row_a, row_b, inp, geom, n_waypoints, mode, max_seeds=8)]
    target_nodes = target_forcing_nodes(inp)
    lam = np.linspace(0.0, 1.0, args.search_points)
    dim = len(base)

    global_best_var = base.copy()
    global_best_joint = np.inf
    global_best_score = np.inf
    global_best_minr = 0.0
    history = []

    for restart in range(args.restarts):
        rng = np.random.default_rng(args.seed + seed_offset + 1009 * restart + 37 * n_waypoints + (0 if mode == "native_aux" else 100000))
        if restart == 0:
            mean = base.copy()
        else:
            mean = global_best_var.copy() if np.isfinite(global_best_joint) else base.copy()
            if len(seed_pool) > 1:
                alt = seed_pool[min(restart, len(seed_pool) - 1)]
                mean = 0.5 * mean + 0.5 * alt
        span = np.maximum(hi - lo, 1e-6)
        sigma = np.minimum(args.initial_sigma, 0.25 * span)
        sigma = np.maximum(sigma, 0.15)

        best_var = mean.copy()
        best_joint = np.inf
        stale = 0
        for gen in range(args.generations):
            pop = rng.normal(mean[None, :], sigma[None, :], size=(args.population, dim))
            pop = np.clip(pop, lo[None, :], hi[None, :])
            pop[0] = np.clip(best_var, lo, hi)
            if gen == 0:
                for k, seedv in enumerate(seed_pool[: min(len(seed_pool), args.population - 1)], start=1):
                    pop[k] = seedv
            elif gen % 5 == 0 and len(seed_pool) > 1:
                for k in range(min(3, len(seed_pool))):
                    pop[-1 - k] = seed_pool[(gen // 5 + k) % len(seed_pool)]

            max_score, cache = evaluate_population_calibration(
                pop, row_a, row_b, lam, inp, geom, n_waypoints, mode, args.calibration_dt_ms
            )
            score_ratio = max_score / inp.cutoff_strict
            pre = score_ratio + np.array([0.002 * max(0.0, c["stretch"] - 1.0) for c in cache])
            eligible = np.where(np.isfinite(pre) & (score_ratio <= max(args.flex_trigger, 2.0)))[0]
            if len(eligible) == 0:
                selected = list(np.argsort(pre)[: args.response_top])
            else:
                order = eligible[np.argsort(pre[eligible])]
                selected = list(order[: args.response_top])
                # Add a small random exploratory subset when room remains.
                remain = [i for i in eligible.tolist() if i not in selected]
                while len(selected) < args.response_top and remain:
                    j = int(rng.integers(0, len(remain)))
                    selected.append(remain.pop(j))

            minr = evaluate_selected_response(cache, selected, target_nodes, args.grid_dt_ms)
            joint = np.full(args.population, np.inf, dtype=float)
            for i in selected:
                r = float(minr.get(i, 0.0))
                if not np.isfinite(r) or r <= 0:
                    continue
                stretch_pen = 0.002 * max(0.0, cache[i]["stretch"] - 1.0)
                joint[i] = max(score_ratio[i], inp.threshold / r) + stretch_pen

            order = np.argsort(joint)
            if np.isfinite(joint[order[0]]) and joint[order[0]] < best_joint:
                best_joint = float(joint[order[0]])
                best_var = pop[order[0]].copy()
                stale = 0
            else:
                stale += 1
            if best_joint < global_best_joint:
                global_best_joint = best_joint
                global_best_var = best_var.copy()
                ib = int(order[0])
                global_best_score = float(max_score[ib]) if np.isfinite(joint[ib]) else global_best_score
                global_best_minr = float(minr.get(ib, global_best_minr))

            elite_idx = [int(i) for i in order if np.isfinite(joint[i])][: max(4, min(7, args.response_top))]
            if elite_idx:
                elite = pop[elite_idx]
                weights = np.linspace(2.0, 1.0, len(elite_idx))
                weights /= weights.sum()
                mean = np.sum(elite * weights[:, None], axis=0)
                esd = elite.std(axis=0, ddof=0)
                sigma = np.maximum(0.08, 0.68 * sigma + 0.32 * esd)
                sigma *= 0.95
            else:
                sigma = np.minimum(1.5, sigma * 1.08)

            history.append({
                "restart": restart,
                "generation": gen,
                "mode": mode,
                "n_waypoints": n_waypoints,
                "best_joint_targeted": float(best_joint),
                "global_best_joint_targeted": float(global_best_joint),
                "stale_generations": stale,
            })
            if global_best_joint <= 0.995 and stale >= 3:
                break

    # One full frozen-grid screen of the best search candidate at search resolution.
    rates, pulse, tau, z, _ = path_from_variables(
        row_a, row_b, global_best_var, lam, geom, n_waypoints, mode
    )
    cal = calibration_eval(rates, pulse, tau, inp, args.calibration_dt_ms)
    full = response_scan_nodes(rates, full_grid_nodes(), args.grid_dt_ms)
    max_score = float(np.nanmax(cal["score"]))
    min_r_full = float(np.nanmin(full["r"]))
    full_joint = max(max_score / inp.cutoff_strict, inp.threshold / max(min_r_full, 1e-300))
    return {
        "variables": global_best_var,
        "search_max_score": max_score,
        "search_min_r_targeted_dt0p1": float(global_best_minr),
        "search_min_r_fullgrid_dt0p1": min_r_full,
        "search_joint_barrier_fullgrid": float(full_joint),
        "search_path_stretch": float(path_stretch(z)),
        "history": history,
    }


def checkpoint_key(a: int, b: int, mode: str, n_waypoints: int) -> str:
    x, y = sorted((int(a), int(b)))
    return f"bridge_{x}_{y}_{mode}_wp{int(n_waypoints)}"


def save_stage_checkpoint(chkdir: Path, key: str, result: dict[str, Any], verify: dict[str, Any] | None, frame: pd.DataFrame | None) -> None:
    d = chkdir / key
    d.mkdir(parents=True, exist_ok=True)
    vars_ = np.asarray(result["variables"], dtype=float)
    np.save(d / "variables.npy", vars_)
    payload = {k: v for k, v in result.items() if k not in {"variables", "history"}}
    payload["history"] = result.get("history", [])
    payload["verify"] = verify
    json_dump(d / "summary.json", payload)
    if frame is not None:
        frame.to_csv(d / "verified_path_points.csv.gz", index=False, compression="gzip")


def load_stage_checkpoint(chkdir: Path, key: str) -> tuple[dict[str, Any], dict[str, Any] | None, pd.DataFrame | None] | None:
    d = chkdir / key
    if not (d / "variables.npy").exists() or not (d / "summary.json").exists():
        return None
    vars_ = np.load(d / "variables.npy")
    payload = json.loads((d / "summary.json").read_text(encoding="utf-8"))
    verify = payload.pop("verify", None)
    frame = pd.read_csv(d / "verified_path_points.csv.gz") if (d / "verified_path_points.csv.gz").exists() else None
    payload["variables"] = vars_
    return payload, verify, frame


def components_table(ids: list[int], edges: pd.DataFrame, pass_col: str, criterion: str) -> pd.DataFrame:
    uf = UnionFind(ids)
    if len(edges):
        for r in edges.itertuples(index=False):
            if bool(getattr(r, pass_col)):
                uf.union(int(r.sample_id_a), int(r.sample_id_b))
    rows = []
    for cid, members in enumerate(uf.groups()):
        for sid in members:
            rows.append({"criterion": criterion, "component_id": cid, "component_size": len(members), "sample_id": sid})
    return pd.DataFrame(rows)


def make_figures(inp: Inputs, final_edges: pd.DataFrame, final_components: pd.DataFrame, verified_paths: pd.DataFrame, attempts: pd.DataFrame, outdir: Path) -> None:
    ids = sorted(inp.confirmed["sample_id"].astype(int).tolist())
    coords = inp.geometry.set_index("sample_id") if "sample_id" in inp.geometry.columns else None

    # Figure A: sampled network in first two available geometry coordinates.
    plt.figure(figsize=(7.4, 6.0))
    if coords is not None:
        numcols = [c for c in coords.columns if np.issubdtype(coords[c].dtype, np.number)]
        xcol = "PC1" if "PC1" in coords.columns else (numcols[0] if numcols else None)
        ycol = "PC2" if "PC2" in coords.columns else (numcols[1] if len(numcols) > 1 else None)
    else:
        xcol = ycol = None
    if xcol and ycol:
        for r in final_edges.itertuples(index=False):
            a, b = int(r.sample_id_a), int(r.sample_id_b)
            if a in coords.index and b in coords.index:
                plt.plot([coords.loc[a, xcol], coords.loc[b, xcol]], [coords.loc[a, ycol], coords.loc[b, ycol]], linewidth=0.8, alpha=0.55)
        comp_strict = final_components[final_components["criterion"].eq("step28_strict")]
        for cid, g in comp_strict.groupby("component_id"):
            xs = [coords.loc[int(s), xcol] for s in g["sample_id"] if int(s) in coords.index]
            ys = [coords.loc[int(s), ycol] for s in g["sample_id"] if int(s) in coords.index]
            plt.scatter(xs, ys, s=34, label=f"component {cid} (n={len(g)})")
        for sid in [6113, 27498, 29959, 37318]:
            if sid in coords.index:
                plt.annotate(str(sid), (coords.loc[sid, xcol], coords.loc[sid, ycol]), xytext=(4, 4), textcoords="offset points", fontsize=8)
        plt.xlabel(xcol)
        plt.ylabel(ycol)
        plt.legend(fontsize=8, frameon=False)
    else:
        plt.text(0.5, 0.5, "Geometry coordinates unavailable", ha="center", va="center")
        plt.axis("off")
    plt.title("Step 28 verified strong-compatible network")
    plt.tight_layout()
    plt.savefig(outdir / "Fig28A_verified_network.png", dpi=220)
    plt.savefig(outdir / "Fig28A_verified_network.pdf")
    plt.close()

    # Figure B: best barrier per endpoint-pair search stage.
    plt.figure(figsize=(8.0, 5.4))
    if len(attempts):
        a = attempts.copy()
        a["pair"] = a["sample_id_a"].astype(str) + "-" + a["sample_id_b"].astype(str)
        a = a.sort_values("search_joint_barrier_fullgrid").head(25)
        xpos = np.arange(len(a))
        plt.bar(xpos, a["search_joint_barrier_fullgrid"].to_numpy(float))
        plt.axhline(1.0, linewidth=1.0, linestyle="--")
        plt.xticks(xpos, a["pair"] + "\n" + a["mode"] + "/wp" + a["n_waypoints"].astype(str), rotation=90, fontsize=7)
        plt.ylabel("search joint barrier (1 = feasible)")
    else:
        plt.text(0.5, 0.5, "No search attempts", ha="center", va="center")
        plt.axis("off")
    plt.title("Residual bridge search barriers")
    plt.tight_layout()
    plt.savefig(outdir / "Fig28B_search_barriers.png", dpi=220)
    plt.savefig(outdir / "Fig28B_search_barriers.pdf")
    plt.close()

    # Figure C: verified bridge profiles.
    plt.figure(figsize=(8.0, 5.4))
    if len(verified_paths):
        for (a, b, mode, wp), g in verified_paths.groupby(["sample_id_a", "sample_id_b", "mode", "n_waypoints"]):
            label = f"{int(a)}-{int(b)} {mode} wp{int(wp)}"
            plt.plot(g["lambda"], g["r_plateau_dt0p05"] / inp.threshold, linewidth=1.3, label=label + " r/r*")
            plt.plot(g["lambda"], g["calibration_score"] / inp.cutoff_strict, linewidth=1.0, linestyle="--", label=label + " S/S*")
        plt.axhline(1.0, linewidth=1.0, linestyle=":")
        plt.ylabel("normalized constraint value")
        plt.xlabel("path coordinate")
        plt.legend(fontsize=7, frameon=False, ncol=1)
    else:
        plt.text(0.5, 0.5, "No new verified bridge", ha="center", va="center")
        plt.axis("off")
    plt.title("Verified Step 28 bridge profiles")
    plt.tight_layout()
    plt.savefig(outdir / "Fig28C_verified_bridge_profiles.png", dpi=220)
    plt.savefig(outdir / "Fig28C_verified_bridge_profiles.pdf")
    plt.close()


def self_test() -> None:
    # Pure geometry test; does not require project files.
    rate_cols = len(eng.RATE_COLS)
    fake = {}
    for i, c in enumerate(eng.RATE_COLS):
        fake[c] = 1.0 + 0.1 * i
    fake["pulse_amplitude_au"] = 0.5
    fake["glutamate_tau_ms"] = 10.0
    fake["sample_id"] = 1
    a = pd.Series(fake)
    b = a.copy()
    b[eng.RATE_COLS] = a[eng.RATE_COLS].to_numpy(float) * 1.5
    b["pulse_amplitude_au"] = 0.9
    b["glutamate_tau_ms"] = 20.0
    b["sample_id"] = 2
    logs = np.vstack([np.log(a[eng.RATE_COLS].to_numpy(float)), np.log(b[eng.RATE_COLS].to_numpy(float))])
    amean = np.log(np.array([0.65, 14.0]))
    geom = {
        "rate_log_mean": logs.mean(axis=0),
        "rate_log_sd": np.ones(rate_cols),
        "aux_log_mean": amean,
        "aux_log_sd": np.ones(2),
        "rate_z_min": np.full(rate_cols, -5.0),
        "rate_z_max": np.full(rate_cols, 5.0),
        "aux_z_min": np.full(2, -5.0),
        "aux_z_max": np.full(2, 5.0),
    }
    for mode in ["native_aux", "flex_aux"]:
        v = linear_internal_variables(a, b, 3, mode, geom)
        lam = np.linspace(0, 1, 13)
        rates, pulse, tau, z, _ = path_from_variables(a, b, v, lam, geom, 3, mode)
        assert np.allclose(rates[0], a[eng.RATE_COLS].to_numpy(float))
        assert np.allclose(rates[-1], b[eng.RATE_COLS].to_numpy(float))
        assert abs(pulse[0] - float(a["pulse_amplitude_au"])) < 1e-12
        assert abs(pulse[-1] - float(b["pulse_amplitude_au"])) < 1e-12
        assert abs(tau[0] - float(a["glutamate_tau_ms"])) < 1e-12
        assert abs(tau[-1] - float(b["glutamate_tau_ms"])) < 1e-12
        assert path_stretch(z) >= 1.0 - 1e-12
    print("Step28 self-test PASS")


def main() -> None:
    args = parse_args()
    if args.self_test:
        self_test()
        return
    if args.threads < 1:
        raise ValueError("--threads must be >=1")
    set_num_threads(args.threads)

    step25 = discover_step25(args.step25)
    step27 = discover_step27(args.step27)
    inp = load_inputs(step25, step27)
    replay = replay_gate(inp, args.calibration_dt_ms, args.confirm_dt_ms)
    bounds, bounds_source = load_full_broad_bounds(inp, args.allow_accepted_bounds_fallback)
    geom = geometry_stats(inp, bounds)
    comp_map = strict_component_map(inp)

    audit = {
        "script_version": SCRIPT_VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "python": sys.version,
        "platform": platform.platform(),
        "numba_threads": int(get_num_threads()),
        "step25": str(step25),
        "step27": str(step27),
        "step25_status": inp.decision25.get("status"),
        "step27_status": inp.decision27.get("status"),
        "n_confirmed_strong": int(len(inp.confirmed)),
        "n_corrected_accepted_broad": int(len(inp.accepted)),
        "step27_strict_components": [sorted(list(x)) for x in EXPECTED_STRICT_COMPONENTS],
        "strong_threshold": inp.threshold,
        "corrected_broad_cutoff_midpoint": inp.cutoff_mid,
        "corrected_broad_cutoff_strict": inp.cutoff_strict,
        "full_broad_bounds_source": bounds_source,
        "primary_acceptance_cutoff": "strict",
        "no_denominator_exclusion_applied": True,
        "search_path_family": "full-12D endpoint-preserving piecewise-log-linear multi-waypoint",
        "secondary_flex_aux": not args.disable_flex_aux,
        "search_points": args.search_points,
        "verify_points": args.verify_points,
        "max_waypoints": args.max_waypoints,
        "population": args.population,
        "generations": args.generations,
        "restarts": args.restarts,
        **replay,
    }
    if args.preflight_only:
        print(json.dumps(audit, indent=2))
        return

    outdir = args.output.expanduser().resolve()
    if outdir.exists():
        if args.force:
            shutil.rmtree(outdir)
        elif not args.resume:
            raise FileExistsError(f"Output exists: {outdir}. Set STEP28_RESUME=1 to continue or STEP28_FORCE=1 to replace it.")
    outdir.mkdir(parents=True, exist_ok=True)
    chkdir = outdir / "_checkpoints"
    chkdir.mkdir(parents=True, exist_ok=True)
    json_dump(outdir / "00_input_audit.json", audit)

    # Freeze Step27 residual structure and priorities.
    strict27 = inp.components27[inp.components27["criterion"].astype(str).eq("verified_strict")].copy()
    strict27.to_csv(outdir / "01_step27_residual_components.csv", index=False)
    priority = component_priority(inp, comp_map)
    priority.to_csv(outdir / "02_component_pair_priority.csv", index=False)

    ids = sorted(inp.confirmed["sample_id"].astype(int).tolist())
    uf = UnionFind(ids)
    for r in inp.network27.itertuples(index=False):
        if bool(r.strict_cutoff_pass):
            uf.union(int(r.sample_id_a), int(r.sample_id_b))
    if [len(x) for x in uf.groups()] != [36, 2, 1, 1]:
        raise RuntimeError(f"Step27 network replay produced unexpected component sizes: {[len(x) for x in uf.groups()]}")

    lookup = inp.confirmed.set_index("sample_id", drop=False)
    attempts_rows = []
    verified_rows = []
    verified_frames: list[pd.DataFrame] = []
    history_rows = []

    component_pairs_tested = 0
    for pr in priority.itertuples(index=False):
        if component_pairs_tested >= args.max_component_pairs:
            break
        ca, cb = int(pr.component_a), int(pr.component_b)
        members_a = [sid for sid, c in comp_map.items() if c == ca]
        members_b = [sid for sid, c in comp_map.items() if c == cb]
        # If already linked indirectly by a new Step28 bridge, no need to search this original pair.
        if uf.find(members_a[0]) == uf.find(members_b[0]):
            continue
        component_pairs_tested += 1
        epairs = endpoint_pairs_for_components(inp, comp_map, ca, cb, args.pair_attempts)

        pair_succeeded = False
        for ep_idx, ep in epairs.iterrows():
            a, b = int(ep["sample_id_a"]), int(ep["sample_id_b"])
            row_a, row_b = lookup.loc[a], lookup.loc[b]
            stages: list[tuple[str, int]] = [("native_aux", 2)]
            if args.max_waypoints >= 3 and ep_idx < args.deep_pair_attempts:
                stages.append(("native_aux", 3))
            if args.max_waypoints >= 4 and ep_idx == 0:
                stages.append(("native_aux", 4))

            native_best = np.inf
            stage_results: list[tuple[str, int, dict[str, Any], dict[str, Any] | None, pd.DataFrame | None]] = []
            for mode, nwp in stages:
                key = checkpoint_key(a, b, mode, nwp)
                loaded = load_stage_checkpoint(chkdir, key) if args.resume else None
                if loaded is None:
                    result = search_waypoint_path(row_a, row_b, inp, geom, args, nwp, mode, seed_offset=component_pairs_tested * 10000 + ep_idx * 100 + nwp)
                    verify_summary = None
                    frame = None
                    if result["search_joint_barrier_fullgrid"] <= args.verify_trigger:
                        verify_summary, frame = verify_path(
                            row_a, row_b, result["variables"], nwp, mode, inp, geom,
                            args.verify_points, args.grid_dt_ms, args.confirm_dt_ms, args.calibration_dt_ms,
                        )
                    save_stage_checkpoint(chkdir, key, result, verify_summary, frame)
                else:
                    result, verify_summary, frame = loaded
                native_best = min(native_best, float(result["search_joint_barrier_fullgrid"]))
                stage_results.append((mode, nwp, result, verify_summary, frame))
                for h in result.get("history", []):
                    history_rows.append({"sample_id_a": a, "sample_id_b": b, **h})
                attempts_rows.append({
                    "component_a_step27": ca,
                    "component_b_step27": cb,
                    "endpoint_rank": int(ep_idx),
                    "endpoint_distance": float(ep["distance"]),
                    "step27_best_joint_barrier": float(ep["step27_best_joint_barrier"]) if np.isfinite(ep["step27_best_joint_barrier"]) else np.nan,
                    "sample_id_a": a,
                    "sample_id_b": b,
                    "mode": mode,
                    "n_waypoints": nwp,
                    "search_max_score": float(result["search_max_score"]),
                    "search_min_r_targeted_dt0p1": float(result["search_min_r_targeted_dt0p1"]),
                    "search_min_r_fullgrid_dt0p1": float(result["search_min_r_fullgrid_dt0p1"]),
                    "search_joint_barrier_fullgrid": float(result["search_joint_barrier_fullgrid"]),
                    "search_path_stretch": float(result["search_path_stretch"]),
                    "verification_run": verify_summary is not None,
                    "verified_strict": bool(verify_summary and verify_summary["strict_cutoff_pass"]),
                    "verified_midpoint": bool(verify_summary and verify_summary["midpoint_cutoff_pass"]),
                    "verified_max_score": float(verify_summary["max_calibration_score"]) if verify_summary else np.nan,
                    "verified_min_r_dt0p05": float(verify_summary["min_r_dt0p05"]) if verify_summary else np.nan,
                })
                if verify_summary and bool(verify_summary["strict_cutoff_pass"]):
                    verified_rows.append({
                        "sample_id_a": a,
                        "sample_id_b": b,
                        "edge_role": "step28_multiwaypoint_bridge",
                        "mode": mode,
                        "n_waypoints": nwp,
                        "midpoint_cutoff_pass": bool(verify_summary["midpoint_cutoff_pass"]),
                        "strict_cutoff_pass": True,
                        "max_calibration_score": float(verify_summary["max_calibration_score"]),
                        "min_r_dt0p05": float(verify_summary["min_r_dt0p05"]),
                        "min_control_plateau_at_support_node_dt0p05": float(verify_summary["min_control_plateau_at_support_node_dt0p05"]),
                        "path_stretch_standardized_log_rate": float(verify_summary["path_stretch_standardized_log_rate"]),
                    })
                    if frame is not None:
                        verified_frames.append(frame.assign(path_role="step28_multiwaypoint_bridge"))
                    uf.union(a, b)
                    pair_succeeded = True
                    break
            if pair_succeeded:
                break

            # Secondary flex-aux only after native failure and only when native search approached feasibility.
            if not args.disable_flex_aux and (native_best <= args.flex_trigger or ep_idx < args.deep_pair_attempts):
                flex_stages: list[tuple[str, int]] = [("flex_aux", 2)]
                if args.max_waypoints >= 3 and ep_idx < args.deep_pair_attempts:
                    flex_stages.append(("flex_aux", 3))
                if args.max_waypoints >= 4 and ep_idx == 0:
                    flex_stages.append(("flex_aux", 4))
                for mode, nwp in flex_stages:
                    key = checkpoint_key(a, b, mode, nwp)
                    loaded = load_stage_checkpoint(chkdir, key) if args.resume else None
                    if loaded is None:
                        result = search_waypoint_path(row_a, row_b, inp, geom, args, nwp, mode, seed_offset=500000 + component_pairs_tested * 10000 + ep_idx * 100 + nwp)
                        verify_summary = None
                        frame = None
                        if result["search_joint_barrier_fullgrid"] <= args.verify_trigger:
                            verify_summary, frame = verify_path(
                                row_a, row_b, result["variables"], nwp, mode, inp, geom,
                                args.verify_points, args.grid_dt_ms, args.confirm_dt_ms, args.calibration_dt_ms,
                            )
                        save_stage_checkpoint(chkdir, key, result, verify_summary, frame)
                    else:
                        result, verify_summary, frame = loaded
                    for h in result.get("history", []):
                        history_rows.append({"sample_id_a": a, "sample_id_b": b, **h})
                    attempts_rows.append({
                        "component_a_step27": ca,
                        "component_b_step27": cb,
                        "endpoint_rank": int(ep_idx),
                        "endpoint_distance": float(ep["distance"]),
                        "step27_best_joint_barrier": float(ep["step27_best_joint_barrier"]) if np.isfinite(ep["step27_best_joint_barrier"]) else np.nan,
                        "sample_id_a": a,
                        "sample_id_b": b,
                        "mode": mode,
                        "n_waypoints": nwp,
                        "search_max_score": float(result["search_max_score"]),
                        "search_min_r_targeted_dt0p1": float(result["search_min_r_targeted_dt0p1"]),
                        "search_min_r_fullgrid_dt0p1": float(result["search_min_r_fullgrid_dt0p1"]),
                        "search_joint_barrier_fullgrid": float(result["search_joint_barrier_fullgrid"]),
                        "search_path_stretch": float(result["search_path_stretch"]),
                        "verification_run": verify_summary is not None,
                        "verified_strict": bool(verify_summary and verify_summary["strict_cutoff_pass"]),
                        "verified_midpoint": bool(verify_summary and verify_summary["midpoint_cutoff_pass"]),
                        "verified_max_score": float(verify_summary["max_calibration_score"]) if verify_summary else np.nan,
                        "verified_min_r_dt0p05": float(verify_summary["min_r_dt0p05"]) if verify_summary else np.nan,
                    })
                    if verify_summary and bool(verify_summary["strict_cutoff_pass"]):
                        verified_rows.append({
                            "sample_id_a": a,
                            "sample_id_b": b,
                            "edge_role": "step28_multiwaypoint_bridge",
                            "mode": mode,
                            "n_waypoints": nwp,
                            "midpoint_cutoff_pass": bool(verify_summary["midpoint_cutoff_pass"]),
                            "strict_cutoff_pass": True,
                            "max_calibration_score": float(verify_summary["max_calibration_score"]),
                            "min_r_dt0p05": float(verify_summary["min_r_dt0p05"]),
                            "min_control_plateau_at_support_node_dt0p05": float(verify_summary["min_control_plateau_at_support_node_dt0p05"]),
                            "path_stretch_standardized_log_rate": float(verify_summary["path_stretch_standardized_log_rate"]),
                        })
                        if frame is not None:
                            verified_frames.append(frame.assign(path_role="step28_multiwaypoint_bridge"))
                        uf.union(a, b)
                        pair_succeeded = True
                        break
            if pair_succeeded:
                break

        if len(uf.groups()) == 1:
            break

    attempts = pd.DataFrame(attempts_rows)
    attempts.to_csv(outdir / "03_search_attempts.csv", index=False)
    if history_rows:
        pd.DataFrame(history_rows).to_csv(outdir / "04_search_history.csv.gz", index=False, compression="gzip")
    else:
        pd.DataFrame().to_csv(outdir / "04_search_history.csv.gz", index=False, compression="gzip")
    verified = pd.DataFrame(verified_rows)
    verified.to_csv(outdir / "05_verified_bridges.csv", index=False)
    all_path_points = pd.concat(verified_frames, ignore_index=True) if verified_frames else pd.DataFrame()
    all_path_points.to_csv(outdir / "06_verified_path_points.csv.gz", index=False, compression="gzip")

    existing = inp.network27.copy()
    if "n_waypoints" not in existing.columns:
        existing["n_waypoints"] = 0
    final_edges = pd.concat([existing, verified], ignore_index=True, sort=False) if len(verified) else existing.copy()
    final_edges.to_csv(outdir / "07_final_network_edges.csv", index=False)
    comp_mid = components_table(ids, final_edges, "midpoint_cutoff_pass", "step28_midpoint")
    comp_strict = components_table(ids, final_edges, "strict_cutoff_pass", "step28_strict")
    native_edges = final_edges[final_edges["mode"].astype(str).eq("native_aux")].copy()
    comp_native_strict = components_table(ids, native_edges, "strict_cutoff_pass", "step28_native_only_strict")
    final_components = pd.concat([comp_mid, comp_strict, comp_native_strict], ignore_index=True)
    final_components.to_csv(outdir / "08_final_components.csv", index=False)

    # Barrier summary by original Step27 component pair.
    barrier_rows = []
    if len(attempts):
        for (ca, cb), g in attempts.groupby(["component_a_step27", "component_b_step27"]):
            best = g.sort_values("search_joint_barrier_fullgrid").iloc[0]
            barrier_rows.append({
                "component_a_step27": int(ca),
                "component_b_step27": int(cb),
                "best_sample_id_a": int(best["sample_id_a"]),
                "best_sample_id_b": int(best["sample_id_b"]),
                "best_mode": str(best["mode"]),
                "best_n_waypoints": int(best["n_waypoints"]),
                "best_search_joint_barrier_fullgrid": float(best["search_joint_barrier_fullgrid"]),
                "best_search_max_score": float(best["search_max_score"]),
                "best_search_min_r_fullgrid_dt0p1": float(best["search_min_r_fullgrid_dt0p1"]),
                "any_verified_strict_bridge": bool(g["verified_strict"].astype(bool).any()),
            })
    barriers = pd.DataFrame(barrier_rows)
    barriers.to_csv(outdir / "09_residual_barrier_summary.csv", index=False)

    strict_groups = comp_strict.groupby("component_id")["sample_id"].apply(list).tolist()
    strict_sizes = sorted([len(x) for x in strict_groups], reverse=True)
    native_groups = comp_native_strict.groupby("component_id")["sample_id"].apply(list).tolist()
    native_strict_sizes = sorted([len(x) for x in native_groups], reverse=True)
    native_new = verified[verified["mode"].eq("native_aux")] if len(verified) else pd.DataFrame()
    flex_new = verified[verified["mode"].eq("flex_aux")] if len(verified) else pd.DataFrame()
    if strict_sizes == [40]:
        if native_strict_sizes == [40]:
            status = "FULL12D_MULTIWAYPOINT_CONNECTIVITY_ESTABLISHED_NATIVE_AUX"
        else:
            status = "FULL14D_MULTIWAYPOINT_CONNECTIVITY_ESTABLISHED_WITH_FLEX_AUX"
    elif len(verified) > 0:
        status = "PARTIAL_RESIDUAL_CONNECTIVITY_AFTER_FULLDIM_STRESS_TEST"
    else:
        status = "RESIDUAL_CONNECTIVITY_NOT_ESTABLISHED_AFTER_FULLDIM_STRESS_TEST"

    decision = {
        "status": status,
        "interpretation_scope": (
            "Finite-budget connectivity stress test of the three residual Step27 blocks using endpoint-preserving "
            "piecewise-log-linear multi-waypoint paths in the full 12-dimensional standardized log-rate space; "
            "secondary flex_aux paths also optimize the two auxiliary calibration coordinates. Failure to connect "
            "is not a proof of global disconnection."
        ),
        "n_strong_candidates": 40,
        "step27_strict_component_sizes": [36, 2, 1, 1],
        "step28_strict_component_sizes": strict_sizes,
        "step28_native_only_strict_component_sizes": native_strict_sizes,
        "component_pairs_tested": int(component_pairs_tested),
        "search_attempts": int(len(attempts)),
        "new_verified_bridges": int(len(verified)),
        "new_verified_native_bridges": int(len(native_new)),
        "new_verified_flex_aux_bridges": int(len(flex_new)),
        "strong_threshold": inp.threshold,
        "midpoint_cutoff": inp.cutoff_mid,
        "strict_cutoff": inp.cutoff_strict,
        "primary_acceptance_cutoff": "strict",
        "search_waypoint_levels": list(range(2, args.max_waypoints + 1)),
        "full_broad_bounds_source": bounds_source,
        "no_denominator_exclusion_applied": True,
        "engine_replay_gate_pass": True,
    }
    json_dump(outdir / "10_scientific_decision.json", decision)

    make_figures(inp, final_edges, final_components, all_path_points, attempts, outdir)

    readme = f"""# Step 28 results - residual full-dimensional connectivity stress test\n\nStatus: **{status}**\n\nStep27 strict component sizes entering Step28: **36 + 2 + 1 + 1**.\nStep28 strict component sizes: **{' + '.join(map(str, strict_sizes))}**.\nNew verified bridges: **{len(verified)}** (native={len(native_new)}, flex_aux={len(flex_new)}).\n\nPrimary paths use freely optimized internal waypoints in the full 12-dimensional standardized log-rate space.\nSecondary flex_aux paths also optimize the two auxiliary calibration coordinates while preserving both endpoints.\nEvery accepted bridge was verified at {args.verify_points} path points on the full frozen 80-node response grid with dt=0.05-ms confirmation.\nThe strict corrected calibration cutoff ({inp.cutoff_strict:.12g}) is the primary acceptance criterion.\nNo denominator exclusion was introduced.\n\nFailure to obtain one connected network within this finite multi-waypoint search budget is not proof of global topological disconnection.\n"""
    (outdir / "README_RESULTS.md").write_text(readme, encoding="utf-8")

    run_summary = {
        **decision,
        "output_directory": str(outdir),
        "numba_threads": int(get_num_threads()),
        "files": sorted([p.name for p in outdir.iterdir() if p.is_file()]),
    }
    json_dump(outdir / "run_summary.json", run_summary)
    print(json.dumps(decision, indent=2))


if __name__ == "__main__":
    main()
