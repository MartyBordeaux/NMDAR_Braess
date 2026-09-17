#!/usr/bin/env python3
"""NMDA Braess Step 27 v1.0 - curved bridge optimization.

Step 26 showed that direct log-linear paths do not establish global connectivity
of the corrected strong/control-compatible ensemble. Step 27 tests the next
scientific question: can the apparently disconnected sampled components be
linked by curved continuous paths in parameter space?

Primary path family:
  - endpoints are fixed at two Step-25 confirmed strong vectors;
  - 12 kinetic rates are represented in standardized log-rate coordinates;
  - the straight endpoint interpolation is bent with a low-dimensional sine
    basis aligned to the leading PCA directions of the corrected accepted broad
    ensemble;
  - historical calibration-drive coordinates remain on their endpoint-preserving
    geometric interpolation ("native_aux").

Secondary path family:
  - the same kinetic curve is used;
  - pulse_amplitude_au and glutamate_tau_ms may also bend smoothly, with fixed
    endpoint values ("flex_aux").

A path is accepted only after a 41-point final verification in which every path
point is control-compatible at dt=0.05 ms and is strong on the frozen 80-node
response grid. Response maxima are first found at dt=0.1 ms and reconfirmed at
0.05 ms, with a full dt=0.05-ms grid fallback for any failing point.

No denominator cutoff is introduced. Small control plateaus are retained as a
reported diagnostic.
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
EXPECTED_STEP26_STATUS = "SAMPLED_CONNECTIVITY_NOT_ESTABLISHED"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Step 27 curved strong-compatible bridge optimization")
    p.add_argument("--step25", default="auto")
    p.add_argument("--step26", default="auto")
    p.add_argument(
        "--output",
        type=Path,
        default=Path("/root/nmda2/step_27/results_step_27_curved_bridge_optimization"),
    )
    p.add_argument("--threads", type=int, default=2)
    p.add_argument("--search-points", type=int, default=13)
    p.add_argument("--verify-points", type=int, default=41)
    p.add_argument("--n-pcs", type=int, default=5)
    p.add_argument("--n-modes", type=int, default=2)
    p.add_argument("--population", type=int, default=24)
    p.add_argument("--generations", type=int, default=18)
    p.add_argument("--joint-generations", type=int, default=10)
    p.add_argument("--max-component-pairs", type=int, default=12)
    p.add_argument("--pair-attempts", type=int, default=3)
    p.add_argument("--coeff-bound", type=float, default=1.5)
    p.add_argument("--seed", type=int, default=20260914)
    p.add_argument("--grid-dt-ms", type=float, default=0.1)
    p.add_argument("--confirm-dt-ms", type=float, default=0.05)
    p.add_argument("--calibration-dt-ms", type=float, default=0.05)
    p.add_argument("--allow-accepted-bounds-fallback", action="store_true")
    p.add_argument("--disable-flex-aux", action="store_true")
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
            p = marker.parent
            if required.issubset({x.name for x in p.iterdir()}):
                return p.resolve()
    raise FileNotFoundError("Could not auto-discover Step-25 results")


def discover_step26(spec: str) -> Path:
    if spec != "auto":
        p = Path(spec).expanduser().resolve()
        if not p.is_dir():
            raise FileNotFoundError(p)
        return p
    candidates = [
        Path("/root/nmda2/step_26/results_step_26_strong_connectivity"),
        Path("/root/nmda2/results_step_26_strong_connectivity"),
        Path("/root/nmda2/NMDAR_Braess_step26_strong_connectivity_v1_0/results_step_26_strong_connectivity"),
    ]
    required = {"03_coarse_edge_summary.csv", "05_coarse_components.csv", "11_scientific_decision.json"}
    for p in candidates:
        if p.is_dir() and required.issubset({x.name for x in p.iterdir()}):
            return p.resolve()
    root = Path("/root/nmda2")
    if root.exists():
        for marker in root.rglob("11_scientific_decision.json"):
            p = marker.parent
            if required.issubset({x.name for x in p.iterdir()}):
                return p.resolve()
    raise FileNotFoundError("Could not auto-discover Step-26 results")


@dataclass
class Inputs:
    step25: Path
    step26: Path
    audit25: dict[str, Any]
    decision25: dict[str, Any]
    decision26: dict[str, Any]
    targets: dict[str, Any]
    confirmed: pd.DataFrame
    accepted: pd.DataFrame
    geometry: pd.DataFrame
    pairwise: pd.DataFrame
    coarse_edges: pd.DataFrame
    coarse_components: pd.DataFrame
    refined_edges26: pd.DataFrame
    threshold: float
    cutoff_mid: float
    cutoff_strict: float


def load_inputs(step25: Path, step26: Path) -> Inputs:
    req25 = [
        "00_input_audit.json",
        "01_corrected_target_reconstruction.json",
        "02_corrected_accepted_cohorts.csv.gz",
        "05_strong_candidates_confirmed.csv",
        "08_geometry_coordinates.csv.gz",
        "09_geometry_pairwise_distances.csv",
        "15_scientific_decision.json",
    ]
    req26 = [
        "03_coarse_edge_summary.csv",
        "05_coarse_components.csv",
        "06_refined_edge_summary.csv",
        "11_scientific_decision.json",
    ]
    for name in req25:
        if not (step25 / name).exists():
            raise FileNotFoundError(step25 / name)
    for name in req26:
        if not (step26 / name).exists():
            raise FileNotFoundError(step26 / name)

    audit25 = json.loads((step25 / "00_input_audit.json").read_text(encoding="utf-8"))
    decision25 = json.loads((step25 / "15_scientific_decision.json").read_text(encoding="utf-8"))
    decision26 = json.loads((step26 / "11_scientific_decision.json").read_text(encoding="utf-8"))
    targets = json.loads((step25 / "01_corrected_target_reconstruction.json").read_text(encoding="utf-8"))
    confirmed = pd.read_csv(step25 / "05_strong_candidates_confirmed.csv")
    accepted = pd.read_csv(step25 / "02_corrected_accepted_cohorts.csv.gz")
    accepted = accepted[accepted["prior"].astype(str).str.lower().eq("broad")].copy()
    geometry = pd.read_csv(step25 / "08_geometry_coordinates.csv.gz")
    pairwise = pd.read_csv(step25 / "09_geometry_pairwise_distances.csv")
    coarse_edges = pd.read_csv(step26 / "03_coarse_edge_summary.csv")
    coarse_components = pd.read_csv(step26 / "05_coarse_components.csv")
    refined_edges26 = pd.read_csv(step26 / "06_refined_edge_summary.csv")

    if decision25.get("status") != EXPECTED_STEP25_STATUS:
        raise RuntimeError(f"Unexpected Step25 status: {decision25.get('status')!r}")
    if decision26.get("status") != EXPECTED_STEP26_STATUS:
        raise RuntimeError(f"Unexpected Step26 status: {decision26.get('status')!r}")
    if len(confirmed) != 40 or not confirmed["strong_confirmed_dt0p05"].astype(bool).all():
        raise RuntimeError("Step25 confirmed strong table is not the expected 40/40 cohort")
    if len(accepted) != 5000:
        raise RuntimeError(f"Expected 5000 corrected accepted broad rows; got {len(accepted)}")

    threshold = float(audit25["strong_threshold"])
    cutoff_mid = float(audit25["corrected_broad_cutoff_midpoint"])
    cutoff_strict = float(audit25["corrected_broad_max_accepted_score"])
    if abs(float(decision26["strong_threshold"]) - threshold) > 1e-12:
        raise RuntimeError("Step25/26 strong-threshold mismatch")
    if abs(float(decision26["midpoint_cutoff"]) - cutoff_mid) > 1e-12:
        raise RuntimeError("Step25/26 midpoint-cutoff mismatch")

    return Inputs(
        step25=step25,
        step26=step26,
        audit25=audit25,
        decision25=decision25,
        decision26=decision26,
        targets=targets,
        confirmed=confirmed,
        accepted=accepted,
        geometry=geometry,
        pairwise=pairwise,
        coarse_edges=coarse_edges,
        coarse_components=coarse_components,
        refined_edges26=refined_edges26,
        threshold=threshold,
        cutoff_mid=cutoff_mid,
        cutoff_strict=cutoff_strict,
    )


def replay_gate(inp: Inputs, calib_dt: float, confirm_dt: float) -> dict[str, Any]:
    ids = [1900, 10052, 27084, 37318, 49029]
    ids = [x for x in ids if x in set(inp.confirmed["sample_id"].astype(int))]
    sub = inp.confirmed.set_index("sample_id").loc[ids].reset_index()
    rates = sub[eng.RATE_COLS].to_numpy(float)
    pulse = sub["pulse_amplitude_au"].to_numpy(float)
    tau = sub["glutamate_tau_ms"].to_numpy(float)
    met = eng.calibration_ensemble(rates, pulse, tau, calib_dt)
    score = eng.corrected_score_from_metrics(met[:, 1], met[:, 2], met[:, 3], met[:, 5], inp.targets)
    score_err = np.abs(score - sub["calibration_score_corrected"].to_numpy(float))
    r_err = []
    for i, row in sub.iterrows():
        out, _ = eng.run_base_condition(rates[i : i + 1], float(row["G_peak"]), float(row["tau_G_ms"]), confirm_dt)
        r_err.append(abs(float(out[0, 2]) - float(row["r_plateau_dt0p05"])))
    r_err = np.asarray(r_err)
    passed = bool(np.max(score_err) < 1e-8 and np.max(r_err) < 5e-4)
    return {
        "sample_ids": ids,
        "max_abs_calibration_score_error": float(np.max(score_err)),
        "max_abs_response_ratio_error": float(np.max(r_err)),
        "pass": passed,
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


def geometry_basis(inp: Inputs, n_pcs: int) -> dict[str, np.ndarray]:
    logs = np.log(inp.accepted[eng.RATE_COLS].to_numpy(float))
    mean = logs.mean(axis=0)
    sd = logs.std(axis=0, ddof=1)
    if np.any(sd <= 0):
        raise RuntimeError("Zero variance in accepted-broad log-rate coordinates")
    z = (logs - mean) / sd
    _, _, vt = np.linalg.svd(z, full_matrices=False)
    n_pcs = min(n_pcs, vt.shape[0])
    aux_logs = np.log(inp.accepted[eng.AUX_COLS].to_numpy(float))
    aux_mean = aux_logs.mean(axis=0)
    aux_sd = aux_logs.std(axis=0, ddof=1)
    if np.any(aux_sd <= 0):
        raise RuntimeError("Zero variance in accepted-broad log-aux coordinates")
    return {
        "rate_log_mean": mean,
        "rate_log_sd": sd,
        "pc_basis": vt[:n_pcs],
        "aux_log_mean": aux_mean,
        "aux_log_sd": aux_sd,
    }


def target_forcing_nodes(inp: Inputs) -> list[tuple[float, float]]:
    nodes = list(dict.fromkeys(
        (float(g), float(t))
        for g, t in inp.confirmed[["G_peak", "tau_G_ms"]].itertuples(index=False, name=None)
    ))
    return nodes


class UnionFind:
    def __init__(self, nodes: list[int]):
        self.parent = {int(x): int(x) for x in nodes}
        self.rank = {int(x): 0 for x in nodes}

    def find(self, x: int) -> int:
        x = int(x)
        p = self.parent[x]
        if p != x:
            self.parent[x] = self.find(p)
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

    def components(self) -> dict[int, list[int]]:
        out: dict[int, list[int]] = {}
        for x in self.parent:
            out.setdefault(self.find(x), []).append(x)
        return out


def native_component_map(inp: Inputs) -> dict[int, int]:
    sub = inp.coarse_components[inp.coarse_components["criterion"].astype(str).eq("native_joint_mid_dt0p1")].copy()
    if len(sub) != len(inp.confirmed):
        raise RuntimeError("Step26 native coarse component table is incomplete")
    return {int(r.sample_id): int(r.component_id) for r in sub.itertuples(index=False)}


def component_priority_table(inp: Inputs, comp_map: dict[int, int]) -> pd.DataFrame:
    rows = []
    pw = inp.pairwise.copy()
    for r in pw.itertuples(index=False):
        a, b = int(r.sample_id_a), int(r.sample_id_b)
        ca, cb = comp_map[a], comp_map[b]
        if ca == cb:
            continue
        c1, c2 = sorted((ca, cb))
        rows.append({
            "component_a": c1,
            "component_b": c2,
            "sample_id_a": a,
            "sample_id_b": b,
            "distance": float(r.standardized_log_rate_distance),
        })
    df = pd.DataFrame(rows)
    out = []
    for (ca, cb), g in df.groupby(["component_a", "component_b"], sort=False):
        g = g.sort_values(["distance", "sample_id_a", "sample_id_b"]).reset_index(drop=True)
        out.append({
            "component_a": int(ca),
            "component_b": int(cb),
            "min_distance": float(g.loc[0, "distance"]),
            "nearest_sample_id_a": int(g.loc[0, "sample_id_a"]),
            "nearest_sample_id_b": int(g.loc[0, "sample_id_b"]),
            "n_endpoint_pairs": int(len(g)),
        })
    return pd.DataFrame(out).sort_values(["min_distance", "component_a", "component_b"]).reset_index(drop=True)


def endpoint_pairs_for_components(inp: Inputs, comp_map: dict[int, int], ca: int, cb: int, n: int) -> pd.DataFrame:
    rows = []
    for r in inp.pairwise.itertuples(index=False):
        a, b = int(r.sample_id_a), int(r.sample_id_b)
        if {comp_map[a], comp_map[b]} == {int(ca), int(cb)}:
            rows.append({
                "sample_id_a": a,
                "sample_id_b": b,
                "distance": float(r.standardized_log_rate_distance),
            })
    return pd.DataFrame(rows).sort_values(["distance", "sample_id_a", "sample_id_b"]).head(n).reset_index(drop=True)


def path_from_coeffs(
    row_a: pd.Series,
    row_b: pd.Series,
    coeffs: np.ndarray,
    lam: np.ndarray,
    geom: dict[str, np.ndarray],
    n_modes: int,
    mode: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    n_pcs = geom["pc_basis"].shape[0]
    rate_dim = n_modes * n_pcs
    if mode == "native_aux":
        expected = rate_dim
    elif mode == "flex_aux":
        expected = rate_dim + n_modes * 2
    else:
        raise ValueError(mode)
    if len(coeffs) != expected:
        raise ValueError(f"Expected {expected} coefficients, got {len(coeffs)}")

    la = np.log(row_a[eng.RATE_COLS].to_numpy(float))
    lb = np.log(row_b[eng.RATE_COLS].to_numpy(float))
    z0 = (la - geom["rate_log_mean"]) / geom["rate_log_sd"]
    z1 = (lb - geom["rate_log_mean"]) / geom["rate_log_sd"]
    z = (1.0 - lam[:, None]) * z0[None, :] + lam[:, None] * z1[None, :]

    cr = coeffs[:rate_dim].reshape(n_modes, n_pcs)
    for m in range(1, n_modes + 1):
        shape = np.sin(math.pi * m * lam)[:, None]
        direction = cr[m - 1] @ geom["pc_basis"]
        z += shape * direction[None, :]
    rate_log = geom["rate_log_mean"][None, :] + z * geom["rate_log_sd"][None, :]
    rates = np.exp(rate_log)

    aux_a = np.log(row_a[eng.AUX_COLS].to_numpy(float))
    aux_b = np.log(row_b[eng.AUX_COLS].to_numpy(float))
    aux_log = (1.0 - lam[:, None]) * aux_a[None, :] + lam[:, None] * aux_b[None, :]
    if mode == "flex_aux":
        ca = coeffs[rate_dim:].reshape(n_modes, 2)
        az = (aux_log - geom["aux_log_mean"][None, :]) / geom["aux_log_sd"][None, :]
        for m in range(1, n_modes + 1):
            az += np.sin(math.pi * m * lam)[:, None] * ca[m - 1][None, :]
        aux_log = geom["aux_log_mean"][None, :] + az * geom["aux_log_sd"][None, :]
    aux = np.exp(aux_log)
    return rates, aux[:, 0], aux[:, 1], rate_log


def within_bounds(rate_log: np.ndarray, pulse: np.ndarray, tau: np.ndarray, bounds: dict[str, np.ndarray]) -> bool:
    aux_log = np.log(np.column_stack([pulse, tau]))
    return bool(
        np.all(rate_log >= bounds["rate_log_min"][None, :] - 1e-12)
        and np.all(rate_log <= bounds["rate_log_max"][None, :] + 1e-12)
        and np.all(aux_log >= bounds["aux_log_min"][None, :] - 1e-12)
        and np.all(aux_log <= bounds["aux_log_max"][None, :] + 1e-12)
    )


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


def full_grid_nodes() -> list[tuple[float, float]]:
    return [(float(g), float(t)) for g in eng.G_PEAK_GRID for t in eng.TAU_GRID]


def confirm_best_nodes(rates: np.ndarray, coarse: dict[str, np.ndarray], dt: float, threshold: float) -> dict[str, np.ndarray]:
    n = len(rates)
    result = {k: np.array(v, copy=True) for k, v in coarse.items()}
    r_confirm = np.full(n, np.nan, dtype=float)
    c_confirm = np.full(n, np.nan, dtype=float)
    b_confirm = np.full(n, np.nan, dtype=float)
    minc = np.full(n, np.nan, dtype=float)
    minb = np.full(n, np.nan, dtype=float)
    unique = sorted({(float(coarse["G_peak"][i]), float(coarse["tau_G_ms"][i])) for i in range(n) if np.isfinite(coarse["G_peak"][i])})
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


def verify_path(
    row_a: pd.Series,
    row_b: pd.Series,
    coeffs: np.ndarray,
    mode: str,
    inp: Inputs,
    geom: dict[str, np.ndarray],
    bounds: dict[str, np.ndarray],
    n_modes: int,
    n_points: int,
    grid_dt: float,
    confirm_dt: float,
    calib_dt: float,
) -> tuple[dict[str, Any], pd.DataFrame]:
    lam = np.linspace(0.0, 1.0, n_points)
    rates, pulse, tau, rate_log = path_from_coeffs(row_a, row_b, coeffs, lam, geom, n_modes, mode)
    in_bounds = within_bounds(rate_log, pulse, tau, bounds)
    calib = calibration_eval(rates, pulse, tau, inp, calib_dt)
    coarse = response_scan_nodes(rates, full_grid_nodes(), grid_dt)
    confirmed = confirm_best_nodes(rates, coarse, confirm_dt, inp.threshold)

    strong = np.isfinite(confirmed["r"]) & (confirmed["r"] >= inp.threshold)
    compat_mid = np.isfinite(calib["score"]) & (calib["score"] <= inp.cutoff_mid)
    compat_strict = np.isfinite(calib["score"]) & (calib["score"] <= inp.cutoff_strict)
    pass_mid = bool(in_bounds and strong.all() and compat_mid.all())
    pass_strict = bool(in_bounds and strong.all() and compat_strict.all())

    sid_a, sid_b = int(row_a["sample_id"]), int(row_b["sample_id"])
    summary = {
        "sample_id_a": sid_a,
        "sample_id_b": sid_b,
        "mode": mode,
        "n_points": int(n_points),
        "within_full_broad_bounds": bool(in_bounds),
        "max_calibration_score": float(np.nanmax(calib["score"])),
        "min_r_dt0p05": float(np.nanmin(confirmed["r"])),
        "lambda_at_max_score": float(lam[int(np.nanargmax(calib["score"]))]),
        "lambda_at_min_r": float(lam[int(np.nanargmin(confirmed["r"]))]),
        "min_control_plateau_at_support_node_dt0p05": float(np.nanmin(confirmed["control"])),
        "midpoint_cutoff_pass": pass_mid,
        "strict_cutoff_pass": pass_strict,
        "strong_pass": bool(strong.all()),
    }
    frame = pd.DataFrame({
        "sample_id_a": sid_a,
        "sample_id_b": sid_b,
        "mode": mode,
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


def batch_calibration_objective(
    proposals: np.ndarray,
    row_a: pd.Series,
    row_b: pd.Series,
    lam: np.ndarray,
    mode: str,
    inp: Inputs,
    geom: dict[str, np.ndarray],
    bounds: dict[str, np.ndarray],
    n_modes: int,
    calib_dt: float,
    curvature_penalty: float = 0.002,
) -> tuple[np.ndarray, list[dict[str, Any]]]:
    pop = len(proposals)
    npts = len(lam)
    objective = np.full(pop, 1e6, dtype=float)
    cache: list[dict[str, Any]] = [{} for _ in range(pop)]
    valid_idx = []
    rates_list = []
    pulse_list = []
    tau_list = []
    for i, c in enumerate(proposals):
        rates, pulse, tau, rate_log = path_from_coeffs(row_a, row_b, c, lam, geom, n_modes, mode)
        ok = within_bounds(rate_log, pulse, tau, bounds)
        cache[i] = {"rates": rates, "pulse": pulse, "tau": tau, "in_bounds": ok}
        if ok:
            valid_idx.append(i)
            rates_list.append(rates)
            pulse_list.append(pulse)
            tau_list.append(tau)
    if valid_idx:
        rates_all = np.concatenate(rates_list, axis=0)
        pulse_all = np.concatenate(pulse_list)
        tau_all = np.concatenate(tau_list)
        cal = calibration_eval(rates_all, pulse_all, tau_all, inp, calib_dt)
        scores = cal["score"].reshape(len(valid_idx), npts)
        for j, i in enumerate(valid_idx):
            barrier = float(np.nanmax(scores[j]) / inp.cutoff_mid)
            objective[i] = barrier + curvature_penalty * float(np.mean(proposals[i] ** 2))
            cache[i]["max_score"] = float(np.nanmax(scores[j]))
            cache[i]["score_profile"] = scores[j].copy()
    return objective, cache


def targeted_min_r(proposals_cache: list[dict[str, Any]], indices: list[int], nodes: list[tuple[float, float]], dt: float) -> dict[int, float]:
    if not indices:
        return {}
    npts = len(proposals_cache[indices[0]]["rates"])
    rates_all = np.concatenate([proposals_cache[i]["rates"] for i in indices], axis=0)
    scan = response_scan_nodes(rates_all, nodes, dt)
    r = scan["r"].reshape(len(indices), npts)
    return {idx: float(np.nanmin(r[j])) for j, idx in enumerate(indices)}


def evolutionary_search(
    row_a: pd.Series,
    row_b: pd.Series,
    mode: str,
    inp: Inputs,
    geom: dict[str, np.ndarray],
    bounds: dict[str, np.ndarray],
    target_nodes: list[tuple[float, float]],
    args: argparse.Namespace,
    seed: int,
    initial: np.ndarray | None = None,
) -> dict[str, Any]:
    n_pcs = geom["pc_basis"].shape[0]
    dim = args.n_modes * n_pcs + (0 if mode == "native_aux" else args.n_modes * 2)
    rng = np.random.default_rng(seed)
    mean = np.zeros(dim, dtype=float) if initial is None else np.asarray(initial, dtype=float).copy()
    sigma = np.full(dim, 0.45 if initial is None else 0.20, dtype=float)
    lam = np.linspace(0.0, 1.0, args.search_points)
    best_coeff = mean.copy()
    best_obj = float("inf")
    best_max_score = float("inf")
    history = []

    for gen in range(args.generations):
        pop = rng.normal(mean, sigma, size=(args.population, dim))
        pop = np.clip(pop, -args.coeff_bound, args.coeff_bound)
        pop[0] = mean
        if gen == 0:
            pop[1] = 0.0
        obj, cache = batch_calibration_objective(
            pop, row_a, row_b, lam, mode, inp, geom, bounds, args.n_modes, args.calibration_dt_ms
        )
        order = np.argsort(obj)
        if obj[order[0]] < best_obj:
            best_obj = float(obj[order[0]])
            best_coeff = pop[order[0]].copy()
            best_max_score = float(cache[order[0]].get("max_score", np.inf))
        elite_n = max(4, min(args.population // 4, 8))
        elite = pop[order[:elite_n]]
        weights = np.linspace(2.0, 1.0, elite_n)
        weights /= weights.sum()
        mean = np.sum(elite * weights[:, None], axis=0)
        elite_sd = elite.std(axis=0, ddof=0)
        sigma = np.maximum(0.06, 0.65 * sigma + 0.35 * elite_sd)
        sigma *= 0.94
        history.append({"generation": gen, "best_calibration_objective": best_obj, "best_max_score": best_max_score})

    # Check true frozen-grid strong support at the search nodes for the calibration-optimal curve.
    rates, pulse, tau, rate_log = path_from_coeffs(row_a, row_b, best_coeff, lam, geom, args.n_modes, mode)
    cal = calibration_eval(rates, pulse, tau, inp, args.calibration_dt_ms)
    response = response_scan_nodes(rates, full_grid_nodes(), args.grid_dt_ms)
    min_r = float(np.nanmin(response["r"]))
    max_score = float(np.nanmax(cal["score"]))
    best_joint = max(max_score / inp.cutoff_mid, inp.threshold / max(min_r, 1e-300))

    # If calibration succeeds but strong support fails, run a smaller joint rescue search.
    if args.joint_generations > 0 and (max_score <= inp.cutoff_mid * 1.25) and (min_r < inp.threshold):
        mean = best_coeff.copy()
        sigma = np.full(dim, 0.18, dtype=float)
        for gen in range(args.joint_generations):
            pop_n = max(12, args.population // 2)
            pop = rng.normal(mean, sigma, size=(pop_n, dim))
            pop = np.clip(pop, -args.coeff_bound, args.coeff_bound)
            pop[0] = mean
            cal_obj, cache = batch_calibration_objective(
                pop, row_a, row_b, lam, mode, inp, geom, bounds, args.n_modes, args.calibration_dt_ms,
                curvature_penalty=0.001,
            )
            order_cal = np.argsort(cal_obj)
            resp_idx = [int(i) for i in order_cal[: min(6, len(order_cal))] if cal_obj[i] < 3.0]
            minr = targeted_min_r(cache, resp_idx, target_nodes, args.grid_dt_ms)
            joint = np.full(pop_n, 1e6, dtype=float)
            for i in resp_idx:
                score_bar = float(cache[i].get("max_score", np.inf) / inp.cutoff_mid)
                rbar = inp.threshold / max(float(minr.get(i, 0.0)), 1e-300)
                joint[i] = max(score_bar, rbar) + 0.001 * float(np.mean(pop[i] ** 2))
            order = np.argsort(joint)
            if joint[order[0]] < best_joint:
                best_joint = float(joint[order[0]])
                best_coeff = pop[order[0]].copy()
            elite_idx = [int(i) for i in order[: max(3, min(5, pop_n))] if np.isfinite(joint[i]) and joint[i] < 1e5]
            if not elite_idx:
                break
            elite = pop[elite_idx]
            mean = elite.mean(axis=0)
            sigma = np.maximum(0.05, 0.70 * sigma + 0.30 * elite.std(axis=0, ddof=0))
            sigma *= 0.92
            history.append({"generation": args.generations + gen, "best_joint_objective": best_joint})

        rates, pulse, tau, rate_log = path_from_coeffs(row_a, row_b, best_coeff, lam, geom, args.n_modes, mode)
        cal = calibration_eval(rates, pulse, tau, inp, args.calibration_dt_ms)
        response = response_scan_nodes(rates, full_grid_nodes(), args.grid_dt_ms)
        min_r = float(np.nanmin(response["r"]))
        max_score = float(np.nanmax(cal["score"]))
        best_joint = max(max_score / inp.cutoff_mid, inp.threshold / max(min_r, 1e-300))

    return {
        "coeffs": best_coeff,
        "search_max_score": max_score,
        "search_min_r_dt0p1": min_r,
        "search_joint_barrier": best_joint,
        "history": history,
    }


def direct_coeffs(mode: str, geom: dict[str, np.ndarray], n_modes: int) -> np.ndarray:
    n_pcs = geom["pc_basis"].shape[0]
    dim = n_modes * n_pcs + (0 if mode == "native_aux" else n_modes * 2)
    return np.zeros(dim, dtype=float)


def checkpoint_key(a: int, b: int, mode: str, label: str) -> str:
    x, y = sorted((int(a), int(b)))
    return f"{label}_{x}_{y}_{mode}"


def save_verify_checkpoint(chkdir: Path, key: str, summary: dict[str, Any], frame: pd.DataFrame, coeffs: np.ndarray) -> None:
    d = chkdir / key
    d.mkdir(parents=True, exist_ok=True)
    json_dump(d / "summary.json", {**summary, "coefficients": [float(x) for x in coeffs]})
    frame.to_csv(d / "path_points.csv.gz", index=False, compression="gzip")


def load_verify_checkpoint(chkdir: Path, key: str) -> tuple[dict[str, Any], pd.DataFrame, np.ndarray] | None:
    d = chkdir / key
    if not (d / "summary.json").exists() or not (d / "path_points.csv.gz").exists():
        return None
    s = json.loads((d / "summary.json").read_text(encoding="utf-8"))
    c = np.asarray(s.pop("coefficients"), dtype=float)
    f = pd.read_csv(d / "path_points.csv.gz")
    return s, f, c


def certify_internal_native_forest(
    inp: Inputs,
    geom: dict[str, np.ndarray],
    bounds: dict[str, np.ndarray],
    comp_map: dict[int, int],
    args: argparse.Namespace,
    chkdir: Path,
) -> tuple[pd.DataFrame, list[pd.DataFrame], UnionFind]:
    ids = sorted(inp.confirmed["sample_id"].astype(int).tolist())
    uf = UnionFind(ids)
    rows = []
    frames: list[pd.DataFrame] = []
    lookup = inp.confirmed.set_index("sample_id", drop=False)

    refined_lookup = {}
    if len(inp.refined_edges26):
        for r in inp.refined_edges26.itertuples(index=False):
            refined_lookup[str(r.edge_id)] = r

    native_edges = inp.coarse_edges[inp.coarse_edges["native_joint_pass_mid_dt0p1"].astype(bool)].copy()
    native_edges = native_edges.sort_values(["distance", "edge_id"]).reset_index(drop=True)

    # Work component by component, selecting only edges required for a spanning forest.
    component_ids = sorted(set(comp_map.values()))
    for cid in component_ids:
        members = sorted([sid for sid, c in comp_map.items() if c == cid])
        if len(members) <= 1:
            continue
        member_set = set(members)
        local = native_edges[
            native_edges["sample_id_a"].astype(int).isin(member_set)
            & native_edges["sample_id_b"].astype(int).isin(member_set)
        ].copy()
        for er in local.itertuples(index=False):
            a, b = int(er.sample_id_a), int(er.sample_id_b)
            if uf.find(a) == uf.find(b):
                continue
            edge_id = str(er.edge_id)
            reused = refined_lookup.get(edge_id)
            if reused is not None and bool(getattr(reused, "native_joint_pass_mid_dt0p05")):
                uf.union(a, b)
                rows.append({
                    "edge_id": edge_id,
                    "sample_id_a": a,
                    "sample_id_b": b,
                    "component_id_step26_coarse": cid,
                    "source": "step26_refined_reuse",
                    "mode": "native_aux",
                    "midpoint_cutoff_pass": True,
                    "strict_cutoff_pass": bool(getattr(reused, "native_joint_pass_strict_dt0p05")),
                    "max_calibration_score": float(getattr(reused, "max_score_interp_aux")),
                    "min_r_dt0p05": float(getattr(reused, "min_r_support_dt0p05")),
                    "min_control_plateau_at_support_node_dt0p05": float(getattr(reused, "min_control_plateau_at_support_node_dt0p05")),
                })
                continue

            key = checkpoint_key(a, b, "native_aux", "internal")
            loaded = load_verify_checkpoint(chkdir, key)
            coeff = direct_coeffs("native_aux", geom, args.n_modes)
            if loaded is None:
                summary, frame = verify_path(
                    lookup.loc[a], lookup.loc[b], coeff, "native_aux", inp, geom, bounds,
                    args.n_modes, args.verify_points, args.grid_dt_ms, args.confirm_dt_ms, args.calibration_dt_ms,
                )
                save_verify_checkpoint(chkdir, key, summary, frame, coeff)
            else:
                summary, frame, coeff = loaded
            frames.append(frame.assign(path_role="internal_native"))
            rows.append({
                "edge_id": edge_id,
                "sample_id_a": a,
                "sample_id_b": b,
                "component_id_step26_coarse": cid,
                "source": "step27_verified",
                "mode": "native_aux",
                **{k: summary[k] for k in ["midpoint_cutoff_pass", "strict_cutoff_pass", "max_calibration_score", "min_r_dt0p05", "min_control_plateau_at_support_node_dt0p05"]},
            })
            if bool(summary["midpoint_cutoff_pass"]):
                uf.union(a, b)
            # Stop early if this Step26 component is now certified internally.
            roots = {uf.find(x) for x in members}
            if len(roots) == 1:
                break
        roots = {uf.find(x) for x in members}
        if len(roots) != 1:
            raise RuntimeError(
                f"Could not certify internal native connectivity of Step26 coarse component {cid}; "
                f"remaining certified subcomponents={len(roots)}"
            )
    return pd.DataFrame(rows), frames, uf


def try_bridge_pair(
    row_a: pd.Series,
    row_b: pd.Series,
    inp: Inputs,
    geom: dict[str, np.ndarray],
    bounds: dict[str, np.ndarray],
    target_nodes: list[tuple[float, float]],
    args: argparse.Namespace,
    chkdir: Path,
    seed_base: int,
) -> tuple[dict[str, Any], pd.DataFrame | None]:
    a, b = int(row_a["sample_id"]), int(row_b["sample_id"])
    attempts = []

    # 1) Direct endpoint-preserving log-linear path.
    mode = "native_aux"
    coeff = direct_coeffs(mode, geom, args.n_modes)
    key = checkpoint_key(a, b, mode, "bridge_direct")
    loaded = load_verify_checkpoint(chkdir, key)
    if loaded is None:
        summary, frame = verify_path(
            row_a, row_b, coeff, mode, inp, geom, bounds,
            args.n_modes, args.verify_points, args.grid_dt_ms, args.confirm_dt_ms, args.calibration_dt_ms,
        )
        save_verify_checkpoint(chkdir, key, summary, frame, coeff)
    else:
        summary, frame, coeff = loaded
    attempts.append({"stage": "direct", **summary})
    if bool(summary["midpoint_cutoff_pass"]):
        return {
            "success": True,
            "success_stage": "direct",
            "success_mode": mode,
            "summary": summary,
            "coeffs": coeff,
            "attempts": attempts,
        }, frame.assign(path_role="bridge_direct")

    # 2) Curved kinetic path, with native endpoint-preserving aux interpolation.
    mode = "native_aux"
    key = checkpoint_key(a, b, mode, "bridge_curved")
    loaded = load_verify_checkpoint(chkdir, key)
    if loaded is None:
        search = evolutionary_search(
            row_a, row_b, mode, inp, geom, bounds, target_nodes, args,
            seed=seed_base + a * 1009 + b * 9173,
        )
        coeff = search["coeffs"]
        summary, frame = verify_path(
            row_a, row_b, coeff, mode, inp, geom, bounds,
            args.n_modes, args.verify_points, args.grid_dt_ms, args.confirm_dt_ms, args.calibration_dt_ms,
        )
        summary.update({
            "search_max_score": float(search["search_max_score"]),
            "search_min_r_dt0p1": float(search["search_min_r_dt0p1"]),
            "search_joint_barrier": float(search["search_joint_barrier"]),
        })
        save_verify_checkpoint(chkdir, key, summary, frame, coeff)
        json_dump(chkdir / key / "search_history.json", search["history"])
    else:
        summary, frame, coeff = loaded
    attempts.append({"stage": "curved_native", **summary})
    if bool(summary["midpoint_cutoff_pass"]):
        return {
            "success": True,
            "success_stage": "curved_native",
            "success_mode": mode,
            "summary": summary,
            "coeffs": coeff,
            "attempts": attempts,
        }, frame.assign(path_role="bridge_curved_native")

    if args.disable_flex_aux:
        return {"success": False, "success_stage": None, "success_mode": None, "attempts": attempts}, None

    # 3) Curved full parameter path. Auxiliary coordinates also bend smoothly,
    # but their endpoint values remain exactly fixed.
    mode = "flex_aux"
    key = checkpoint_key(a, b, mode, "bridge_curved")
    loaded = load_verify_checkpoint(chkdir, key)
    if loaded is None:
        native_coeff = np.asarray(coeff, dtype=float)
        init = np.concatenate([native_coeff, np.zeros(args.n_modes * 2, dtype=float)])
        search = evolutionary_search(
            row_a, row_b, mode, inp, geom, bounds, target_nodes, args,
            seed=seed_base + a * 1217 + b * 7919 + 31,
            initial=init,
        )
        coeff_f = search["coeffs"]
        summary, frame = verify_path(
            row_a, row_b, coeff_f, mode, inp, geom, bounds,
            args.n_modes, args.verify_points, args.grid_dt_ms, args.confirm_dt_ms, args.calibration_dt_ms,
        )
        summary.update({
            "search_max_score": float(search["search_max_score"]),
            "search_min_r_dt0p1": float(search["search_min_r_dt0p1"]),
            "search_joint_barrier": float(search["search_joint_barrier"]),
        })
        save_verify_checkpoint(chkdir, key, summary, frame, coeff_f)
        json_dump(chkdir / key / "search_history.json", search["history"])
        coeff = coeff_f
    else:
        summary, frame, coeff = loaded
    attempts.append({"stage": "curved_flex_aux", **summary})
    if bool(summary["midpoint_cutoff_pass"]):
        return {
            "success": True,
            "success_stage": "curved_flex_aux",
            "success_mode": mode,
            "summary": summary,
            "coeffs": coeff,
            "attempts": attempts,
        }, frame.assign(path_role="bridge_curved_flex_aux")
    return {"success": False, "success_stage": None, "success_mode": None, "attempts": attempts}, None


def components_table(ids: list[int], edges: pd.DataFrame, pass_col: str, criterion: str) -> pd.DataFrame:
    uf = UnionFind(ids)
    if len(edges):
        for r in edges.itertuples(index=False):
            if bool(getattr(r, pass_col)):
                uf.union(int(r.sample_id_a), int(r.sample_id_b))
    comps = sorted(uf.components().values(), key=lambda x: (-len(x), min(x)))
    rows = []
    for cid, members in enumerate(comps):
        for sid in sorted(members):
            rows.append({"criterion": criterion, "component_id": cid, "component_size": len(members), "sample_id": sid})
    return pd.DataFrame(rows)


def make_figures(inp: Inputs, network_edges: pd.DataFrame, path_points: pd.DataFrame, attempts: pd.DataFrame, outdir: Path) -> None:
    strong_ids = set(inp.confirmed["sample_id"].astype(int))
    geo = inp.geometry[inp.geometry["sample_id"].astype(int).isin(strong_ids)].copy()
    pos = {int(r.sample_id): (float(r.PC1), float(r.PC2)) for r in geo.itertuples(index=False)}

    fig, ax = plt.subplots(figsize=(7.2, 6.0))
    ax.scatter(geo["PC1"], geo["PC2"], s=35, zorder=3)
    if len(network_edges):
        for r in network_edges.itertuples(index=False):
            a, b = int(r.sample_id_a), int(r.sample_id_b)
            if a not in pos or b not in pos:
                continue
            x = [pos[a][0], pos[b][0]]
            y = [pos[a][1], pos[b][1]]
            lw = 2.2 if str(r.edge_role).startswith("bridge") else 0.8
            alpha = 0.9 if str(r.edge_role).startswith("bridge") else 0.45
            ax.plot(x, y, linewidth=lw, alpha=alpha, zorder=1)
    for sid in [1900, 10052, 27084, 37318]:
        if sid in pos:
            ax.text(pos[sid][0], pos[sid][1], str(sid), fontsize=8)
    ax.set_xlabel("PC1 of standardized log-rate space")
    ax.set_ylabel("PC2 of standardized log-rate space")
    ax.set_title("Step 27 verified connectivity network")
    fig.tight_layout()
    fig.savefig(outdir / "Fig27A_verified_network.pdf")
    fig.savefig(outdir / "Fig27A_verified_network.png", dpi=180)
    plt.close(fig)

    if len(attempts):
        a = attempts.copy()
        a["normalized_score_barrier"] = a["max_calibration_score"] / inp.cutoff_mid
        a["normalized_strong_barrier"] = inp.threshold / a["min_r_dt0p05"].clip(lower=1e-300)
        fig, ax = plt.subplots(figsize=(7.2, 5.2))
        x = np.arange(len(a))
        ax.scatter(x, a["normalized_score_barrier"], label="calibration barrier")
        ax.scatter(x, a["normalized_strong_barrier"], label="strong barrier")
        ax.axhline(1.0, linestyle="--", linewidth=1.0)
        ax.set_xlabel("Bridge attempt")
        ax.set_ylabel("Normalized constraint barrier")
        ax.set_title("Straight and optimized bridge bottlenecks")
        ax.legend()
        fig.tight_layout()
        fig.savefig(outdir / "Fig27B_bridge_barriers.pdf")
        fig.savefig(outdir / "Fig27B_bridge_barriers.png", dpi=180)
        plt.close(fig)

    if len(path_points):
        bridges = path_points[path_points["path_role"].astype(str).str.startswith("bridge")].copy()
        if len(bridges):
            fig, ax = plt.subplots(figsize=(7.2, 5.2))
            for (a, b, role), g in bridges.groupby(["sample_id_a", "sample_id_b", "path_role"]):
                ax.plot(g["lambda"], g["calibration_score"] / inp.cutoff_mid, label=f"{a}-{b} score")
                ax.plot(g["lambda"], inp.threshold / g["r_plateau_dt0p05"].clip(lower=1e-300), linestyle="--", label=f"{a}-{b} strong")
            ax.axhline(1.0, linestyle=":", linewidth=1.0)
            ax.set_xlabel("Path coordinate lambda")
            ax.set_ylabel("Normalized constraint barrier")
            ax.set_title("Verified curved bridge profiles")
            if bridges[["sample_id_a", "sample_id_b"]].drop_duplicates().shape[0] <= 5:
                ax.legend(fontsize=7)
            fig.tight_layout()
            fig.savefig(outdir / "Fig27C_verified_bridge_profiles.pdf")
            fig.savefig(outdir / "Fig27C_verified_bridge_profiles.png", dpi=180)
            plt.close(fig)


def self_test() -> None:
    rng = np.random.default_rng(1)
    # Test endpoint preservation of the sine path parameterization.
    cols = eng.RATE_COLS + eng.AUX_COLS + ["sample_id"]
    a = pd.Series({c: 1.0 + 0.1 * i for i, c in enumerate(eng.RATE_COLS)})
    b = pd.Series({c: 1.5 + 0.05 * i for i, c in enumerate(eng.RATE_COLS)})
    a["pulse_amplitude_au"] = 0.5; b["pulse_amplitude_au"] = 1.2
    a["glutamate_tau_ms"] = 5.0; b["glutamate_tau_ms"] = 20.0
    a["sample_id"] = 1; b["sample_id"] = 2
    geom = {
        "rate_log_mean": np.zeros(12),
        "rate_log_sd": np.ones(12),
        "pc_basis": np.eye(12)[:3],
        "aux_log_mean": np.zeros(2),
        "aux_log_sd": np.ones(2),
    }
    lam = np.linspace(0, 1, 9)
    coeff = rng.normal(0, 0.2, size=2 * 3 + 2 * 2)
    rates, pulse, tau, _ = path_from_coeffs(a, b, coeff, lam, geom, 2, "flex_aux")
    assert np.allclose(rates[0], a[eng.RATE_COLS].to_numpy(float))
    assert np.allclose(rates[-1], b[eng.RATE_COLS].to_numpy(float))
    assert abs(pulse[0] - float(a["pulse_amplitude_au"])) < 1e-12
    assert abs(pulse[-1] - float(b["pulse_amplitude_au"])) < 1e-12
    assert abs(tau[0] - float(a["glutamate_tau_ms"])) < 1e-12
    assert abs(tau[-1] - float(b["glutamate_tau_ms"])) < 1e-12
    uf = UnionFind([1, 2, 3])
    assert uf.union(1, 2)
    assert len(uf.components()) == 2
    print("SELF_TEST_PASS")


def main() -> None:
    args = parse_args()
    if args.self_test:
        self_test()
        return
    if args.search_points < 5 or args.verify_points < args.search_points:
        raise ValueError("Require verify_points >= search_points >= 5")
    if args.n_pcs < 1 or args.n_modes < 1:
        raise ValueError("n_pcs and n_modes must be >=1")
    if args.population < 8:
        raise ValueError("population must be >=8")

    set_num_threads(args.threads)
    step25 = discover_step25(args.step25)
    step26 = discover_step26(args.step26)
    inp = load_inputs(step25, step26)
    replay = replay_gate(inp, args.calibration_dt_ms, args.confirm_dt_ms)
    if not replay["pass"]:
        raise RuntimeError(f"Engine replay gate failed: {replay}")
    bounds, bounds_source = load_full_broad_bounds(inp, args.allow_accepted_bounds_fallback)
    geom = geometry_basis(inp, args.n_pcs)
    args.n_pcs = int(geom["pc_basis"].shape[0])
    target_nodes = target_forcing_nodes(inp)
    comp_map = native_component_map(inp)
    component_priority = component_priority_table(inp, comp_map)

    preflight = {
        "step25": str(step25),
        "step26": str(step26),
        "step25_status": inp.decision25.get("status"),
        "step26_status": inp.decision26.get("status"),
        "n_confirmed_strong": int(len(inp.confirmed)),
        "n_step26_native_coarse_components": int(len(set(comp_map.values()))),
        "strong_threshold": inp.threshold,
        "midpoint_cutoff": inp.cutoff_mid,
        "strict_cutoff": inp.cutoff_strict,
        "bounds_source": bounds_source,
        "target_forcing_nodes": [[g, t] for g, t in target_nodes],
        "engine_replay": replay,
    }
    if args.preflight_only:
        print(json.dumps(preflight, indent=2))
        return

    outdir = args.output.expanduser().resolve()
    if outdir.exists():
        if args.force:
            shutil.rmtree(outdir)
        elif not args.resume:
            raise FileExistsError(f"Output exists: {outdir}. Set STEP27_RESUME=1 to continue or STEP27_FORCE=1 to replace it.")
    outdir.mkdir(parents=True, exist_ok=True)
    chkdir = outdir / "_checkpoints"
    chkdir.mkdir(parents=True, exist_ok=True)
    json_dump(outdir / "00_input_audit.json", preflight)

    # Save Step26 native component membership and bridge priority before new optimization.
    comp_rows = []
    for sid in sorted(comp_map):
        comp_rows.append({"sample_id": sid, "step26_native_coarse_component": comp_map[sid]})
    pd.DataFrame(comp_rows).to_csv(outdir / "01_step26_native_components.csv", index=False)
    component_priority.to_csv(outdir / "02_component_pair_priority.csv", index=False)

    # Exact dt=0.05 certification of a native spanning forest within each Step26 coarse component.
    internal_edges, internal_frames, uf = certify_internal_native_forest(inp, geom, bounds, comp_map, args, chkdir)
    internal_edges.to_csv(outdir / "03_internal_native_certification.csv", index=False)

    lookup = inp.confirmed.set_index("sample_id", drop=False)
    bridge_attempt_rows = []
    bridge_edge_rows = []
    path_frames = list(internal_frames)
    tested_component_pairs = 0

    # Kruskal-like component bridge search. Each original Step26 component pair is
    # visited by increasing endpoint distance; already merged groups are skipped.
    for cp in component_priority.itertuples(index=False):
        if len(uf.components()) == 1:
            break
        if tested_component_pairs >= args.max_component_pairs:
            break
        ca, cb = int(cp.component_a), int(cp.component_b)
        members_a = [sid for sid, c in comp_map.items() if c == ca]
        members_b = [sid for sid, c in comp_map.items() if c == cb]
        if all(uf.find(a) == uf.find(b) for a in members_a for b in members_b):
            continue
        tested_component_pairs += 1
        pairs = endpoint_pairs_for_components(inp, comp_map, ca, cb, args.pair_attempts)
        pair_success = False
        for pr in pairs.itertuples(index=False):
            a, b = int(pr.sample_id_a), int(pr.sample_id_b)
            if uf.find(a) == uf.find(b):
                continue
            result, frame = try_bridge_pair(
                lookup.loc[a], lookup.loc[b], inp, geom, bounds, target_nodes, args, chkdir,
                seed_base=args.seed + tested_component_pairs * 100003,
            )
            for at in result["attempts"]:
                bridge_attempt_rows.append({
                    "component_a": ca,
                    "component_b": cb,
                    "endpoint_distance": float(pr.distance),
                    "sample_id_a": a,
                    "sample_id_b": b,
                    **at,
                })
            if result["success"]:
                s = result["summary"]
                uf.union(a, b)
                bridge_edge_rows.append({
                    "sample_id_a": a,
                    "sample_id_b": b,
                    "component_a_step26": ca,
                    "component_b_step26": cb,
                    "endpoint_distance": float(pr.distance),
                    "success_stage": result["success_stage"],
                    "success_mode": result["success_mode"],
                    "midpoint_cutoff_pass": bool(s["midpoint_cutoff_pass"]),
                    "strict_cutoff_pass": bool(s["strict_cutoff_pass"]),
                    "max_calibration_score": float(s["max_calibration_score"]),
                    "min_r_dt0p05": float(s["min_r_dt0p05"]),
                    "min_control_plateau_at_support_node_dt0p05": float(s["min_control_plateau_at_support_node_dt0p05"]),
                })
                if frame is not None:
                    path_frames.append(frame)
                pair_success = True
                print(f"Verified bridge {a}-{b}: {result['success_stage']}", flush=True)
                break
            print(f"Bridge attempt failed {a}-{b}", flush=True)
        if not pair_success:
            print(f"No verified bridge for Step26 components {ca}-{cb} within endpoint budget", flush=True)

    bridge_attempts = pd.DataFrame(bridge_attempt_rows)
    bridge_edges = pd.DataFrame(bridge_edge_rows)
    bridge_attempts.to_csv(outdir / "04_bridge_attempts.csv", index=False)
    bridge_edges.to_csv(outdir / "05_verified_bridges.csv", index=False)

    if path_frames:
        path_points = pd.concat(path_frames, ignore_index=True)
    else:
        path_points = pd.DataFrame()
    path_points.to_csv(outdir / "06_verified_path_points.csv.gz", index=False, compression="gzip")

    # Build one verified network table. Internal edges are native. Bridges may be native or flex_aux.
    net_rows = []
    for r in internal_edges.itertuples(index=False):
        net_rows.append({
            "sample_id_a": int(r.sample_id_a),
            "sample_id_b": int(r.sample_id_b),
            "edge_role": "internal_native",
            "mode": "native_aux",
            "midpoint_cutoff_pass": bool(r.midpoint_cutoff_pass),
            "strict_cutoff_pass": bool(r.strict_cutoff_pass),
            "max_calibration_score": float(r.max_calibration_score),
            "min_r_dt0p05": float(r.min_r_dt0p05),
            "min_control_plateau_at_support_node_dt0p05": float(r.min_control_plateau_at_support_node_dt0p05),
        })
    if len(bridge_edges):
        for r in bridge_edges.itertuples(index=False):
            net_rows.append({
                "sample_id_a": int(r.sample_id_a),
                "sample_id_b": int(r.sample_id_b),
                "edge_role": "bridge_" + str(r.success_stage),
                "mode": str(r.success_mode),
                "midpoint_cutoff_pass": bool(r.midpoint_cutoff_pass),
                "strict_cutoff_pass": bool(r.strict_cutoff_pass),
                "max_calibration_score": float(r.max_calibration_score),
                "min_r_dt0p05": float(r.min_r_dt0p05),
                "min_control_plateau_at_support_node_dt0p05": float(r.min_control_plateau_at_support_node_dt0p05),
            })
    network_edges = pd.DataFrame(net_rows)
    network_edges.to_csv(outdir / "07_verified_network_edges.csv", index=False)

    ids = sorted(inp.confirmed["sample_id"].astype(int).tolist())
    full_comp = components_table(ids, network_edges, "midpoint_cutoff_pass", "verified_midpoint")
    strict_comp = components_table(ids, network_edges, "strict_cutoff_pass", "verified_strict")
    comps = pd.concat([full_comp, strict_comp], ignore_index=True)
    comps.to_csv(outdir / "08_final_components.csv", index=False)

    ncomp_mid = int(full_comp["component_id"].nunique())
    ncomp_strict = int(strict_comp["component_id"].nunique())
    all_native = bool(len(bridge_edges) == 0 or bridge_edges["success_mode"].astype(str).eq("native_aux").all())
    if ncomp_mid == 1 and all_native:
        status = "CURVED_KINETIC_CONNECTIVITY_ESTABLISHED_NATIVE_AUX"
    elif ncomp_mid == 1:
        status = "CURVED_FULL_PARAMETER_CONNECTIVITY_ESTABLISHED"
    elif len(bridge_edges):
        status = "PARTIAL_CURVED_CONNECTIVITY_ONLY"
    else:
        status = "CURVED_CONNECTIVITY_NOT_ESTABLISHED_WITHIN_BUDGET"

    bottleneck = network_edges.copy()
    if len(bottleneck):
        bottleneck["calibration_margin"] = inp.cutoff_mid - bottleneck["max_calibration_score"]
        bottleneck["strong_margin"] = bottleneck["min_r_dt0p05"] - inp.threshold
        bottleneck = bottleneck.sort_values(["calibration_margin", "strong_margin"])
    bottleneck.to_csv(outdir / "09_network_bottlenecks.csv", index=False)

    decision = {
        "status": status,
        "interpretation_scope": (
            "Finite-resolution connectivity of the 40 corrected strong vectors under endpoint-preserving curved paths. "
            "Native paths bend only the 12 kinetic rates while auxiliary calibration coordinates follow geometric "
            "endpoint interpolation. Secondary flex_aux paths also bend the two auxiliary coordinates but keep their "
            "endpoint values fixed. Failure within this search budget is not a proof of global disconnection."
        ),
        "n_strong_candidates": 40,
        "step26_native_coarse_components": int(len(set(comp_map.values()))),
        "tested_component_pairs": int(tested_component_pairs),
        "verified_bridges": int(len(bridge_edges)),
        "verified_native_bridges": int((bridge_edges["success_mode"] == "native_aux").sum()) if len(bridge_edges) else 0,
        "verified_flex_aux_bridges": int((bridge_edges["success_mode"] == "flex_aux").sum()) if len(bridge_edges) else 0,
        "n_components_midpoint": ncomp_mid,
        "n_components_strict": ncomp_strict,
        "midpoint_cutoff": inp.cutoff_mid,
        "strict_cutoff": inp.cutoff_strict,
        "strong_threshold": inp.threshold,
        "no_denominator_exclusion_applied": True,
        "engine_replay_gate_pass": bool(replay["pass"]),
        "full_broad_bounds_source": bounds_source,
    }
    json_dump(outdir / "10_scientific_decision.json", decision)

    make_figures(inp, network_edges, path_points, bridge_attempts, outdir)

    readme_lines = [
        "# Step 27 results - curved bridge optimization",
        "",
        f"Status: **{status}**",
        "",
        f"Step26 native coarse components entering Step27: **{len(set(comp_map.values()))}**.",
        f"Verified bridges found: **{len(bridge_edges)}**.",
        f"Final midpoint-cutoff components: **{ncomp_mid}**.",
        f"Final strict-cutoff components: **{ncomp_strict}**.",
        "",
        "Primary paths preserve endpoint auxiliary values and bend only kinetic coordinates. Secondary flex_aux paths",
        "also bend auxiliary coordinates smoothly while preserving both endpoints. No fixed-A/fixed-B whole-edge schedule",
        "is used as a topological connector because such a schedule does not terminate at both stored endpoint auxiliary values.",
        "",
        "Every accepted Step27 path was verified on the full frozen 80-node response grid with dt=0.05-ms confirmation.",
        "No denominator threshold was introduced.",
        "",
        "Failure to establish one connected network within the finite bridge-search budget is not proof of global disconnection.",
    ]
    (outdir / "README_RESULTS.md").write_text("\n".join(readme_lines) + "\n", encoding="utf-8")

    run_summary = {
        "pipeline_step": "27_curved_bridge_optimization",
        "script_version": SCRIPT_VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "numba_threads": get_num_threads(),
        "step25": str(step25),
        "step26": str(step26),
        "search_points": args.search_points,
        "verify_points": args.verify_points,
        "n_pcs": args.n_pcs,
        "n_modes": args.n_modes,
        "population": args.population,
        "generations": args.generations,
        "joint_generations": args.joint_generations,
        "max_component_pairs": args.max_component_pairs,
        "pair_attempts": args.pair_attempts,
        "decision": decision,
    }
    json_dump(outdir / "run_summary.json", run_summary)
    print(json.dumps(decision, indent=2), flush=True)


if __name__ == "__main__":
    main()
