#!/usr/bin/env python3
"""NMDA Braess Step 26 v1.0 — strong-compatible path connectivity.

Step 25 established 40 numerically confirmed strong, corrected-calibration-compatible
broad-prior vectors, but point-cloud clustering was unstable. Step 26 therefore tests
an explicit question that clustering cannot answer: can the sampled strong vectors be
linked by continuous log-linear paths that remain both strong and corrected-control
compatible?

The analysis is deliberately finite and conservative. A successful spanning network
supports connectivity of the sampled region along the tested paths; failure to find one
does not prove global disconnection.

Primary path definition
-----------------------
Kinetic rates (12 coordinates) are interpolated geometrically (linearly in log-rate
space). Corrected control compatibility is evaluated under three prespecified,
continuous auxiliary schedules for the historical calibration drive:
  1. interp_aux: geometric interpolation of pulse amplitude and decay time;
  2. fixed_A_aux: endpoint-A auxiliary values held fixed along the path;
  3. fixed_B_aux: endpoint-B auxiliary values held fixed along the path.

"native" connectivity uses only interp_aux. "any_schedule" connectivity allows an
edge only when one *single* prespecified schedule passes at every point on that edge;
pointwise switching between schedules is not allowed.

Strong-response support is always tested with the frozen final Base forcing grid. No
new denominator cutoff is introduced.
"""

from __future__ import annotations

import argparse
import hashlib
import heapq
import json
import math
import os
import platform
import shutil
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from numba import get_num_threads, set_num_threads

import step25_engine_snapshot as eng


SCRIPT_VERSION = "1.0.0"
EXPECTED_STEP25_STATUS = "KINETIC_REGION_STRUCTURE_UNRESOLVED"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Step 26 strong-compatible path connectivity")
    p.add_argument("--step25", default="auto", help="Step-25 results directory or 'auto'")
    p.add_argument(
        "--output",
        type=Path,
        default=Path("/root/nmda2/step_26/results_step_26_strong_connectivity"),
    )
    p.add_argument("--threads", type=int, default=2)
    p.add_argument("--coarse-points", type=int, default=21)
    p.add_argument("--refine-points", type=int, default=41)
    p.add_argument("--max-edges", type=int, default=180)
    p.add_argument("--knn-k", type=int, default=4)
    p.add_argument("--grid-dt-ms", type=float, default=0.1)
    p.add_argument("--confirm-dt-ms", type=float, default=0.05)
    p.add_argument("--calibration-dt-ms", type=float, default=0.05)
    p.add_argument("--allow-step25-status", action="store_true")
    p.add_argument("--force", action="store_true")
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
    required = {"05_strong_candidates_confirmed.csv", "15_scientific_decision.json", "00_input_audit.json"}
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


@dataclass
class Inputs:
    step25: Path
    audit: dict[str, Any]
    decision: dict[str, Any]
    targets: dict[str, Any]
    confirmed: pd.DataFrame
    accepted: pd.DataFrame
    reps: pd.DataFrame
    pca: pd.DataFrame
    pairwise: pd.DataFrame
    mst: pd.DataFrame
    threshold: float
    cutoff_mid: float
    cutoff_strict: float


def load_inputs(step25: Path, allow_status: bool) -> Inputs:
    names = [
        "00_input_audit.json",
        "01_corrected_target_reconstruction.json",
        "02_corrected_accepted_cohorts.csv.gz",
        "05_strong_candidates_confirmed.csv",
        "08_geometry_coordinates.csv.gz",
        "09_geometry_pairwise_distances.csv",
        "11C_strong_mst_edges.csv",
        "12_representative_centres.csv",
        "15_scientific_decision.json",
        "run_summary.json",
    ]
    for name in names:
        if not (step25 / name).exists():
            raise FileNotFoundError(step25 / name)

    audit = json.loads((step25 / "00_input_audit.json").read_text(encoding="utf-8"))
    decision = json.loads((step25 / "15_scientific_decision.json").read_text(encoding="utf-8"))
    targets = json.loads((step25 / "01_corrected_target_reconstruction.json").read_text(encoding="utf-8"))
    confirmed = pd.read_csv(step25 / "05_strong_candidates_confirmed.csv")
    accepted = pd.read_csv(step25 / "02_corrected_accepted_cohorts.csv.gz")
    reps = pd.read_csv(step25 / "12_representative_centres.csv")
    pca = pd.read_csv(step25 / "08_geometry_coordinates.csv.gz")
    pairwise = pd.read_csv(step25 / "09_geometry_pairwise_distances.csv")
    mst = pd.read_csv(step25 / "11C_strong_mst_edges.csv")

    status = str(decision.get("status"))
    if status != EXPECTED_STEP25_STATUS and not allow_status:
        raise RuntimeError(
            f"Step25 status is {status!r}; expected {EXPECTED_STEP25_STATUS!r}. "
            "Use --allow-step25-status only for diagnostic reuse."
        )
    if not bool(decision.get("numerical_gate_pass", False)):
        raise RuntimeError("Step25 numerical replay gate did not pass")
    if not bool(decision.get("no_denominator_exclusion_applied", False)):
        raise RuntimeError("Step25 unexpectedly reports a denominator exclusion")
    if len(confirmed) != int(decision.get("confirmed_strong_n_dt0p05", -1)):
        raise RuntimeError("Confirmed strong count does not match Step25 decision")
    if not confirmed["strong_confirmed_dt0p05"].astype(bool).all():
        raise RuntimeError("Not all Step25 confirmed candidates are strong at dt=0.05 ms")
    if confirmed["sample_id"].duplicated().any():
        raise RuntimeError("Duplicate sample_id in Step25 confirmed table")

    accepted_broad = accepted[accepted["prior"].astype(str).str.lower().eq("broad")].copy()
    if len(accepted_broad) != 5000:
        raise RuntimeError(f"Expected 5000 corrected accepted broad rows; got {len(accepted_broad)}")

    threshold = float(audit["strong_threshold"])
    cutoff_mid = float(audit["corrected_broad_cutoff_midpoint"])
    cutoff_strict = float(audit["corrected_broad_max_accepted_score"])

    return Inputs(
        step25=step25,
        audit=audit,
        decision=decision,
        targets=targets,
        confirmed=confirmed,
        accepted=accepted_broad,
        reps=reps,
        pca=pca,
        pairwise=pairwise,
        mst=mst,
        threshold=threshold,
        cutoff_mid=cutoff_mid,
        cutoff_strict=cutoff_strict,
    )


def replay_gate(inp: Inputs, calib_dt: float, confirm_dt: float) -> dict[str, Any]:
    ids = list(dict.fromkeys([37318, *inp.reps["sample_id"].astype(int).tolist()]))
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
        out, _ = eng.run_base_condition(
            rates[i : i + 1], float(row["G_peak"]), float(row["tau_G_ms"]), confirm_dt
        )
        r_err.append(abs(float(out[0, 2]) - float(row["r_plateau_dt0p05"])))
    r_err = np.asarray(r_err)

    gate = bool(np.max(score_err) < 1e-8 and np.max(r_err) < 5e-4)
    return {
        "sample_ids": ids,
        "max_abs_calibration_score_error": float(np.max(score_err)),
        "max_abs_response_ratio_error": float(np.max(r_err)),
        "pass": gate,
    }


def build_edge_pool(inp: Inputs, knn_k: int) -> pd.DataFrame:
    ids = sorted(inp.confirmed["sample_id"].astype(int).tolist())
    idset = set(ids)
    pw = inp.pairwise.copy()
    pw = pw[pw["sample_id_a"].astype(int).isin(idset) & pw["sample_id_b"].astype(int).isin(idset)].copy()
    pw["sample_id_a"] = pw["sample_id_a"].astype(int)
    pw["sample_id_b"] = pw["sample_id_b"].astype(int)
    pw["edge_key"] = pw.apply(lambda r: tuple(sorted((int(r["sample_id_a"]), int(r["sample_id_b"])))), axis=1)
    if len(pw) != len(ids) * (len(ids) - 1) // 2:
        raise RuntimeError("Step25 pairwise table is not complete for the 40 strong candidates")

    mst_keys = {
        tuple(sorted((int(a), int(b))))
        for a, b in inp.mst[["sample_id_a", "sample_id_b"]].itertuples(index=False, name=None)
    }

    dist = {row.edge_key: float(row.standardized_log_rate_distance) for row in pw.itertuples(index=False)}
    knn_keys: set[tuple[int, int]] = set()
    for a in ids:
        nbrs = []
        for b in ids:
            if a == b:
                continue
            key = tuple(sorted((a, b)))
            nbrs.append((dist[key], b, key))
        nbrs.sort()
        for _, _, key in nbrs[:knn_k]:
            knn_keys.add(key)

    rep_ids = [int(x) for x in inp.reps["sample_id"].tolist() if int(x) in idset]
    rep_keys = {tuple(sorted((a, b))) for i, a in enumerate(rep_ids) for b in rep_ids[i + 1 :]}

    rows = []
    for row in pw.itertuples(index=False):
        key = row.edge_key
        priority = 3
        tags = []
        if key in mst_keys:
            priority = min(priority, 0)
            tags.append("step25_mst")
        if key in knn_keys:
            priority = min(priority, 1)
            tags.append(f"knn{k_nn_label(knn_k)}")
        if key in rep_keys:
            priority = min(priority, 2)
            tags.append("representative_pair")
        rows.append(
            {
                "edge_id": f"E{key[0]}_{key[1]}",
                "sample_id_a": key[0],
                "sample_id_b": key[1],
                "distance": float(row.standardized_log_rate_distance),
                "priority": priority,
                "tags": ";".join(tags) if tags else "distance_pool",
            }
        )
    out = pd.DataFrame(rows).sort_values(["priority", "distance", "sample_id_a", "sample_id_b"]).reset_index(drop=True)
    out["pool_rank"] = np.arange(1, len(out) + 1)
    return out


def k_nn_label(k: int) -> str:
    return f"_{k}"


def log_interp(a: np.ndarray, b: np.ndarray, lam: np.ndarray) -> np.ndarray:
    if np.any(a <= 0) or np.any(b <= 0):
        raise ValueError("Log interpolation requires positive endpoints")
    return np.exp((1.0 - lam[:, None]) * np.log(a)[None, :] + lam[:, None] * np.log(b)[None, :])


def scalar_log_interp(a: float, b: float, lam: np.ndarray) -> np.ndarray:
    if a <= 0 or b <= 0:
        raise ValueError("Log interpolation requires positive endpoints")
    return np.exp((1.0 - lam) * math.log(a) + lam * math.log(b))


def calibration_profile(
    rates: np.ndarray,
    row_a: pd.Series,
    row_b: pd.Series,
    lam: np.ndarray,
    targets: dict[str, Any],
    dt_ms: float,
) -> dict[str, np.ndarray]:
    pa = float(row_a["pulse_amplitude_au"])
    pb = float(row_b["pulse_amplitude_au"])
    ta = float(row_a["glutamate_tau_ms"])
    tb = float(row_b["glutamate_tau_ms"])

    schedules = {
        "interp_aux": (scalar_log_interp(pa, pb, lam), scalar_log_interp(ta, tb, lam)),
        "fixed_A_aux": (np.full(len(lam), pa), np.full(len(lam), ta)),
        "fixed_B_aux": (np.full(len(lam), pb), np.full(len(lam), tb)),
    }
    out: dict[str, np.ndarray] = {}
    for name, (pulse, tau) in schedules.items():
        met = eng.calibration_ensemble(rates, pulse.astype(float), tau.astype(float), dt_ms)
        score = eng.corrected_score_from_metrics(met[:, 1], met[:, 2], met[:, 3], met[:, 5], targets)
        valid = met[:, 8] > 0.5
        score[~valid] = np.inf
        out[f"score_{name}"] = score
        out[f"H_{name}"] = met[:, 1]
        out[f"E_{name}"] = met[:, 2]
        out[f"L_{name}"] = met[:, 3]
        out[f"Jcorr_{name}"] = met[:, 5]
        out[f"valid_{name}"] = valid
    return out


def forcing_cache(dt_ms: float) -> dict[tuple[float, float], tuple[np.ndarray, np.ndarray, np.ndarray, float]]:
    return {
        (float(g), float(t)): eng.forcing_arrays(float(g), float(t), dt_ms)
        for g in eng.G_PEAK_GRID
        for t in eng.TAU_GRID
    }


def response_scan(
    rates: np.ndarray,
    dt_ms: float,
    target_nodes: list[tuple[float, float]],
    full_grid: bool,
    cache: dict[tuple[float, float], tuple[np.ndarray, np.ndarray, np.ndarray, float]],
    threshold: float,
) -> dict[str, np.ndarray]:
    n = len(rates)
    best_r = np.full(n, -np.inf, dtype=float)
    best_c = np.full(n, np.nan, dtype=float)
    best_b = np.full(n, np.nan, dtype=float)
    best_g = np.full(n, np.nan, dtype=float)
    best_tau = np.full(n, np.nan, dtype=float)
    min_state = np.full(n, np.nan, dtype=float)
    evaluated = np.zeros(n, dtype=int)
    full_eval = np.zeros(n, dtype=bool)

    all_nodes = [(float(g), float(t)) for g in eng.G_PEAK_GRID for t in eng.TAU_GRID]
    target_unique = list(dict.fromkeys(target_nodes))
    target_set = set(target_unique)

    def eval_nodes(indices: np.ndarray, nodes: Iterable[tuple[float, float]]) -> None:
        if len(indices) == 0:
            return
        rr = rates[indices]
        for g, tau in nodes:
            g0, gm, g1, _ = cache[(g, tau)]
            out = eng.base_ratio_ensemble(rr, g0, gm, g1, dt_ms)
            r = out[:, 2]
            better = np.isfinite(r) & (r > best_r[indices])
            if np.any(better):
                ii = indices[better]
                oo = out[better]
                best_r[ii] = oo[:, 2]
                best_c[ii] = oo[:, 0]
                best_b[ii] = oo[:, 1]
                best_g[ii] = g
                best_tau[ii] = tau
                min_state[ii] = np.minimum(oo[:, 3], oo[:, 4])
            evaluated[indices] += 1

    if full_grid:
        idx = np.arange(n, dtype=int)
        eval_nodes(idx, all_nodes)
        full_eval[:] = True
    else:
        idx = np.arange(n, dtype=int)
        eval_nodes(idx, target_unique)
        unresolved = np.where(~(best_r >= threshold))[0]
        if len(unresolved):
            remaining = [x for x in all_nodes if x not in target_set]
            eval_nodes(unresolved, remaining)
            full_eval[unresolved] = True

    return {
        "r": best_r,
        "control": best_c,
        "blocked": best_b,
        "G": best_g,
        "tau": best_tau,
        "minimum_state": min_state,
        "n_nodes_evaluated": evaluated,
        "full_grid_evaluated": full_eval,
    }


def test_edge(
    edge: pd.Series,
    cand: pd.DataFrame,
    n_points: int,
    stage: str,
    inp: Inputs,
    grid_dt: float,
    calib_dt: float,
    target_nodes: list[tuple[float, float]],
    cache_grid: dict[tuple[float, float], tuple[np.ndarray, np.ndarray, np.ndarray, float]],
    exact_full_grid: bool,
) -> tuple[dict[str, Any], pd.DataFrame]:
    a = int(edge["sample_id_a"])
    b = int(edge["sample_id_b"])
    row_a = cand.loc[a]
    row_b = cand.loc[b]
    lam = np.linspace(0.0, 1.0, n_points)
    rates = log_interp(
        row_a[eng.RATE_COLS].to_numpy(float), row_b[eng.RATE_COLS].to_numpy(float), lam
    )
    cal = calibration_profile(rates, row_a, row_b, lam, inp.targets, calib_dt)
    rsp = response_scan(rates, grid_dt, target_nodes, exact_full_grid, cache_grid, inp.threshold)

    strong = np.isfinite(rsp["r"]) & (rsp["r"] >= inp.threshold)
    comp_mid = {}
    comp_strict = {}
    for s in ["interp_aux", "fixed_A_aux", "fixed_B_aux"]:
        comp_mid[s] = np.isfinite(cal[f"score_{s}"]) & (cal[f"score_{s}"] <= inp.cutoff_mid)
        comp_strict[s] = np.isfinite(cal[f"score_{s}"]) & (cal[f"score_{s}"] <= inp.cutoff_strict)

    schedule_pass_mid = {s: bool(np.all(strong & comp_mid[s])) for s in comp_mid}
    schedule_pass_strict = {s: bool(np.all(strong & comp_strict[s])) for s in comp_strict}
    native_mid = schedule_pass_mid["interp_aux"]
    native_strict = schedule_pass_strict["interp_aux"]
    any_mid = any(schedule_pass_mid.values())
    any_strict = any(schedule_pass_strict.values())
    passing_schedule_mid = next((s for s in ["interp_aux", "fixed_A_aux", "fixed_B_aux"] if schedule_pass_mid[s]), "")
    passing_schedule_strict = next((s for s in ["interp_aux", "fixed_A_aux", "fixed_B_aux"] if schedule_pass_strict[s]), "")

    pointwise_best_score = np.minimum.reduce(
        [cal["score_interp_aux"], cal["score_fixed_A_aux"], cal["score_fixed_B_aux"]]
    )
    path = pd.DataFrame(
        {
            "edge_id": str(edge["edge_id"]),
            "sample_id_a": a,
            "sample_id_b": b,
            "distance": float(edge["distance"]),
            "stage": stage,
            "lambda": lam,
            "r_support_dt0p1": rsp["r"],
            "control_plateau_at_support_node_dt0p1": rsp["control"],
            "blocked_plateau_at_support_node_dt0p1": rsp["blocked"],
            "G_peak_support_node": rsp["G"],
            "tau_G_ms_support_node": rsp["tau"],
            "n_forcing_nodes_evaluated": rsp["n_nodes_evaluated"],
            "full_grid_evaluated": rsp["full_grid_evaluated"],
            "strong_dt0p1": strong,
            "score_interp_aux": cal["score_interp_aux"],
            "score_fixed_A_aux": cal["score_fixed_A_aux"],
            "score_fixed_B_aux": cal["score_fixed_B_aux"],
            "pointwise_best_score": pointwise_best_score,
            "compatible_interp_aux_mid": comp_mid["interp_aux"],
            "compatible_fixed_A_aux_mid": comp_mid["fixed_A_aux"],
            "compatible_fixed_B_aux_mid": comp_mid["fixed_B_aux"],
        }
    )

    summary = {
        "edge_id": str(edge["edge_id"]),
        "sample_id_a": a,
        "sample_id_b": b,
        "distance": float(edge["distance"]),
        "priority": int(edge["priority"]),
        "tags": str(edge["tags"]),
        "stage": stage,
        "n_path_points": int(n_points),
        "strong_path_pass_dt0p1": bool(np.all(strong)),
        "native_joint_pass_mid_dt0p1": native_mid,
        "native_joint_pass_strict_dt0p1": native_strict,
        "any_schedule_joint_pass_mid_dt0p1": any_mid,
        "any_schedule_joint_pass_strict_dt0p1": any_strict,
        "passing_schedule_mid": passing_schedule_mid,
        "passing_schedule_strict": passing_schedule_strict,
        "min_r_support_dt0p1": float(np.nanmin(rsp["r"])),
        "lambda_at_min_r": float(lam[int(np.nanargmin(rsp["r"]))]),
        "max_score_interp_aux": float(np.nanmax(cal["score_interp_aux"])),
        "max_score_fixed_A_aux": float(np.nanmax(cal["score_fixed_A_aux"])),
        "max_score_fixed_B_aux": float(np.nanmax(cal["score_fixed_B_aux"])),
        "min_control_plateau_at_support_node": float(np.nanmin(rsp["control"])),
        "median_control_plateau_at_support_node": float(np.nanmedian(rsp["control"])),
        "all_points_full_grid": bool(np.all(rsp["full_grid_evaluated"])),
    }
    return summary, path


class UnionFind:
    def __init__(self, items: Iterable[int]):
        self.parent = {int(x): int(x) for x in items}
        self.rank = {int(x): 0 for x in items}

    def find(self, x: int) -> int:
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

    def components(self) -> list[list[int]]:
        d: dict[int, list[int]] = {}
        for x in self.parent:
            d.setdefault(self.find(x), []).append(x)
        return [sorted(v) for v in d.values()]


def component_table(ids: list[int], edges: pd.DataFrame, pass_col: str, label: str) -> pd.DataFrame:
    uf = UnionFind(ids)
    if len(edges):
        for row in edges.itertuples(index=False):
            if bool(getattr(row, pass_col)):
                uf.union(int(row.sample_id_a), int(row.sample_id_b))
    comps = sorted(uf.components(), key=lambda x: (-len(x), x[0]))
    rows = []
    for ci, comp in enumerate(comps):
        for sid in comp:
            rows.append({"criterion": label, "component_id": ci, "component_size": len(comp), "sample_id": sid})
    return pd.DataFrame(rows)


def is_connected(ids: list[int], edges: pd.DataFrame, pass_col: str) -> bool:
    t = component_table(ids, edges, pass_col, pass_col)
    return t["component_id"].nunique() == 1


def kruskal_tree(ids: list[int], edges: pd.DataFrame, pass_col: str) -> list[str]:
    uf = UnionFind(ids)
    chosen: list[str] = []
    sub = edges[edges[pass_col].astype(bool)].sort_values(["distance", "edge_id"])
    for row in sub.itertuples(index=False):
        if uf.union(int(row.sample_id_a), int(row.sample_id_b)):
            chosen.append(str(row.edge_id))
            if len(chosen) == len(ids) - 1:
                break
    if len(uf.components()) != 1:
        return []
    return chosen


def exact_refine_with_dt_confirmation(
    edge: pd.Series,
    cand: pd.DataFrame,
    n_points: int,
    inp: Inputs,
    grid_dt: float,
    confirm_dt: float,
    calib_dt: float,
    target_nodes: list[tuple[float, float]],
    cache_grid: dict[tuple[float, float], tuple[np.ndarray, np.ndarray, np.ndarray, float]],
    cache_confirm: dict[tuple[float, float], tuple[np.ndarray, np.ndarray, np.ndarray, float]],
) -> tuple[dict[str, Any], pd.DataFrame]:
    summary, path = test_edge(
        edge, cand, n_points, "refined", inp, grid_dt, calib_dt, target_nodes, cache_grid, True
    )
    a = int(edge["sample_id_a"])
    b = int(edge["sample_id_b"])
    row_a, row_b = cand.loc[a], cand.loc[b]
    lam = path["lambda"].to_numpy(float)
    rates = log_interp(row_a[eng.RATE_COLS].to_numpy(float), row_b[eng.RATE_COLS].to_numpy(float), lam)

    r05 = np.full(len(path), np.nan)
    c05 = np.full(len(path), np.nan)
    b05 = np.full(len(path), np.nan)
    g05 = np.full(len(path), np.nan)
    t05 = np.full(len(path), np.nan)
    fallback_full = np.zeros(len(path), dtype=bool)

    # Confirm every refined path point at its dt=0.1 maximum node.
    for i in range(len(path)):
        g = float(path.iloc[i]["G_peak_support_node"])
        tau = float(path.iloc[i]["tau_G_ms_support_node"])
        g0, gm, g1, _ = cache_confirm[(g, tau)]
        out = eng.base_ratio_ensemble(rates[i : i + 1], g0, gm, g1, confirm_dt)[0]
        rr = float(out[2])
        best = (rr, float(out[0]), float(out[1]), g, tau)
        if not np.isfinite(rr) or rr < inp.threshold:
            # A dt-sensitive failure is rescanned over the complete frozen grid.
            fallback_full[i] = True
            rr_best = -np.inf
            best = (np.nan, np.nan, np.nan, np.nan, np.nan)
            for gg in eng.G_PEAK_GRID:
                for tt in eng.TAU_GRID:
                    g0, gm, g1, _ = cache_confirm[(float(gg), float(tt))]
                    oo = eng.base_ratio_ensemble(rates[i : i + 1], g0, gm, g1, confirm_dt)[0]
                    if np.isfinite(oo[2]) and float(oo[2]) > rr_best:
                        rr_best = float(oo[2])
                        best = (rr_best, float(oo[0]), float(oo[1]), float(gg), float(tt))
        r05[i], c05[i], b05[i], g05[i], t05[i] = best

    path["r_support_dt0p05"] = r05
    path["control_plateau_at_support_node_dt0p05"] = c05
    path["blocked_plateau_at_support_node_dt0p05"] = b05
    path["G_peak_support_node_dt0p05"] = g05
    path["tau_G_ms_support_node_dt0p05"] = t05
    path["dt0p05_full_grid_fallback"] = fallback_full
    strong05 = np.isfinite(r05) & (r05 >= inp.threshold)
    path["strong_dt0p05"] = strong05

    comp_mid = {
        "interp_aux": path["compatible_interp_aux_mid"].to_numpy(bool),
        "fixed_A_aux": path["compatible_fixed_A_aux_mid"].to_numpy(bool),
        "fixed_B_aux": path["compatible_fixed_B_aux_mid"].to_numpy(bool),
    }
    # Strict compatibility is recalculated directly from stored scores.
    comp_strict = {
        "interp_aux": path["score_interp_aux"].to_numpy(float) <= inp.cutoff_strict,
        "fixed_A_aux": path["score_fixed_A_aux"].to_numpy(float) <= inp.cutoff_strict,
        "fixed_B_aux": path["score_fixed_B_aux"].to_numpy(float) <= inp.cutoff_strict,
    }
    sched_mid = {s: bool(np.all(strong05 & comp_mid[s])) for s in comp_mid}
    sched_strict = {s: bool(np.all(strong05 & comp_strict[s])) for s in comp_strict}

    summary.update(
        {
            "strong_path_pass_dt0p05": bool(np.all(strong05)),
            "native_joint_pass_mid_dt0p05": sched_mid["interp_aux"],
            "native_joint_pass_strict_dt0p05": sched_strict["interp_aux"],
            "any_schedule_joint_pass_mid_dt0p05": any(sched_mid.values()),
            "any_schedule_joint_pass_strict_dt0p05": any(sched_strict.values()),
            "passing_schedule_mid_dt0p05": next((s for s in ["interp_aux", "fixed_A_aux", "fixed_B_aux"] if sched_mid[s]), ""),
            "passing_schedule_strict_dt0p05": next((s for s in ["interp_aux", "fixed_A_aux", "fixed_B_aux"] if sched_strict[s]), ""),
            "min_r_support_dt0p05": float(np.nanmin(r05)),
            "lambda_at_min_r_dt0p05": float(lam[int(np.nanargmin(r05))]),
            "max_abs_dt_ratio_change": float(np.nanmax(np.abs(r05 - path["r_support_dt0p1"].to_numpy(float)))),
            "n_dt0p05_full_grid_fallback": int(fallback_full.sum()),
            "min_control_plateau_at_support_node_dt0p05": float(np.nanmin(c05)),
            "median_control_plateau_at_support_node_dt0p05": float(np.nanmedian(c05)),
        }
    )
    return summary, path


def shortest_path(ids: list[int], edges: pd.DataFrame, pass_col: str, src: int, dst: int) -> tuple[float, list[int]]:
    adj: dict[int, list[tuple[float, int]]] = {x: [] for x in ids}
    for row in edges.itertuples(index=False):
        if bool(getattr(row, pass_col)):
            a, b, w = int(row.sample_id_a), int(row.sample_id_b), float(row.distance)
            adj[a].append((w, b))
            adj[b].append((w, a))
    pq = [(0.0, src)]
    dist = {src: 0.0}
    prev: dict[int, int] = {}
    while pq:
        d, u = heapq.heappop(pq)
        if d != dist.get(u):
            continue
        if u == dst:
            break
        for w, v in adj[u]:
            nd = d + w
            if nd < dist.get(v, float("inf")):
                dist[v] = nd
                prev[v] = u
                heapq.heappush(pq, (nd, v))
    if dst not in dist:
        return float("inf"), []
    path = [dst]
    while path[-1] != src:
        path.append(prev[path[-1]])
    path.reverse()
    return float(dist[dst]), path


def make_figures(
    inp: Inputs,
    refined_edges: pd.DataFrame,
    refined_points: pd.DataFrame,
    outdir: Path,
) -> None:
    coords = inp.pca.set_index("sample_id")
    ids = inp.confirmed["sample_id"].astype(int).tolist()

    fig, ax = plt.subplots(figsize=(8, 6.5))
    for row in refined_edges.itertuples(index=False):
        a, b = int(row.sample_id_a), int(row.sample_id_b)
        if a not in coords.index or b not in coords.index:
            continue
        x = [coords.loc[a, "PC1"], coords.loc[b, "PC1"]]
        y = [coords.loc[a, "PC2"], coords.loc[b, "PC2"]]
        if bool(getattr(row, "native_joint_pass_mid_dt0p05", False)):
            ax.plot(x, y, linewidth=1.1, alpha=0.75)
        elif bool(getattr(row, "any_schedule_joint_pass_mid_dt0p05", False)):
            ax.plot(x, y, linewidth=0.9, linestyle="--", alpha=0.65)
    ax.scatter(coords.loc[ids, "PC1"], coords.loc[ids, "PC2"], s=28, zorder=3)
    rep_ids = set(inp.reps["sample_id"].astype(int))
    for sid in rep_ids:
        if sid in coords.index:
            ax.text(coords.loc[sid, "PC1"], coords.loc[sid, "PC2"], str(sid), fontsize=7)
    ax.set_xlabel("PC1 of standardized log-rate coordinates")
    ax.set_ylabel("PC2")
    ax.set_title("Step 26 tested strong-compatible bridge network")
    fig.tight_layout()
    fig.savefig(outdir / "Fig26A_connectivity.pdf")
    fig.savefig(outdir / "Fig26A_connectivity.png", dpi=180)
    plt.close(fig)

    if len(refined_edges):
        fig, ax = plt.subplots(figsize=(7.5, 5.5))
        strong_margin = refined_edges["min_r_support_dt0p05"] - inp.threshold
        score_margin = inp.cutoff_mid - refined_edges[["max_score_interp_aux", "max_score_fixed_A_aux", "max_score_fixed_B_aux"]].min(axis=1)
        ax.scatter(strong_margin, score_margin, s=34)
        ax.axvline(0.0, linestyle="--", linewidth=1)
        ax.axhline(0.0, linestyle="--", linewidth=1)
        ax.set_xlabel("Minimum strong-response margin along refined edge")
        ax.set_ylabel("Best whole-edge calibration margin")
        ax.set_title("Refined bridge bottlenecks")
        fig.tight_layout()
        fig.savefig(outdir / "Fig26B_bridge_bottlenecks.pdf")
        fig.savefig(outdir / "Fig26B_bridge_bottlenecks.png", dpi=180)
        plt.close(fig)

    if len(refined_points):
        fig, ax = plt.subplots(figsize=(7.5, 5.5))
        ax.scatter(
            refined_points["control_plateau_at_support_node_dt0p05"],
            refined_points["r_support_dt0p05"],
            s=11,
            alpha=0.6,
        )
        ax.axhline(inp.threshold, linestyle="--", linewidth=1)
        ax.set_xscale("log")
        ax.set_xlabel("Control plateau at response-support node")
        ax.set_ylabel("Base route-removal ratio")
        ax.set_title("Denominator diagnostic along refined paths")
        fig.tight_layout()
        fig.savefig(outdir / "Fig26C_path_denominator.pdf")
        fig.savefig(outdir / "Fig26C_path_denominator.png", dpi=180)
        plt.close(fig)

        # Show direct tested representative-pair paths when available.
        rep_ids = set(inp.reps["sample_id"].astype(int))
        rep_edge_ids = []
        for r in refined_edges.itertuples(index=False):
            if int(r.sample_id_a) in rep_ids and int(r.sample_id_b) in rep_ids:
                rep_edge_ids.append(str(r.edge_id))
        rep_edge_ids = rep_edge_ids[:6]
        if rep_edge_ids:
            fig, ax = plt.subplots(figsize=(8, 5.8))
            for eid in rep_edge_ids:
                d = refined_points[refined_points["edge_id"] == eid]
                ax.plot(d["lambda"], d["r_support_dt0p05"] / inp.threshold, label=eid)
            ax.axhline(1.0, linestyle="--", linewidth=1)
            ax.set_xlabel("Path coordinate, lambda")
            ax.set_ylabel("r / frozen strong threshold")
            ax.set_title("Direct representative-pair bridge profiles")
            ax.legend(fontsize=7)
            fig.tight_layout()
            fig.savefig(outdir / "Fig26D_representative_paths.pdf")
            fig.savefig(outdir / "Fig26D_representative_paths.png", dpi=180)
            plt.close(fig)


def self_test() -> None:
    lam = np.linspace(0, 1, 5)
    a = np.array([1.0, 4.0])
    b = np.array([4.0, 1.0])
    x = log_interp(a, b, lam)
    assert np.allclose(x[0], a)
    assert np.allclose(x[-1], b)
    assert np.allclose(x[2], np.array([2.0, 2.0]))
    uf = UnionFind([1, 2, 3])
    assert uf.union(1, 2)
    assert len(uf.components()) == 2
    assert uf.union(2, 3)
    assert len(uf.components()) == 1
    print("Step26 self-test PASS")


def main() -> None:
    args = parse_args()
    if args.self_test:
        self_test()
        return
    if args.coarse_points < 3 or args.refine_points < args.coarse_points:
        raise ValueError("Require refine_points >= coarse_points >= 3")
    if args.max_edges < 39:
        raise ValueError("max_edges must be at least 39")
    set_num_threads(max(1, args.threads))

    step25 = discover_step25(args.step25)
    inp = load_inputs(step25, args.allow_step25_status)
    replay = replay_gate(inp, args.calibration_dt_ms, args.confirm_dt_ms)
    if not replay["pass"]:
        raise RuntimeError(f"Bundled Step25-engine replay gate failed: {replay}")

    outdir = args.output.expanduser().resolve()
    if outdir.exists():
        if not args.force:
            raise FileExistsError(f"Output exists: {outdir}; use --force to replace it")
        shutil.rmtree(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    input_audit = {
        "pipeline_step": "26_strong_compatible_path_connectivity",
        "script_version": SCRIPT_VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "step25": str(step25),
        "step25_status": inp.decision.get("status"),
        "confirmed_strong_n": int(len(inp.confirmed)),
        "strong_threshold": inp.threshold,
        "corrected_broad_cutoff_midpoint": inp.cutoff_mid,
        "corrected_broad_max_accepted_score": inp.cutoff_strict,
        "coarse_points_per_edge": args.coarse_points,
        "refine_points_per_edge": args.refine_points,
        "max_edges": args.max_edges,
        "knn_k": args.knn_k,
        "grid_dt_ms": args.grid_dt_ms,
        "confirmation_dt_ms": args.confirm_dt_ms,
        "calibration_dt_ms": args.calibration_dt_ms,
        "threads": get_num_threads(),
        "engine_replay_gate": replay,
        "no_new_denominator_cutoff": True,
        "path_definition": "geometric interpolation of 12 kinetic rates; three prespecified continuous auxiliary schedules",
        "decision_scope": "finite sampled-path connectivity, not a global topological proof",
    }
    json_dump(outdir / "00_input_audit.json", input_audit)
    inp.confirmed.to_csv(outdir / "01_strong_candidates.csv", index=False)

    edge_pool = build_edge_pool(inp, args.knn_k)
    edge_pool.to_csv(outdir / "02_edge_pool.csv", index=False)

    if args.preflight_only:
        json_dump(outdir / "preflight.json", {"status": "PASS", **input_audit})
        print(f"Step26 preflight PASS: {outdir}")
        return

    cand = inp.confirmed.set_index("sample_id", drop=False)
    ids = sorted(cand.index.astype(int).tolist())
    target_nodes = list(
        dict.fromkeys(
            (float(g), float(t))
            for g, t in inp.confirmed[["G_peak", "tau_G_ms"]].itertuples(index=False, name=None)
        )
    )
    cache_grid = forcing_cache(args.grid_dt_ms)
    cache_confirm = forcing_cache(args.confirm_dt_ms)

    coarse_summaries: list[dict[str, Any]] = []
    coarse_points: list[pd.DataFrame] = []
    max_edges = min(args.max_edges, len(edge_pool))
    for idx, edge in edge_pool.iloc[:max_edges].iterrows():
        summary, points = test_edge(
            edge,
            cand,
            args.coarse_points,
            "coarse",
            inp,
            args.grid_dt_ms,
            args.calibration_dt_ms,
            target_nodes,
            cache_grid,
            False,
        )
        summary["pool_rank"] = int(edge["pool_rank"])
        coarse_summaries.append(summary)
        coarse_points.append(points)
        if (len(coarse_summaries) % 10) == 0:
            print(f"Coarse edges tested: {len(coarse_summaries)}/{max_edges}", flush=True)

    coarse_edges = pd.DataFrame(coarse_summaries)
    coarse_pts = pd.concat(coarse_points, ignore_index=True)
    coarse_edges.to_csv(outdir / "03_coarse_edge_summary.csv", index=False)
    coarse_pts.to_csv(outdir / "04_coarse_path_points.csv.gz", index=False, compression="gzip")

    coarse_comp = pd.concat(
        [
            component_table(ids, coarse_edges, "strong_path_pass_dt0p1", "strong_only_dt0p1"),
            component_table(ids, coarse_edges, "native_joint_pass_mid_dt0p1", "native_joint_mid_dt0p1"),
            component_table(ids, coarse_edges, "any_schedule_joint_pass_mid_dt0p1", "any_schedule_joint_mid_dt0p1"),
        ],
        ignore_index=True,
    )
    coarse_comp.to_csv(outdir / "05_coarse_components.csv", index=False)

    # Select coarse passing spanning trees, then refine only edges needed to certify connectivity.
    wanted: set[str] = set()
    native_tree = kruskal_tree(ids, coarse_edges, "native_joint_pass_mid_dt0p1")
    any_tree = kruskal_tree(ids, coarse_edges, "any_schedule_joint_pass_mid_dt0p1")
    wanted.update(native_tree)
    wanted.update(any_tree)

    # If no full tree exists for a criterion, refine the passing edges of the Step25 MST and representative pairs
    # so the failure modes remain inspectable without refining every coarse edge.
    if not wanted:
        inspect = coarse_edges[
            coarse_edges["tags"].str.contains("step25_mst|representative_pair", regex=True)
            & coarse_edges["strong_path_pass_dt0p1"].astype(bool)
        ]
        wanted.update(inspect.sort_values("distance").head(60)["edge_id"].astype(str))

    refined_map: dict[str, dict[str, Any]] = {}
    refined_points_map: dict[str, pd.DataFrame] = {}

    def refine_edge_id(eid: str) -> None:
        if eid in refined_map:
            return
        erow = edge_pool[edge_pool["edge_id"] == eid].iloc[0]
        summary, points = exact_refine_with_dt_confirmation(
            erow,
            cand,
            args.refine_points,
            inp,
            args.grid_dt_ms,
            args.confirm_dt_ms,
            args.calibration_dt_ms,
            target_nodes,
            cache_grid,
            cache_confirm,
        )
        summary["pool_rank"] = int(erow["pool_rank"])
        refined_map[eid] = summary
        refined_points_map[eid] = points
        print(f"Refined edge {eid}", flush=True)

    for eid in sorted(wanted):
        refine_edge_id(eid)

    # If a coarse-connected criterion lost connectivity after refinement, refine alternative coarse-passing
    # edges that connect current refined components until connectivity is restored or candidates are exhausted.
    def restore_connectivity(coarse_pass_col: str, refined_pass_col: str) -> None:
        candidates = coarse_edges[coarse_edges[coarse_pass_col].astype(bool)].sort_values(["distance", "edge_id"])
        while True:
            redf = pd.DataFrame(refined_map.values()) if refined_map else pd.DataFrame()
            if len(redf) and is_connected(ids, redf, refined_pass_col):
                return
            # Determine current components under refined passing edges.
            uf = UnionFind(ids)
            if len(redf):
                for r in redf.itertuples(index=False):
                    if bool(getattr(r, refined_pass_col)):
                        uf.union(int(r.sample_id_a), int(r.sample_id_b))
            next_eid = None
            for r in candidates.itertuples(index=False):
                eid = str(r.edge_id)
                if eid in refined_map:
                    continue
                if uf.find(int(r.sample_id_a)) != uf.find(int(r.sample_id_b)):
                    next_eid = eid
                    break
            if next_eid is None:
                return
            refine_edge_id(next_eid)

    if native_tree:
        restore_connectivity("native_joint_pass_mid_dt0p1", "native_joint_pass_mid_dt0p05")
    if any_tree:
        restore_connectivity("any_schedule_joint_pass_mid_dt0p1", "any_schedule_joint_pass_mid_dt0p05")

    refined_edges = pd.DataFrame(refined_map.values()) if refined_map else pd.DataFrame()
    refined_pts = pd.concat(refined_points_map.values(), ignore_index=True) if refined_points_map else pd.DataFrame()
    refined_edges.to_csv(outdir / "06_refined_edge_summary.csv", index=False)
    if len(refined_pts):
        refined_pts.to_csv(outdir / "07_refined_path_points.csv.gz", index=False, compression="gzip")
    else:
        pd.DataFrame().to_csv(outdir / "07_refined_path_points.csv.gz", index=False, compression="gzip")

    if len(refined_edges):
        refined_comp = pd.concat(
            [
                component_table(ids, refined_edges, "strong_path_pass_dt0p05", "strong_only_dt0p05"),
                component_table(ids, refined_edges, "native_joint_pass_mid_dt0p05", "native_joint_mid_dt0p05"),
                component_table(ids, refined_edges, "any_schedule_joint_pass_mid_dt0p05", "any_schedule_joint_mid_dt0p05"),
                component_table(ids, refined_edges, "native_joint_pass_strict_dt0p05", "native_joint_strict_dt0p05"),
                component_table(ids, refined_edges, "any_schedule_joint_pass_strict_dt0p05", "any_schedule_joint_strict_dt0p05"),
            ],
            ignore_index=True,
        )
    else:
        refined_comp = pd.DataFrame(columns=["criterion", "component_id", "component_size", "sample_id"])
    refined_comp.to_csv(outdir / "08_refined_components.csv", index=False)

    native_connected = bool(len(refined_edges) and is_connected(ids, refined_edges, "native_joint_pass_mid_dt0p05"))
    any_connected = bool(len(refined_edges) and is_connected(ids, refined_edges, "any_schedule_joint_pass_mid_dt0p05"))
    strong_connected = bool(len(refined_edges) and is_connected(ids, refined_edges, "strong_path_pass_dt0p05"))
    native_strict_connected = bool(len(refined_edges) and is_connected(ids, refined_edges, "native_joint_pass_strict_dt0p05"))
    any_strict_connected = bool(len(refined_edges) and is_connected(ids, refined_edges, "any_schedule_joint_pass_strict_dt0p05"))

    # Representative-to-representative connectivity within the refined tested graph.
    rep_ids = [int(x) for x in inp.reps["sample_id"].tolist() if int(x) in ids]
    rep_rows = []
    for i, a in enumerate(rep_ids):
        for b in rep_ids[i + 1 :]:
            for crit, col in [
                ("native_mid", "native_joint_pass_mid_dt0p05"),
                ("any_schedule_mid", "any_schedule_joint_pass_mid_dt0p05"),
            ]:
                if len(refined_edges):
                    dist, path = shortest_path(ids, refined_edges, col, a, b)
                else:
                    dist, path = float("inf"), []
                rep_rows.append(
                    {
                        "criterion": crit,
                        "sample_id_a": a,
                        "sample_id_b": b,
                        "connected": bool(path),
                        "path_distance": dist if np.isfinite(dist) else np.nan,
                        "path_sample_ids": "->".join(map(str, path)),
                        "n_edges": max(0, len(path) - 1),
                    }
                )
    pd.DataFrame(rep_rows).to_csv(outdir / "09_representative_connectivity.csv", index=False)

    if len(refined_edges):
        bottleneck = refined_edges[
            [
                "edge_id", "sample_id_a", "sample_id_b", "distance",
                "min_r_support_dt0p05", "lambda_at_min_r_dt0p05",
                "max_score_interp_aux", "max_score_fixed_A_aux", "max_score_fixed_B_aux",
                "passing_schedule_mid_dt0p05", "min_control_plateau_at_support_node_dt0p05",
                "native_joint_pass_mid_dt0p05", "any_schedule_joint_pass_mid_dt0p05",
            ]
        ].copy()
        bottleneck["strong_margin"] = bottleneck["min_r_support_dt0p05"] - inp.threshold
        bottleneck["best_whole_edge_score"] = bottleneck[["max_score_interp_aux", "max_score_fixed_A_aux", "max_score_fixed_B_aux"]].min(axis=1)
        bottleneck["best_whole_edge_calibration_margin"] = inp.cutoff_mid - bottleneck["best_whole_edge_score"]
    else:
        bottleneck = pd.DataFrame()
    bottleneck.to_csv(outdir / "10_bottleneck_summary.csv", index=False)

    ncomp_native = int(refined_comp.loc[refined_comp["criterion"] == "native_joint_mid_dt0p05", "component_id"].nunique()) if len(refined_comp) else len(ids)
    ncomp_any = int(refined_comp.loc[refined_comp["criterion"] == "any_schedule_joint_mid_dt0p05", "component_id"].nunique()) if len(refined_comp) else len(ids)
    ncomp_strong = int(refined_comp.loc[refined_comp["criterion"] == "strong_only_dt0p05", "component_id"].nunique()) if len(refined_comp) else len(ids)

    if native_connected:
        status = "ONE_SAMPLED_CONNECTED_REGION_SUPPORTED_NATIVE_INTERPOLATION"
    elif any_connected:
        status = "ONE_SAMPLED_KINETIC_REGION_SUPPORTED_WITH_PRESPECIFIED_AUX_SCHEDULES"
    else:
        status = "SAMPLED_CONNECTIVITY_NOT_ESTABLISHED"

    decision = {
        "status": status,
        "interpretation_scope": (
            "Connectivity of the 40 Step25 strong corrected-compatible vectors along explicit finite-resolution "
            "log-linear paths. This is not a proof of global basin topology or physiological population structure."
        ),
        "n_strong_candidates": len(ids),
        "coarse_edges_tested": int(len(coarse_edges)),
        "refined_edges_tested": int(len(refined_edges)),
        "strong_only_connected_refined": strong_connected,
        "native_midpoint_cutoff_connected_refined": native_connected,
        "any_prespecified_aux_schedule_midpoint_cutoff_connected_refined": any_connected,
        "native_strict_cutoff_connected_refined": native_strict_connected,
        "any_prespecified_aux_schedule_strict_cutoff_connected_refined": any_strict_connected,
        "n_components_strong_only_refined": ncomp_strong,
        "n_components_native_mid_refined": ncomp_native,
        "n_components_any_schedule_mid_refined": ncomp_any,
        "midpoint_cutoff": inp.cutoff_mid,
        "strict_cutoff": inp.cutoff_strict,
        "strong_threshold": inp.threshold,
        "no_denominator_exclusion_applied": True,
        "engine_replay_gate_pass": replay["pass"],
    }
    json_dump(outdir / "11_scientific_decision.json", decision)

    make_figures(inp, refined_edges, refined_pts, outdir)

    readme = [
        "# Step 26 strong-compatible path connectivity",
        "",
        f"Step25 input status: **{inp.decision.get('status')}**.",
        f"Strong confirmed candidates: **{len(ids)}**.",
        f"Coarse tested edges: **{len(coarse_edges)}**; refined edges: **{len(refined_edges)}**.",
        f"Strong-only refined components: **{ncomp_strong}**.",
        f"Native interpolated-aux refined components: **{ncomp_native}**.",
        f"Any prespecified whole-edge auxiliary schedule refined components: **{ncomp_any}**.",
        f"Decision: **{status}**.",
        "",
        "No new denominator cutoff was used. Control plateaus along refined paths are retained as diagnostics.",
        "A passing edge requires every sampled path point to remain strong at dt=0.05 ms and corrected-control compatible under one entire prespecified auxiliary schedule.",
        "The result is finite sampled-path evidence, not a global proof of connectedness or disconnectedness.",
    ]
    (outdir / "README_RESULTS.md").write_text("\n".join(readme) + "\n", encoding="utf-8")

    run_summary = {
        "pipeline_step": "26_strong_compatible_path_connectivity",
        "script_version": SCRIPT_VERSION,
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "numba_threads": get_num_threads(),
        "step25": str(step25),
        "step25_decision": inp.decision,
        "engine_replay": replay,
        "strong_threshold": inp.threshold,
        "cutoff_midpoint": inp.cutoff_mid,
        "cutoff_strict": inp.cutoff_strict,
        "coarse_points": args.coarse_points,
        "refine_points": args.refine_points,
        "max_edges": args.max_edges,
        "knn_k": args.knn_k,
        "grid_dt_ms": args.grid_dt_ms,
        "confirm_dt_ms": args.confirm_dt_ms,
        "calibration_dt_ms": args.calibration_dt_ms,
        "targeted_forcing_nodes_in_coarse_screen": len(target_nodes),
        "decision": decision,
    }
    json_dump(outdir / "run_summary.json", run_summary)
    print(f"Step26 results written to: {outdir}", flush=True)


if __name__ == "__main__":
    main()
