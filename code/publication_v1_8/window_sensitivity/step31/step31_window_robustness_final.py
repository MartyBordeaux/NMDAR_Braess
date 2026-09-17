#!/usr/bin/env python3
"""NMDAR Braess Step 31: final analysis-window robustness + manuscript geometry.

Production goals
----------------
1. Re-extract the five frozen experimental control series from raw ABF files and
   verify that the baseline 140-200 / 140-250 / 200-500 / 140-500 definitions
   reproduce the frozen calibration targets.
2. Vary the three otherwise operational response-window boundaries on a
   prespecified 3x3x3 grid while keeping the post-train start fixed at 140 ms.
3. Recompute experiment-derived calibration targets/tolerances and replay the
   control calibration for all 100,000 Step-10 candidates using one model
   integration per candidate plus cumulative-window readout.
4. Quantify accepted-cohort stability and retention of the final 40
   experiment-scale amplification solutions without introducing a new response
   threshold or denominator exclusion.
5. Produce a manuscript-facing parameter-space geometry figure and a
   non-Bayesian sampling->calibration->experiment-scale marginal summary.

This pipeline contains no historical-error narrative. It treats the frozen
140-500-ms implementation as the baseline analysis and asks only whether the
conclusions depend strongly on nearby operational window choices.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import math
import os
import re
import shutil
import sys
import time
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd

try:
    from numba import njit, prange, set_num_threads
except Exception as exc:  # pragma: no cover
    raise SystemExit("numba is required for Step31") from exc

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

SCRIPT_VERSION = "1.0.2"
PULSE_COUNT = 25
PULSE_INTERVAL_MS = 5.1
STRONG_THRESHOLD = 1.6359295439505144
BASELINE_START_MS = 140.0
SPLIT_VALUES = np.array([180.0, 200.0, 220.0])
PRIMARY_END_VALUES = np.array([230.0, 250.0, 270.0])
LATE_END_VALUES = np.array([450.0, 500.0, 550.0])
MODEL_BOUNDARIES_MS = np.array([140.0, 180.0, 200.0, 220.0, 230.0, 250.0, 270.0, 450.0, 500.0, 550.0])
BASELINE_SCHEME = (200.0, 250.0, 500.0)

RATE_COLS = [
    "kon_A_per_ms", "kon_B_per_ms", "koff_A_per_ms", "koff_B_per_ms",
    "kf_plus_per_ms", "ks_plus_per_ms", "kf_minus_per_ms", "ks_minus_per_ms",
    "kd1_plus_per_ms", "kd1_minus_per_ms", "kd2_plus_per_ms", "kd2_minus_per_ms",
]
AUX_COLS = ["pulse_amplitude_au", "glutamate_tau_ms"]

# Frozen calibration series-level ratios used only as a hard replay check.
# They are not fitted or altered by this pipeline.
BASELINE_SERIES_EXPECTED = {
    "2012?": None,  # sentinel never used; prevents accidental positional logic
    "12.10.2022": (1.167, 0.615, 254.440),
    "17.10.2022 cell1": (2.110, 0.396, 7.539),
    "02.12.2021 cell1": (0.948, 0.795, 295.305),
    "17.10.2022 cell2": (1.100, 0.638, 257.220),
    "02.12.2021 cell2": (0.972, 1.002, 358.844),
}
BASELINE_SERIES_EXPECTED.pop("2012?")

BASELINE_TARGET_EXPECTED = np.array([1.099894, 0.637550, 257.219945], dtype=float)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Step31 final window robustness and parameter-space geometry")
    p.add_argument("--step10", default="auto")
    p.add_argument("--step24", default="auto")
    p.add_argument("--step29", default="auto")
    p.add_argument("--step27", default="auto")
    p.add_argument("--step28", default="auto")
    p.add_argument("--step2", default="auto", help="Step-2 results archive/directory containing slow_metrics_files.csv")
    p.add_argument("--raw-root", default="auto", help="Raw ABF root")
    p.add_argument("--output", type=Path, default=Path("/root/nmda2/step_31/results_step_31_window_robustness_final"))
    p.add_argument("--threads", type=int, default=int(os.environ.get("STEP31_THREADS", "2")))
    p.add_argument("--chunk-size", type=int, default=int(os.environ.get("STEP31_CHUNK_SIZE", "5000")))
    p.add_argument("--model-dt-ms", type=float, default=0.05)
    p.add_argument("--geometry-only", action="store_true")
    p.add_argument("--force", action="store_true")
    p.add_argument("--resume", action="store_true")
    p.add_argument("--self-test", action="store_true")
    p.add_argument("--preflight-only", action="store_true")
    return p.parse_args()


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_json(path: Path, obj: Any) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False, sort_keys=True) + "\n", encoding="utf-8")


def discover_dir(spec: str, candidates: list[Path], marker: str) -> Path:
    if spec != "auto":
        p = Path(spec).expanduser().resolve()
        if not p.exists():
            raise FileNotFoundError(p)
        return p
    for p in candidates:
        if p.is_dir() and (p / marker).exists():
            return p.resolve()
    for root in [Path("/root/nmda2"), Path("/root/nmda")]:
        if root.exists():
            for m in root.rglob(marker):
                return m.parent.resolve()
    raise FileNotFoundError(f"Could not auto-discover directory containing {marker}")


def discover_step24(spec: str) -> Path:
    return discover_dir(spec, [
        Path("/root/nmda2/step_24/results_step_24_calibration_window_correction"),
        Path("/root/nmda2/results_step_24_calibration_window_correction"),
    ], "02_corrected_candidate_scores.csv.gz")


def discover_step29(spec: str) -> Path:
    return discover_dir(spec, [
        Path("/root/nmda2/step_29/results_step_29_corrected_publication_consolidation"),
        Path("/root/nmda2/results_step_29_corrected_publication_consolidation"),
    ], "02_final_strong40_topology.csv")


def discover_optional_step(spec: str, candidates: list[Path], marker: str) -> Path | None:
    if spec != "auto":
        p = Path(spec).expanduser().resolve()
        if not p.exists():
            raise FileNotFoundError(p)
        return p
    for p in candidates:
        if p.is_dir() and (p / marker).exists():
            return p.resolve()
    return None


def discover_step10(spec: str, step24: Path) -> tuple[Path, str | None]:
    if spec != "auto":
        p = Path(spec).expanduser().resolve()
        if not p.exists():
            raise FileNotFoundError(p)
        return p, None
    audit_path = step24 / "00_input_audit.json"
    if audit_path.exists():
        audit = json.loads(audit_path.read_text(encoding="utf-8"))
        src = str(audit.get("step10_source", ""))
        if "::" in src:
            left, member = src.split("::", 1)
            pp = Path(left)
            if pp.exists():
                return pp.resolve(), member
    for p in [Path("/root/nmda/results_step10.zip"), Path("/root/nmda/results_step10_tables.zip")]:
        if p.exists():
            return p.resolve(), None
    raise FileNotFoundError("Could not auto-discover results_step10.zip")


def read_step10(path: Path, preferred_member: str | None) -> tuple[pd.DataFrame, str]:
    cols = ["prior", "sample_id", *RATE_COLS, *AUX_COLS]
    if path.suffix.lower() == ".zip":
        with zipfile.ZipFile(path) as zf:
            names = zf.namelist()
            member = preferred_member if preferred_member in names else None
            if member is None:
                matches = [n for n in names if n.endswith("matched_parameter_all_samples.csv")]
                if not matches:
                    matches = [n for n in names if n.endswith(".csv") and "all_samples" in n]
                if not matches:
                    raise FileNotFoundError(f"No full candidate table inside {path}")
                member = sorted(matches, key=len)[0]
            with zf.open(member) as f:
                df = pd.read_csv(f, usecols=cols)
        return df, f"{path}::{member}"
    if path.is_file():
        return pd.read_csv(path, usecols=cols), str(path)
    matches = list(path.rglob("matched_parameter_all_samples.csv"))
    if not matches:
        raise FileNotFoundError(f"No matched_parameter_all_samples.csv below {path}")
    pp = sorted(matches, key=lambda x: len(str(x)))[0]
    return pd.read_csv(pp, usecols=cols), str(pp)


def discover_step2(spec: str) -> Path:
    if spec != "auto":
        p = Path(spec).expanduser().resolve()
        if not p.exists():
            raise FileNotFoundError(p)
        return p
    candidates = [
        Path("/root/nmda/iv_nmda_pipeline_step2/results_step2_tables.zip"),
        Path("/root/nmda/IV_NMDA/iv_nmda_pipeline_step2/results_step2_tables.zip"),
        Path("/root/nmda/results_step2_tables.zip"),
    ]
    for p in candidates:
        if p.exists():
            return p.resolve()
    for root in [Path("/root/nmda"), Path("/root/nmda2")]:
        if root.exists():
            ms = list(root.rglob("results_step2_tables.zip"))
            if ms:
                return sorted(ms, key=lambda x: len(str(x)))[0].resolve()
    raise FileNotFoundError("Could not auto-discover results_step2_tables.zip")


def read_slow_metrics(step2: Path) -> tuple[pd.DataFrame, str]:
    if step2.suffix.lower() == ".zip":
        with zipfile.ZipFile(step2) as zf:
            ms = [n for n in zf.namelist() if n.endswith("slow_metrics_files.csv")]
            if not ms:
                raise FileNotFoundError(f"slow_metrics_files.csv not found inside {step2}")
            member = sorted(ms, key=len)[0]
            with zf.open(member) as f:
                return pd.read_csv(f), f"{step2}::{member}"
    if step2.is_file():
        return pd.read_csv(step2), str(step2)
    ms = list(step2.rglob("slow_metrics_files.csv"))
    if not ms:
        raise FileNotFoundError(f"slow_metrics_files.csv not found below {step2}")
    pp = sorted(ms, key=lambda x: len(str(x)))[0]
    return pd.read_csv(pp), str(pp)


def discover_raw_root(spec: str) -> Path:
    if spec != "auto":
        p = Path(spec).expanduser().resolve()
        if not p.is_dir():
            raise FileNotFoundError(p)
        return p
    for p in [Path("/root/nmda/IV_NMDA"), Path("/home/vlad/nmda/db_nmda")]:
        if p.is_dir():
            return p.resolve()
    raise FileNotFoundError("Could not auto-discover raw ABF root")


def norm_text(x: Any) -> str:
    s = str(x).lower().replace("ё", "е")
    return re.sub(r"[^a-zа-я0-9]+", "", s)


def row_blob(row: pd.Series) -> str:
    vals = []
    for c, v in row.items():
        if pd.isna(v):
            continue
        if any(k in c.lower() for k in ["series", "date", "cell", "file", "path", "name", "folder", "condition", "drug", "nbqx"]):
            vals.append(str(v))
    return " | ".join(vals)


def date_tokens(label: str) -> list[str]:
    m = re.search(r"(\d{2})\.(\d{2})\.(\d{4})", label)
    if not m:
        return [norm_text(label)]
    d, mo, y = m.groups()
    return [f"{d}{mo}{y}", f"{y}{mo}{d}"]


def alias_match(blob: str, alias: str) -> bool:
    nb = norm_text(blob)
    toks = date_tokens(alias)
    if not any(t in nb for t in toks):
        return False
    cm = re.search(r"cell\s*([0-9]+)", alias.lower())
    if cm and f"cell{cm.group(1)}" not in nb:
        return False
    return True


def detect_voltage_mask(df: pd.DataFrame) -> np.ndarray:
    candidates = [c for c in df.columns if "volt" in c.lower()]
    for c in candidates:
        vals = pd.to_numeric(df[c], errors="coerce")
        if vals.notna().sum() > 0:
            # support V or mV encodings
            arr = vals.to_numpy(float)
            if np.nanmedian(np.abs(arr[np.isfinite(arr)])) < 2:
                return np.isclose(arr, -0.04, atol=0.005)
            return np.isclose(arr, -40.0, atol=2.0)
    # last resort: textual -40 marker
    return np.array([("-40" in row_blob(r)) for _, r in df.iterrows()], dtype=bool)


def detect_nbqx_mask(df: pd.DataFrame) -> np.ndarray:
    cols = [c for c in df.columns if "nbqx" in c.lower()]
    if not cols:
        return np.ones(len(df), dtype=bool)
    c = cols[0]
    v = df[c]
    if pd.api.types.is_bool_dtype(v):
        return v.fillna(False).to_numpy(bool)
    s = v.astype(str).str.lower().str.strip()
    return s.isin(["true", "1", "yes", "y", "nbqx"]).to_numpy(bool)


def file_value_from_row(row: pd.Series) -> str | None:
    preferred = [c for c in row.index if any(k in c.lower() for k in ["abf", "file", "path"])]
    for c in preferred:
        v = row[c]
        if pd.isna(v):
            continue
        s = str(v)
        if ".abf" in s.lower():
            m = re.search(r"([^/\\]+\.abf)", s, flags=re.I)
            return m.group(1) if m else s
    for v in row.values:
        if isinstance(v, str) and ".abf" in v.lower():
            m = re.search(r"([^/\\]+\.abf)", v, flags=re.I)
            return m.group(1) if m else v
    return None


def build_abf_index(root: Path) -> dict[str, list[Path]]:
    out: dict[str, list[Path]] = {}
    for p in root.rglob("*.abf"):
        out.setdefault(p.name.lower(), []).append(p)
    return out


def _local_path_components(h: Path, n_parents: int = 3) -> list[str]:
    """Return basename plus the nearest protocol-level path components.

    The archive itself may live below a generic directory literally named
    ``Ro25``.  Looking at the full absolute path therefore produces false drug
    positives.  Conversely, some recordings have an extra voltage/subfolder
    below the protocol folder, so checking only the immediate parent misses the
    actual condition.  The basename plus the three nearest parents covers both
    layouts while keeping the archive-level ``Ro25`` directory out of normal
    date/cell/protocol paths.
    """
    parts = [h.name]
    q = h.parent
    for _ in range(n_parents):
        if q == q.parent:
            break
        parts.append(q.name)
        q = q.parent
    return parts


def _condition_is_ro25_path(h: Path) -> bool:
    for comp in _local_path_components(h):
        raw = comp.lower().replace(" ", "")
        norm = norm_text(comp)
        # A component exactly named Ro25 is an archive/container name, not a
        # treatment label.  Protocol folders such as SlowEPSC_IV+NBQX+Ro25 are
        # treatment labels and must be excluded.
        if norm == "ro25":
            continue
        if ("ro25" in norm) or ("ro-25" in raw):
            return True
    return False


def _condition_is_mem_path(h: Path) -> bool:
    for comp in _local_path_components(h):
        norm = norm_text(comp)
        if ("memantine" in norm) or ("memslow" in norm) or ("nbqxmem" in norm):
            return True
    return False


def _is_slowepsc_path(h: Path) -> bool:
    # Protocol folders may be the immediate parent or one level higher (e.g.
    # protocol/voltage/file).  Check the same local path context used for the
    # treatment filters.
    for comp in _local_path_components(h):
        norm = norm_text(comp)
        if "slowepsc" in norm or ("slow" in norm and "epsc" in norm):
            return True
    return False


def _alias_path_match(h: Path, alias: str) -> bool:
    hp = norm_text(str(h))
    if not any(tok in hp for tok in date_tokens(alias)):
        return False
    cm = re.search(r"cell\s*([0-9]+)", alias.lower())
    if cm and f"cell{cm.group(1)}" not in hp:
        return False
    return True


def locate_abf(name_or_path: str, raw_root: Path, index: dict[str, list[Path]], row: pd.Series, alias: str | None = None, require_reference: bool = False) -> Path:
    p = Path(name_or_path)
    base = p.name.lower()

    # Never trust an existing path blindly for calibration. Step-2 tables can
    # contain paths to a condition file (Ro25/memantine) or to another protocol
    # that shares the same basename. Existing paths must pass the same identity
    # and condition/protocol filters as basename-resolved candidates.
    direct = p.resolve() if p.exists() else None
    if direct is not None:
        direct_ok = True
        if alias and not _alias_path_match(direct, alias):
            direct_ok = False
        if require_reference:
            if _condition_is_ro25_path(direct) or _condition_is_mem_path(direct):
                direct_ok = False
            if not _is_slowepsc_path(direct):
                direct_ok = False
        if direct_ok:
            return direct

    hits = list(index.get(base, []))
    if not hits:
        raise FileNotFoundError(f"Raw ABF {p.name!r} not found under {raw_root}")

    # First use the frozen calibration-series identity (date and, where present,
    # explicit cell number).  This is substantially safer than basename matching
    # because identical ABF basenames recur in different cells.
    if alias:
        ah = [h for h in hits if _alias_path_match(h, alias)]
        if ah:
            hits = ah

    # Calibration targets are derived from the reference (+NBQX, no +Ro25)
    # recordings.  Filter on the immediate condition folder/basename only; the
    # archive root may itself be called Ro25.
    if require_reference:
        rh = [h for h in hits if (not _condition_is_ro25_path(h)) and (not _condition_is_mem_path(h))]
        if rh:
            hits = rh
        else:
            raise RuntimeError(
                f"No untreated reference candidate remained for {base!r}, alias {alias!r}; "
                "all basename matches were Ro25 or memantine conditions."
            )

        # Step 2 slow-metric calibration must come from SlowEPSC recordings.
        # Identical basenames can also occur in PF_IV+NBQX directories within
        # the same date/cell. Prefer the SlowEPSC protocol explicitly rather
        # than resolving this by an arbitrary basename tie-break.
        sh = [h for h in hits if _is_slowepsc_path(h)]
        if sh:
            hits = sh
        else:
            raise RuntimeError(
                f"No SlowEPSC untreated reference candidate remained for {base!r}, alias {alias!r}; "
                f"candidates after treatment filtering were {[str(h) for h in hits[:8]]}"
            )

    if len(hits) == 1:
        return hits[0]

    blob = norm_text(row_blob(row))
    scored = []
    for h in hits:
        hp = norm_text(str(h.parent))
        score = sum(1 for token in re.findall(r"[a-zа-я0-9]{4,}", blob) if token in hp)
        scored.append((score, h))
    scored.sort(key=lambda z: (-z[0], len(str(z[1]))))
    if len(scored) >= 2 and scored[0][0] == scored[1][0]:
        hs = [sha256_file(x[1]) for x in scored[:2]]
        if hs[0] != hs[1]:
            raise RuntimeError(
                f"Ambiguous non-identical ABF basename {base} for alias {alias!r}: "
                f"{[str(x[1]) for x in scored[:6]]}"
            )
    return scored[0][1]


def current_to_nA(y: np.ndarray, unit: str) -> np.ndarray:
    u = unit.lower().replace("µ", "u").strip()
    if "pa" in u:
        return y.astype(float) / 1000.0
    if "na" in u:
        return y.astype(float)
    if "ua" in u:
        return y.astype(float) * 1000.0
    if u == "a" or u.endswith(" amp"):
        return y.astype(float) * 1e9
    raise RuntimeError(f"Unsupported current unit {unit!r}")


def select_current_channel(abf) -> int:
    units = [str(x) for x in getattr(abf, "adcUnits", [])]
    names = [str(x) for x in getattr(abf, "adcNames", [])]
    scores = []
    for i in range(getattr(abf, "channelCount", len(units))):
        unit = units[i] if i < len(units) else ""
        name = names[i] if i < len(names) else ""
        s = 0
        lu = unit.lower().replace("µ", "u")
        if any(k in lu for k in ["pa", "na", "ua"]): s += 10
        if any(k in name.lower() for k in ["current", "im", "amp", "i"]): s += 2
        scores.append((s, i, unit, name))
    scores.sort(reverse=True)
    if not scores or scores[0][0] < 10:
        raise RuntimeError(f"No current-like ADC channel found; channels={list(zip(names, units))}")
    return scores[0][1]


def detrend_sweep(t_s: np.ndarray, i_nA: np.ndarray) -> np.ndarray:
    m0 = (t_s >= 0.0) & (t_s < 0.008)
    mf = (t_s >= 4.5) & (t_s <= 5.0)
    if m0.sum() < 3 or mf.sum() < 3:
        raise RuntimeError("Trace does not contain required baseline windows 0-8 ms and 4.5-5.0 s")
    b0 = float(np.mean(i_nA[m0]))
    bf = float(np.mean(i_nA[mf]))
    baseline = b0 + (bf - b0) / (4.750 - 0.004) * (t_s - 0.004)
    return i_nA - baseline


def window_mean(t_ms: np.ndarray, y: np.ndarray, a: float, b: float, closed_right: bool = False) -> float:
    if closed_right:
        m = (t_ms >= a) & (t_ms <= b)
    else:
        m = (t_ms >= a) & (t_ms < b)
    if m.sum() < 2:
        return float("nan")
    return float(np.mean(y[m]))


def window_charge(t_s: np.ndarray, y_nA: np.ndarray, a_ms: float, b_ms: float, method: str) -> float:
    t_ms = t_s * 1000.0
    if method == "rect_halfopen":
        m = (t_ms >= a_ms) & (t_ms < b_ms)
        if m.sum() < 2: return float("nan")
        dt_s = float(np.median(np.diff(t_s[m])))
        return float(np.sum(y_nA[m]) * dt_s * 1000.0)  # nA*s -> pC
    if method == "trap_halfopen":
        m = (t_ms >= a_ms) & (t_ms < b_ms)
    elif method == "trap_closed":
        m = (t_ms >= a_ms) & (t_ms <= b_ms)
    else:
        raise ValueError(method)
    if m.sum() < 2: return float("nan")
    return float(np.trapezoid(y_nA[m], t_s[m]) * 1000.0)


def read_abf_sweep_cache(path: Path) -> list[tuple[np.ndarray, np.ndarray]]:
    try:
        import pyabf
    except Exception as exc:  # pragma: no cover
        raise RuntimeError("pyabf is required for experimental window sensitivity") from exc
    abf = pyabf.ABF(str(path))
    ch = select_current_channel(abf)
    unit = str(abf.adcUnits[ch])
    out = []
    for sw in abf.sweepList:
        abf.setSweep(int(sw), channel=ch)
        t = np.asarray(abf.sweepX, dtype=float)
        y = current_to_nA(np.asarray(abf.sweepY, dtype=float), unit)
        out.append((t, detrend_sweep(t, y)))
    if not out:
        raise RuntimeError(f"No sweeps in {path}")
    return out


def file_metrics_from_cache(cache: list[tuple[np.ndarray, np.ndarray]], split_ms: float, primary_end_ms: float,
                            late_end_ms: float, integration_method: str) -> dict[str, float]:
    rows = []
    for t_s, yd in cache:
        t_ms = t_s * 1000.0
        e = window_mean(t_ms, yd, 140.0, split_ms)
        h = window_mean(t_ms, yd, 140.0, primary_end_ms)
        l = window_mean(t_ms, yd, split_ms, late_end_ms)
        q = window_charge(t_s, yd, 140.0, late_end_ms, integration_method)
        rows.append((e, h, l, q))
    a = np.asarray(rows, dtype=float)
    return {"E_nA": float(np.mean(a[:,0])), "H_nA": float(np.mean(a[:,1])), "L_nA": float(np.mean(a[:,2])), "Q_pC": float(np.mean(a[:,3])), "n_sweeps": int(len(rows))}


def select_calibration_rows(slow: pd.DataFrame) -> dict[str, pd.DataFrame]:
    vmask = detect_voltage_mask(slow)
    nmask = detect_nbqx_mask(slow)
    out = {}
    series_cols = [c for c in slow.columns if c.lower() in {"series_id", "series", "cell_series"} or "series_id" in c.lower()]
    for alias in BASELINE_SERIES_EXPECTED:
        exact = np.zeros(len(slow), dtype=bool)
        # Prefer an explicit series identifier. This prevents a date-only alias from
        # accidentally absorbing multiple cell series recorded on the same day.
        for c in series_cols:
            vals = slow[c].astype(str)
            for ii, v in enumerate(vals):
                nv = norm_text(v); na = norm_text(alias)
                same = (nv == na)
                if not same:
                    # Allow YYYYMMDD versus DDMMYYYY date order while still requiring
                    # the same explicit cell token when one is present.
                    same = alias_match(v, alias) and (re.search(r"cell\s*[0-9]+", alias.lower()) is not None or "cell" not in nv)
                if same: exact[ii] = True
        if exact.any():
            mask = exact & vmask & nmask
        else:
            mask = np.array([alias_match(row_blob(r), alias) for _, r in slow.iterrows()], dtype=bool) & vmask & nmask
        sub = slow.loc[mask].copy()
        if sub.empty:
            raise RuntimeError(f"No Step-2 rows matched calibration series {alias!r} at -40 mV")
        # Keep the reference recording, not the +Ro25 recording.  Do this from
        # the actual ABF filename/path field rather than row_blob(): metadata may
        # legitimately contain the drug name even for the paired reference row.
        keep = []
        for _, r in sub.iterrows():
            fv = file_value_from_row(r)
            leaf = (fv or "").lower().replace(" ", "")
            is_drug = ("ro25" in leaf) or ("ro-25" in leaf) or ("+mem" in leaf) or ("memantine" in leaf)
            keep.append(not is_drug)
        if any(keep):
            sub = sub.loc[np.asarray(keep, dtype=bool)].copy()
        else:
            raise RuntimeError(
                f"No untreated reference Step-2 rows remained for calibration series {alias!r}. "
                "Refusing to substitute Ro25 or memantine rows."
            )
        out[alias] = sub
    return out


def choose_integration_method(series_files: dict[str, list[tuple[pd.Series, Path, list[tuple[np.ndarray, np.ndarray]]]]]) -> tuple[str, pd.DataFrame]:
    methods = ["rect_halfopen", "trap_halfopen", "trap_closed"]
    audits = []
    for method in methods:
        err = 0.0
        rows = []
        for alias, entries in series_files.items():
            fm = [file_metrics_from_cache(cache, 200.0, 250.0, 500.0, method) for _, _, cache in entries]
            e = float(np.median([x["E_nA"] for x in fm])); h = float(np.median([x["H_nA"] for x in fm]))
            l = float(np.median([x["L_nA"] for x in fm])); q = float(np.median([x["Q_pC"] for x in fm]))
            ratios = np.array([abs(e)/abs(h), abs(l)/abs(h), abs(q)/abs(h)], dtype=float)
            exp = np.asarray(BASELINE_SERIES_EXPECTED[alias], dtype=float)
            rel = np.abs(np.log(ratios / exp))
            err += float(np.sum(rel*rel))
            rows.append({"method":method,"series":alias,"E_over_H":ratios[0],"L_over_H":ratios[1],"Q_over_H_ms":ratios[2],"expected_E_over_H":exp[0],"expected_L_over_H":exp[1],"expected_Q_over_H_ms":exp[2],"log_error_norm":float(np.linalg.norm(rel))})
        audits.extend(rows)
        audits.append({"method":method,"series":"__TOTAL__","E_over_H":np.nan,"L_over_H":np.nan,"Q_over_H_ms":np.nan,"expected_E_over_H":np.nan,"expected_L_over_H":np.nan,"expected_Q_over_H_ms":np.nan,"log_error_norm":math.sqrt(err)})
    audit = pd.DataFrame(audits)
    totals = audit[audit.series == "__TOTAL__"].sort_values("log_error_norm")
    method = str(totals.iloc[0].method)
    return method, audit


def robust_tol(vals: np.ndarray, floor: float) -> float:
    x = np.log(np.asarray(vals, dtype=float))
    q25, q75 = np.quantile(x, [0.25,0.75])
    return float(max(floor, (q75-q25)/1.349))


def experimental_targets_for_scheme(series_files, split_ms, primary_end_ms, late_end_ms, method) -> tuple[dict[str,float], list[dict[str,Any]]]:
    rows=[]
    er=[]; lr=[]; qr=[]
    for alias, entries in series_files.items():
        fm=[file_metrics_from_cache(cache,split_ms,primary_end_ms,late_end_ms,method) for _,_,cache in entries]
        e=float(np.median([x["E_nA"] for x in fm])); h=float(np.median([x["H_nA"] for x in fm])); l=float(np.median([x["L_nA"] for x in fm])); q=float(np.median([x["Q_pC"] for x in fm]))
        a,b,c=abs(e)/abs(h),abs(l)/abs(h),abs(q)/abs(h)
        er.append(a); lr.append(b); qr.append(c)
        rows.append({"series":alias,"E_nA":e,"H_nA":h,"L_nA":l,"Q_pC":q,"E_over_H":a,"L_over_H":b,"Q_over_H_ms":c,"n_files":len(fm),"n_sweeps_total":sum(int(x["n_sweeps"]) for x in fm)})
    target={
        "target_E_over_H":float(np.median(er)),
        "target_L_over_H":float(np.median(lr)),
        "target_Q_over_H_ms":float(np.median(qr)),
        "tol_E":robust_tol(np.asarray(er),0.35),
        "tol_L":robust_tol(np.asarray(lr),0.35),
        "tol_Q":robust_tol(np.asarray(qr),0.50),
    }
    return target,rows


def anchored_experimental_targets_for_scheme(series_files, split_ms, primary_end_ms, late_end_ms, method, baseline_raw_ratios):
    """Return raw and baseline-anchored series ratios for one window scheme.

    The five frozen series-level baseline ratios define the absolute calibration
    reference. Raw ABFs supply only the within-series change produced by moving
    a window boundary. This prevents condition/file-layout details or tiny
    implementation differences in raw replay from redefining the frozen baseline.
    """
    raw_target, raw_rows = experimental_targets_for_scheme(
        series_files, split_ms, primary_end_ms, late_end_ms, method
    )
    anchored_rows=[]; er=[]; lr=[]; qr=[]
    for rr in raw_rows:
        alias=rr["series"]
        base=baseline_raw_ratios[alias]
        exp=np.asarray(BASELINE_SERIES_EXPECTED[alias], dtype=float)
        raw=np.asarray([rr["E_over_H"], rr["L_over_H"], rr["Q_over_H_ms"]], dtype=float)
        den=np.asarray(base, dtype=float)
        if np.any(~np.isfinite(raw)) or np.any(~np.isfinite(den)) or np.any(raw<=0) or np.any(den<=0):
            raise RuntimeError(f"Invalid raw ratio for experimental anchoring: {alias}")
        anc=exp*(raw/den)
        er.append(float(anc[0])); lr.append(float(anc[1])); qr.append(float(anc[2]))
        anchored_rows.append({
            **rr,
            "raw_E_over_H": float(raw[0]), "raw_L_over_H": float(raw[1]), "raw_Q_over_H_ms": float(raw[2]),
            "anchored_E_over_H": float(anc[0]), "anchored_L_over_H": float(anc[1]), "anchored_Q_over_H_ms": float(anc[2]),
            "frozen_baseline_E_over_H": float(exp[0]), "frozen_baseline_L_over_H": float(exp[1]), "frozen_baseline_Q_over_H_ms": float(exp[2]),
        })
    target={
        "target_E_over_H":float(np.median(er)),
        "target_L_over_H":float(np.median(lr)),
        "target_Q_over_H_ms":float(np.median(qr)),
        "tol_E":robust_tol(np.asarray(er),0.35),
        "tol_L":robust_tol(np.asarray(lr),0.35),
        "tol_Q":robust_tol(np.asarray(qr),0.50),
        "raw_target_E_over_H":float(raw_target["target_E_over_H"]),
        "raw_target_L_over_H":float(raw_target["target_L_over_H"]),
        "raw_target_Q_over_H_ms":float(raw_target["target_Q_over_H_ms"]),
    }
    return target, anchored_rows


def recover_frozen_targets_from_scores(scores: pd.DataFrame) -> dict[str, float]:
    h=scores["control_plateau_PO_140_250"].to_numpy(float)
    e=scores["control_early_PO_140_200"].to_numpy(float)
    l=scores["control_late_PO_200_500"].to_numpy(float)
    j=scores["control_integrated_open_state_140_500_ms"].to_numpy(float)
    ss=scores["calibration_score_corrected"].to_numpy(float)
    m=np.isfinite(h)&np.isfinite(e)&np.isfinite(l)&np.isfinite(j)&np.isfinite(ss)&(h>0)&(e>0)&(l>0)&(j>0)
    x=np.column_stack([np.log(e[m]/h[m]),np.log(l[m]/h[m]),np.log(j[m]/h[m])])
    w=np.array([1/0.35**2,1/0.35**2,0.5/0.50**2],float)
    y=ss[m]-np.sum(w*x*x,axis=1)
    beta=np.linalg.lstsq(np.column_stack([x,np.ones(len(x))]),y,rcond=None)[0]
    mu=-beta[:3]/(2*w)
    replay=np.sum(w*(x-mu)**2,axis=1)
    err=float(np.max(np.abs(replay-ss[m])))
    if err>1e-8: raise RuntimeError(f"Could not reconstruct frozen calibration targets: {err:g}")
    return {"target_E_over_H":float(np.exp(mu[0])),"target_L_over_H":float(np.exp(mu[1])),"target_Q_over_H_ms":float(np.exp(mu[2])),"tol_E":0.35,"tol_L":0.35,"tol_Q":0.50,"score_reconstruction_max_error":err}


# ----------------------- model calibration replay -------------------------
@njit(cache=True)
def _calib_rhs(state: np.ndarray, g: float, rates: np.ndarray) -> np.ndarray:
    p0,p1a,p1b,p2,pd1,pd2,pc1,pc2,po,pro=state
    kon_a,kon_b,koff_a,koff_b=rates[0],rates[1],rates[2],rates[3]
    kfp,ksp,kfm,ksm=rates[4],rates[5],rates[6],rates[7]
    kd1p,kd1m,kd2p,kd2m=rates[8],rates[9],rates[10],rates[11]
    f0a=kon_a*g*p0; f0b=kon_b*g*p0; r_a=koff_a*p1a; r_b=koff_b*p1b
    f_ab=kon_b*g*p1a; r_2a=koff_b*p2; f_ba=kon_a*g*p1b; r_2b=koff_a*p2
    f_d1=kd1p*p2; r_d1=kd1m*pd1; f_d2=kd2p*p2; r_d2=kd2m*pd2
    f_c1=kfp*p2; r_c1=kfm*pc1; f_c2=ksp*p2; r_c2=ksm*pc2
    c1_o=ksp*pc1; o_c1=ksm*po; c2_o=kfp*pc2; o_c2=kfm*po
    out=np.empty(10,np.float64)
    out[0]=-f0a-f0b+r_a+r_b
    out[1]=f0a-r_a-f_ab+r_2a
    out[2]=f0b-r_b-f_ba+r_2b
    out[3]=f_ab+f_ba-r_2a-r_2b-f_d1-f_d2-f_c1-f_c2+r_d1+r_d2+r_c1+r_c2
    out[4]=f_d1-r_d1; out[5]=f_d2-r_d2
    out[6]=f_c1-r_c1-c1_o+o_c1; out[7]=f_c2-r_c2-c2_o+o_c2
    out[8]=c1_o+c2_o-o_c1-o_c2; out[9]=0.0
    return out

@njit(cache=True)
def _prefix_one(rates: np.ndarray, pulse_amp: float, tau_ms: float, dt_ms: float, boundaries: np.ndarray) -> np.ndarray:
    duration_ms=600.0
    steps=int(round(duration_ms/dt_ms))
    state=np.zeros(10,np.float64); state[0]=1.0
    g=0.0; decay=math.exp(-dt_ms/tau_ms); next_pulse=0
    nb=len(boundaries)
    prefix=np.zeros(nb,np.float64)
    cum=0.0
    bi=0
    mass_error=0.0; minimum_state=0.0
    for step in range(steps+1):
        t=step*dt_ms
        while next_pulse<PULSE_COUNT and t+0.5*dt_ms>=next_pulse*PULSE_INTERVAL_MS:
            g+=pulse_amp; next_pulse+=1
        # prefix at boundary contains samples t < boundary
        while bi<nb and t>=boundaries[bi]-1e-12:
            prefix[bi]=cum
            bi+=1
        po=state[8]
        if step<steps:
            cum += po*dt_ms
        mass=0.0
        for k in range(10):
            mass+=state[k]
            if state[k]<minimum_state: minimum_state=state[k]
        err=abs(mass-1.0)
        if err>mass_error: mass_error=err
        if step==steps: break
        g_eff=g/(1.0+g)
        d1=_calib_rhs(state,g_eff,rates)
        midpoint=state+0.5*dt_ms*d1
        g_next=g*decay; g_next_eff=g_next/(1.0+g_next)
        dmid=_calib_rhs(midpoint,0.5*(g_eff+g_next_eff),rates)
        state=state+dt_ms*dmid; g=g_next
    while bi<nb:
        prefix[bi]=cum; bi+=1
    out=np.empty(nb+3,np.float64)
    out[:nb]=prefix
    out[nb]=mass_error; out[nb+1]=minimum_state; out[nb+2]=1.0 if mass_error<1e-6 and minimum_state>-1e-7 else 0.0
    return out

@njit(parallel=True,cache=True)
def _prefix_chunk(rates: np.ndarray, amps: np.ndarray, taus: np.ndarray, dt_ms: float, boundaries: np.ndarray) -> np.ndarray:
    n=rates.shape[0]; out=np.empty((n,len(boundaries)+3),np.float64)
    for i in prange(n):
        out[i,:]=_prefix_one(rates[i],amps[i],taus[i],dt_ms,boundaries)
    return out


def model_prefix_all(full: pd.DataFrame, outdir: Path, dt_ms: float, chunk_size: int, resume: bool) -> pd.DataFrame:
    cdir=outdir/"_checkpoints"/"model_prefix_chunks"
    cdir.mkdir(parents=True,exist_ok=True)
    frames=[]; n=len(full)
    t0=time.time()
    for start in range(0,n,chunk_size):
        stop=min(n,start+chunk_size)
        fp=cdir/f"chunk_{start:06d}_{stop:06d}.csv.gz"
        if resume and fp.exists():
            q=pd.read_csv(fp)
        else:
            sub=full.iloc[start:stop]
            arr=_prefix_chunk(sub[RATE_COLS].to_numpy(float),sub["pulse_amplitude_au"].to_numpy(float),sub["glutamate_tau_ms"].to_numpy(float),dt_ms,MODEL_BOUNDARIES_MS)
            q=sub[["prior","sample_id"]].reset_index(drop=True).copy()
            for j,b in enumerate(MODEL_BOUNDARIES_MS): q[f"prefix_PO_ms_{int(b)}"]=arr[:,j]
            q["mass_error"]=arr[:,len(MODEL_BOUNDARIES_MS)]
            q["minimum_state"]=arr[:,len(MODEL_BOUNDARIES_MS)+1]
            q["integration_ok"]=arr[:,len(MODEL_BOUNDARIES_MS)+2]>0.5
            q.to_csv(fp,index=False,compression="gzip")
        frames.append(q)
        elapsed=(time.time()-t0)/3600
        print(f"Model prefix replay: {stop}/{n} ({100*stop/n:.1f}%) elapsed={elapsed:.2f} h",flush=True)
    return pd.concat(frames,ignore_index=True)


def prefix_col(t: float) -> str: return f"prefix_PO_ms_{int(t)}"


def model_metrics(prefix: pd.DataFrame, split_ms: float, primary_end_ms: float, late_end_ms: float) -> tuple[np.ndarray,np.ndarray,np.ndarray,np.ndarray]:
    p140=prefix[prefix_col(140.0)].to_numpy(float)
    ps=prefix[prefix_col(split_ms)].to_numpy(float)
    pp=prefix[prefix_col(primary_end_ms)].to_numpy(float)
    pl=prefix[prefix_col(late_end_ms)].to_numpy(float)
    e=(ps-p140)/(split_ms-140.0)
    h=(pp-p140)/(primary_end_ms-140.0)
    l=(pl-ps)/(late_end_ms-split_ms)
    j=pl-p140
    return h,e,l,j


def score_from_metrics(h,e,l,j,target):
    out=np.full(len(h),np.inf,float)
    v=np.isfinite(h)&np.isfinite(e)&np.isfinite(l)&np.isfinite(j)&(h>0)&(e>0)&(l>0)&(j>0)
    er=e[v]/h[v]; lr=l[v]/h[v]; qr=j[v]/h[v]
    out[v]=(np.log(er/target["target_E_over_H"])/target["tol_E"])**2+(np.log(lr/target["target_L_over_H"])/target["tol_L"])**2+0.5*(np.log(qr/target["target_Q_over_H_ms"])/target["tol_Q"])**2
    return out


def scheme_id(split,primary,late): return f"split{int(split)}_primary{int(primary)}_late{int(late)}"


def jaccard(a:set,b:set)->float:
    return len(a&b)/len(a|b) if (a|b) else float("nan")


def make_geometry(step29: Path, full: pd.DataFrame, step27: Path|None, step28: Path|None, outdir: Path) -> dict[str,Any]:
    accepted=pd.read_csv(step29/"01_frozen_corrected_accepted_cohorts.csv.gz")
    broad=accepted[accepted.prior.astype(str).str.lower()=="broad"].copy()
    strong=pd.read_csv(step29/"02_final_strong40_topology.csv")
    edges=pd.read_csv(step29/"03_final_verified_network_edges.csv")
    reps=pd.read_csv(step29/"05_publication_representatives.csv")
    x=np.log(broad[RATE_COLS].to_numpy(float))
    sc=StandardScaler(); z=sc.fit_transform(x); pca=PCA(n_components=2,random_state=0); xy=pca.fit_transform(z)
    coords=broad[["sample_id"]].copy(); coords["PC1"]=xy[:,0]; coords["PC2"]=xy[:,1]; coords["is_experiment_scale"]=coords.sample_id.isin(set(strong.sample_id.astype(int)))
    regime_map=strong.set_index("sample_id")["final_regime"].to_dict(); coords["final_regime"]=coords.sample_id.map(regime_map).fillna("control_compatible_other")
    coords.to_csv(outdir/"07_parameter_space_coordinates.csv.gz",index=False,compression="gzip")
    edges.to_csv(outdir/"08_geometry_edges.csv",index=False)

    pathdf=None; path_source=None
    candidates=[]
    if step27 and (step27/"06_verified_path_points.csv.gz").exists(): candidates.append(step27/"06_verified_path_points.csv.gz")
    if step28 and (step28/"06_verified_path_points.csv.gz").exists(): candidates.append(step28/"06_verified_path_points.csv.gz")
    for p in candidates:
        d=pd.read_csv(p)
        # Prefer path from extreme experimental match to primary typical representative.
        m=((d.sample_id_a==27084)&(d.sample_id_b==30902))|((d.sample_id_a==30902)&(d.sample_id_b==27084))
        if m.any(): pathdf=d.loc[m].copy(); path_source=str(p); break
    if pathdf is None and candidates:
        d=pd.read_csv(candidates[0]); key=d[["sample_id_a","sample_id_b","mode"]].drop_duplicates().iloc[0]
        m=(d.sample_id_a==key.sample_id_a)&(d.sample_id_b==key.sample_id_b)&(d["mode"]==key["mode"])
        pathdf=d.loc[m].copy(); path_source=str(candidates[0])
    if pathdf is not None:
        pathdf.to_csv(outdir/"09_representative_path_profile.csv",index=False)

    # Figure: four manuscript-facing panels.
    fig,axs=plt.subplots(2,2,figsize=(11.0,8.5))
    ax=axs[0,0]
    ax.scatter(coords.PC1,coords.PC2,s=6,alpha=.18,label="control-compatible models")
    reg_order=["dominant_native_connected","dominant_flex_only_accession","secondary_native_connected","historical_witness_singleton"]
    labels={"dominant_native_connected":"dominant connected region","dominant_flex_only_accession":"aux-flex accession","secondary_native_connected":"smaller connected region","historical_witness_singleton":"isolated sampled solution"}
    markers={"dominant_native_connected":"o","dominant_flex_only_accession":"s","secondary_native_connected":"^","historical_witness_singleton":"X"}
    cmap={}
    for reg in reg_order:
        s=coords[coords.final_regime==reg]
        if len(s): ax.scatter(s.PC1,s.PC2,s=42,marker=markers[reg],label=labels[reg])
    ax.set_xlabel("PC1 of standardized log-rates"); ax.set_ylabel("PC2")
    ax.set_title("A  Control-compatible parameter space")
    ax.legend(frameon=False,fontsize=7,loc="best")

    ax=axs[0,1]
    cxy=coords.set_index("sample_id")[["PC1","PC2"]]
    for _,e in edges.iterrows():
        a=int(e.sample_id_a); b=int(e.sample_id_b)
        if a in cxy.index and b in cxy.index:
            ax.plot([cxy.loc[a,"PC1"],cxy.loc[b,"PC1"]],[cxy.loc[a,"PC2"],cxy.loc[b,"PC2"]],lw=.7,alpha=.35,zorder=1,color="0.45")
    for reg in reg_order:
        s=coords[coords.final_regime==reg]
        if len(s): ax.scatter(s.PC1,s.PC2,s=45,marker=markers[reg],label=labels[reg],zorder=2)
    ax.set_xlabel("PC1"); ax.set_ylabel("PC2"); ax.set_title("B  Verified connectivity of experiment-scale solutions")
    ax.legend(frameon=False,fontsize=7,loc="best")

    ax=axs[1,0]
    for reg in reg_order:
        s=strong[strong.final_regime==reg]
        if len(s): ax.scatter(s.calibration_score_corrected,s.r_plateau_dt0p05,s=45,marker=markers[reg],label=labels[reg])
    ax.axhline(STRONG_THRESHOLD,ls="--",lw=1,label="experiment-scale threshold")
    cutoff=float(strong.calibration_score_corrected.max())
    # true calibration cutoff from step29 headline if available
    head=pd.read_csv(step29/"07_manuscript_headline_numbers.csv")
    rr=head[head.quantity.astype(str).str.contains("cutoff",case=False,regex=True)]
    if len(rr): cutoff=float(rr.iloc[0].value)
    ax.axvline(cutoff,ls=":",lw=1,label="control-compatibility cutoff")
    ax.set_xlabel("control-calibration score"); ax.set_ylabel("maximum amplification ratio")
    ax.set_yscale("log"); ax.set_title("C  Compatibility and amplification")
    ax.legend(frameon=False,fontsize=7,loc="best")

    ax=axs[1,1]
    if pathdf is not None:
        p=pathdf.sort_values("lambda")
        xlam=p["lambda"].to_numpy(float)
        s=p["calibration_score"].to_numpy(float); r=p["r_plateau_dt0p05"].to_numpy(float)
        ax.plot(xlam,s,label="calibration score")
        ax.axhline(cutoff,ls=":",lw=1)
        ax.set_xlabel("position along verified path"); ax.set_ylabel("calibration score")
        ax2=ax.twinx(); ax2.plot(xlam,r,ls="--",label="amplification ratio"); ax2.axhline(STRONG_THRESHOLD,ls=":",lw=1); ax2.set_ylabel("amplification ratio")
        ax.set_title("D  Meaning of a verified path")
        lines=ax.get_lines()[:1]+ax2.get_lines()[:1]; ax.legend(lines,[l.get_label() for l in lines],frameon=False,fontsize=7,loc="best")
    else:
        ax.text(.5,.5,"Verified path-point file not available",ha="center",va="center",transform=ax.transAxes); ax.set_axis_off()
    fig.tight_layout()
    fig.savefig(outdir/"Fig31B_parameter_space_geometry.pdf",bbox_inches="tight")
    fig.savefig(outdir/"Fig31B_parameter_space_geometry.png",dpi=220,bbox_inches="tight")
    plt.close(fig)

    # Non-Bayesian marginal filtering figure.
    pre=full[full.prior.astype(str).str.lower()=="broad"].copy()
    strong_ids=set(strong.sample_id.astype(int))
    fig,axs=plt.subplots(3,4,figsize=(12,8.5))
    for ax,col in zip(axs.flat,RATE_COLS):
        arrays=[np.log10(pre[col].to_numpy(float)),np.log10(broad[col].to_numpy(float)),np.log10(strong[col].to_numpy(float))]
        lo=min(np.nanmin(a) for a in arrays); hi=max(np.nanmax(a) for a in arrays); bins=np.linspace(lo,hi,30)
        for arr,lab,lw in [(arrays[0],"sampled broad",1.0),(arrays[1],"control-compatible",1.2)]:
            hist,ed=np.histogram(arr,bins=bins,density=True); cen=.5*(ed[:-1]+ed[1:]); ax.plot(cen,hist,label=lab,lw=lw)
        ymax=ax.get_ylim()[1]
        ax.scatter(arrays[2], np.full_like(arrays[2], -0.035*ymax), marker="|", s=36, label="experiment-scale (n=40)", clip_on=False)
        ax.set_ylim(bottom=-0.08*ymax)
        ax.set_title(col.replace("_per_ms",""),fontsize=9); ax.set_xlabel("log10 rate",fontsize=7); ax.tick_params(labelsize=7)
    handles,labels2=axs.flat[0].get_legend_handles_labels(); fig.legend(handles,labels2,loc="upper center",ncol=3,frameon=False)
    fig.suptitle("Kinetic-rate filtering across sampled, control-compatible, and experiment-scale sets",y=.995)
    fig.tight_layout(rect=[0,0,.99,.96])
    fig.savefig(outdir/"Fig31C_rate_filtering_marginals.pdf",bbox_inches="tight")
    fig.savefig(outdir/"Fig31C_rate_filtering_marginals.png",dpi=220,bbox_inches="tight")
    plt.close(fig)

    return {"pca_explained_variance_ratio":pca.explained_variance_ratio_.tolist(),"path_source":path_source,"n_broad_accepted":int(len(broad)),"n_experiment_scale":int(len(strong))}


def plot_window_sensitivity(summary: pd.DataFrame, targets: pd.DataFrame, outdir: Path) -> None:
    b=summary[(summary.prior=="broad")].copy()
    fig,axs=plt.subplots(1,3,figsize=(12,3.6))
    ax=axs[0]
    ax.scatter(np.arange(len(b)),b.jaccard_vs_baseline,s=22)
    ax.axhline(1,ls=":",lw=1); ax.set_ylim(0,1.03); ax.set_xlabel("prespecified window scheme"); ax.set_ylabel("Jaccard vs baseline top-5000"); ax.set_title("A  Calibration-cohort stability")
    ax=axs[1]
    ax.scatter(np.arange(len(b)),b.final40_retained_n,s=22,label="all 40")
    ax.scatter(np.arange(len(b)),b.dominant36_retained_n,s=22,label="dominant 36")
    ax.set_xlabel("prespecified window scheme"); ax.set_ylabel("retained experiment-scale solutions"); ax.set_ylim(0,42); ax.set_title("B  Core-solution retention"); ax.legend(frameon=False,fontsize=7)
    ax=axs[2]
    base=targets[targets.is_baseline].iloc[0]
    ax.scatter(np.arange(len(targets)),targets.target_E_over_H/base.target_E_over_H,s=18,label="E/H")
    ax.scatter(np.arange(len(targets)),targets.target_L_over_H/base.target_L_over_H,s=18,label="L/H")
    ax.scatter(np.arange(len(targets)),targets.target_Q_over_H_ms/base.target_Q_over_H_ms,s=18,label="Q/H")
    ax.axhline(1,ls=":",lw=1); ax.set_xlabel("prespecified window scheme"); ax.set_ylabel("target / baseline target"); ax.set_title("C  Experimental target stability"); ax.legend(frameon=False,fontsize=7)
    fig.tight_layout(); fig.savefig(outdir/"Fig31A_window_sensitivity.pdf",bbox_inches="tight"); fig.savefig(outdir/"Fig31A_window_sensitivity.png",dpi=220,bbox_inches="tight"); plt.close(fig)


def run_self_test() -> None:
    rates=np.array([0.4,0.3,0.1,0.1,1.0,0.8,0.3,0.2,0.05,0.01,0.03,0.01],float)
    a=_prefix_one(rates,0.2,10.0,0.1,MODEL_BOUNDARIES_MS)
    assert np.all(np.isfinite(a)) and a[-1] in (0.0,1.0)
    h=np.array([1.0]);e=np.array([1.1]);l=np.array([.64]);j=np.array([257.0])
    t={"target_E_over_H":1.1,"target_L_over_H":.64,"target_Q_over_H_ms":257.0,"tol_E":.35,"tol_L":.35,"tol_Q":.5}
    assert abs(score_from_metrics(h,e,l,j,t)[0])<1e-12
    assert scheme_id(200,250,500)=="split200_primary250_late500"
    print("Step31 self-test PASS")


def main() -> None:
    args=parse_args()
    if args.self_test:
        run_self_test(); return
    set_num_threads(max(1,args.threads))
    out=args.output.expanduser().resolve()
    if out.exists():
        if args.force:
            shutil.rmtree(out)
        elif not args.resume:
            raise FileExistsError(f"Output exists: {out}. Set STEP31_RESUME=1 or STEP31_FORCE=1.")
    out.mkdir(parents=True,exist_ok=True)

    step24=discover_step24(args.step24); step29=discover_step29(args.step29)
    step27=discover_optional_step(args.step27,[Path("/root/nmda2/step_27/results_step_27_curved_bridge_optimization")],"06_verified_path_points.csv.gz")
    step28=discover_optional_step(args.step28,[Path("/root/nmda2/step_28/results_step_28_residual_connectivity_stress")],"06_verified_path_points.csv.gz")
    step10_path,member=discover_step10(args.step10,step24); full,step10_source=read_step10(step10_path,member)
    scores24=pd.read_csv(step24/"02_corrected_candidate_scores.csv.gz")
    frozen_target=recover_frozen_targets_from_scores(scores24)
    frozen_obs_cols=[
        "control_plateau_PO_140_250",
        "control_early_PO_140_200",
        "control_late_PO_200_500",
        "control_integrated_open_state_140_500_ms",
    ]
    full=full.merge(
        scores24[["prior","sample_id","calibration_score_corrected","accepted_corrected",*frozen_obs_cols]],
        on=["prior","sample_id"],how="left",validate="one_to_one"
    )
    if full.calibration_score_corrected.isna().any(): raise RuntimeError("Step10-Step24 join incomplete")

    audit={"script_version":SCRIPT_VERSION,"step10_source":step10_source,"step24":str(step24),"step29":str(step29),"step27":str(step27) if step27 else None,"step28":str(step28) if step28 else None,"threads":args.threads,"model_dt_ms":args.model_dt_ms,"frozen_target_reconstructed_from_step24":frozen_target,"window_grid":{"post_train_start_ms":140.0,"split_ms":SPLIT_VALUES.tolist(),"primary_end_ms":PRIMARY_END_VALUES.tolist(),"late_end_ms":LATE_END_VALUES.tolist()}}
    geom=make_geometry(step29,full,step27,step28,out); audit["geometry"]=geom
    if args.geometry_only:
        write_json(out/"00_input_audit.json",audit); print(f"Geometry-only Step31 complete: {out}"); return

    step2=discover_step2(args.step2); raw_root=discover_raw_root(args.raw_root); slow,slow_source=read_slow_metrics(step2)
    audit["step2_source"]=slow_source; audit["raw_root"]=str(raw_root)
    selected=select_calibration_rows(slow); idx=build_abf_index(raw_root)
    series_files={}; mapping=[]
    cache_by_path={}
    for alias,sub in selected.items():
        ent=[]
        for ridx,row in sub.iterrows():
            fn=file_value_from_row(row)
            if fn is None: raise RuntimeError(f"No ABF filename in selected Step2 row {ridx} for {alias}")
            p=locate_abf(fn,raw_root,idx,row,alias=alias,require_reference=True)
            if p not in cache_by_path: cache_by_path[p]=read_abf_sweep_cache(p)
            ent.append((row,p,cache_by_path[p]))
            mapping.append({"series":alias,"step2_row_index":int(ridx),"abf_basename":p.name,"abf_path":str(p),"n_sweeps":len(cache_by_path[p]),"sha256":sha256_file(p),"is_ro25_path":_condition_is_ro25_path(p),"is_mem_path":_condition_is_mem_path(p),"is_slowepsc_path":_is_slowepsc_path(p),"local_path_context":" | ".join(_local_path_components(p))})
        # deduplicate exact same ABF rows
        seen=set(); uniq=[]
        for e in ent:
            if e[1] not in seen: uniq.append(e); seen.add(e[1])
        series_files[alias]=uniq
    mapdf=pd.DataFrame(mapping)
    mapdf.to_csv(out/"01_experimental_file_mapping.csv",index=False)
    badmap = mapdf.loc[
        mapdf["is_ro25_path"].astype(bool)
        | mapdf["is_mem_path"].astype(bool)
        | (~mapdf["is_slowepsc_path"].astype(bool))
    ].copy() if not mapdf.empty else mapdf.copy()
    badmap.to_csv(out/"01A_mapping_gate_diagnostics.csv", index=False)
    if mapdf.empty or not badmap.empty:
        audit["mapping_gate_offending_rows"] = int(len(badmap))
        audit["mapping_gate_offending_paths"] = badmap["abf_path"].astype(str).tolist() if not badmap.empty else []
        write_json(out/"00_input_audit.json",audit)
        details = "; ".join(
            f"{r.series}: {r.abf_path} [ro25={r.is_ro25_path}, mem={r.is_mem_path}, slow={r.is_slowepsc_path}]"
            for r in badmap.itertuples(index=False)
        )
        raise RuntimeError(
            "Experimental mapping gate failed after deterministic resolution. "
            f"Offending rows: {details or 'mapping table is empty'}"
        )

    integ_method,replay_methods=choose_integration_method(series_files); replay_methods.to_csv(out/"02A_integration_method_replay.csv",index=False)
    raw_base_target,raw_base_rows=experimental_targets_for_scheme(series_files,200,250,500,integ_method)
    raw_bdf=pd.DataFrame(raw_base_rows)
    replay=[]; baseline_raw_ratios={}
    for rr in raw_base_rows:
        alias=rr["series"]; exp=np.asarray(BASELINE_SERIES_EXPECTED[alias],float)
        raw=np.asarray([rr["E_over_H"],rr["L_over_H"],rr["Q_over_H_ms"]],float)
        baseline_raw_ratios[alias]=raw.tolist()
        replay.append({**rr,"expected_E_over_H":exp[0],"expected_L_over_H":exp[1],"expected_Q_over_H_ms":exp[2],
                       "raw_logerr_E":abs(float(np.log(raw[0]/exp[0]))),"raw_logerr_L":abs(float(np.log(raw[1]/exp[1]))),"raw_logerr_Q":abs(float(np.log(raw[2]/exp[2])))})
    replay_df=pd.DataFrame(replay)
    replay_df.to_csv(out/"02_baseline_series_raw_replay.csv",index=False)
    audit["baseline_raw_series_max_abs_log_error_vs_frozen"] = float(
        replay_df[["raw_logerr_E","raw_logerr_L","raw_logerr_Q"]].to_numpy(float).max()
    )

    # Hard five-series baseline replay after anchoring: every series must equal its
    # frozen baseline ratio, not merely reproduce the overall median target.
    anchored_base_target,anchored_base_rows=anchored_experimental_targets_for_scheme(series_files,200,250,500,integ_method,baseline_raw_ratios)
    abr=pd.DataFrame(anchored_base_rows)
    abr.to_csv(out/"02B_baseline_series_anchored_replay.csv",index=False)
    per_series_err=[]
    for _,rr in abr.iterrows():
        exp=np.asarray(BASELINE_SERIES_EXPECTED[rr["series"]],float)
        anc=np.asarray([rr["anchored_E_over_H"],rr["anchored_L_over_H"],rr["anchored_Q_over_H_ms"]],float)
        per_series_err.extend(np.abs(np.log(anc/exp)).tolist())
    max_series_anchor_err=float(np.max(per_series_err)) if per_series_err else float("inf")
    audit["integration_method"]=integ_method
    audit["baseline_raw_experimental_target"]=raw_base_target
    audit["baseline_anchored_experimental_target"]=anchored_base_target
    audit["baseline_five_series_anchor_max_abs_log_error"]=max_series_anchor_err
    audit["experimental_anchor_policy"]="frozen series-level baseline ratios x raw within-series window-change ratio"
    if max_series_anchor_err>1e-12:
        write_json(out/"00_input_audit.json",audit)
        raise RuntimeError(f"Five-series anchored baseline replay failed: max abs log error {max_series_anchor_err:g}")

    schemes=[]; target_rows=[]; series_rows=[]
    for sp in SPLIT_VALUES:
        for pe in PRIMARY_END_VALUES:
            for le in LATE_END_VALUES:
                if not (140<sp<pe<le): continue
                sid=scheme_id(sp,pe,le); isbase=(sp,pe,le)==BASELINE_SCHEME
                tar,srows=anchored_experimental_targets_for_scheme(series_files,float(sp),float(pe),float(le),integ_method,baseline_raw_ratios)
                schemes.append({"scheme_id":sid,"split_ms":sp,"primary_end_ms":pe,"late_end_ms":le,"E_width_ms":sp-140,"H_width_ms":pe-140,"L_width_ms":le-sp,"is_baseline":isbase})
                target_rows.append({"scheme_id":sid,"split_ms":sp,"primary_end_ms":pe,"late_end_ms":le,"is_baseline":isbase,**tar})
                for rr in srows: series_rows.append({"scheme_id":sid,"split_ms":sp,"primary_end_ms":pe,"late_end_ms":le,**rr})
    schemes_df=pd.DataFrame(schemes); targets_df=pd.DataFrame(target_rows); series_df=pd.DataFrame(series_rows)
    schemes_df.to_csv(out/"03_window_schemes.csv",index=False); targets_df.to_csv(out/"04_experimental_targets_by_window.csv",index=False); series_df.to_csv(out/"04A_experimental_series_metrics_by_window.csv.gz",index=False,compression="gzip")

    prefix=model_prefix_all(full,out,args.model_dt_ms,args.chunk_size,args.resume)
    if not prefix.integration_ok.all(): raise RuntimeError(f"Model integration validity failed for {(~prefix.integration_ok).sum()} candidates")
    full2=full[[
        "prior","sample_id","calibration_score_corrected","accepted_corrected",
        "control_plateau_PO_140_250","control_early_PO_140_200",
        "control_late_PO_200_500","control_integrated_open_state_140_500_ms"
    ]].merge(prefix,on=["prior","sample_id"],validate="one_to_one")

    # The prefix integrator is used only to estimate within-candidate changes when
    # the operational window boundaries move. Absolute baseline observables are
    # anchored to the exact frozen Step-24 values. This paired-ratio construction
    # removes small solver/discretization offsets while preserving each model's
    # relative response to changing the analysis window.
    h0a,e0a,l0a,j0a=model_metrics(full2,200.0,250.0,500.0)
    h0f=full2["control_plateau_PO_140_250"].to_numpy(float)
    e0f=full2["control_early_PO_140_200"].to_numpy(float)
    l0f=full2["control_late_PO_200_500"].to_numpy(float)
    j0f=full2["control_integrated_open_state_140_500_ms"].to_numpy(float)
    valid_anchor=(np.isfinite(h0a)&np.isfinite(e0a)&np.isfinite(l0a)&np.isfinite(j0a)&
                  np.isfinite(h0f)&np.isfinite(e0f)&np.isfinite(l0f)&np.isfinite(j0f)&
                  (h0a>0)&(e0a>0)&(l0a>0)&(j0a>0)&(h0f>0)&(e0f>0)&(l0f>0)&(j0f>0))
    if not np.all(valid_anchor):
        raise RuntimeError(f"Baseline anchoring invalid for {(~valid_anchor).sum()} candidates")
    anchor_log_errors=np.column_stack([
        np.abs(np.log(h0a/h0f)),np.abs(np.log(e0a/e0f)),
        np.abs(np.log(l0a/l0f)),np.abs(np.log(j0a/j0f))
    ])
    audit["prefix_vs_frozen_baseline_observable_max_abs_log_error"]=float(np.max(anchor_log_errors))
    audit["prefix_vs_frozen_baseline_observable_median_abs_log_error"]=float(np.median(anchor_log_errors))
    # A gross mismatch would indicate the wrong model/forcing implementation.
    if np.max(anchor_log_errors)>0.05:
        write_json(out/"00_input_audit.json",audit)
        raise RuntimeError(
            f"Prefix integrator differs too much from frozen baseline observables: "
            f"max abs log error {np.max(anchor_log_errors):g}"
        )
    topo=pd.read_csv(step29/"02_final_strong40_topology.csv"); ids40=set(topo.sample_id.astype(int)); ids36=set(topo.loc[topo.strict_component_size==36,"sample_id"].astype(int)); ids35=set(topo.loc[topo.native_component_size==35,"sample_id"].astype(int))
    reps=pd.read_csv(step29/"05_publication_representatives.csv"); primary_id=int(reps.loc[reps.role=="primary_typical_dominant_native","sample_id"].iloc[0])

    summary=[]; retention=[]; baseline_sets={}
    score_matrix={}
    for _,sr in schemes_df.iterrows():
        sid=sr.scheme_id; raw_tar=targets_df[targets_df.scheme_id==sid].iloc[0].to_dict()
        tar=frozen_target if bool(sr.is_baseline) else raw_tar
        ha,ea,la,ja=model_metrics(full2,float(sr.split_ms),float(sr.primary_end_ms),float(sr.late_end_ms))
        # Paired-ratio anchoring: baseline is exactly the frozen Step-24 observable;
        # alternate windows inherit only the relative change measured by the same
        # prefix integrator. At the baseline scheme all four ratios equal 1 exactly.
        h=h0f*(ha/h0a); e=e0f*(ea/e0a); l=l0f*(la/l0a); j=j0f*(ja/j0a)
        s=score_from_metrics(h,e,l,j,tar); score_matrix[sid]=s
        if bool(sr.is_baseline):
            err=np.max(np.abs(s-full2.calibration_score_corrected.to_numpy(float)))
            audit["baseline_model_score_max_abs_error"]=float(err)
            if err>2e-8:
                write_json(out/"00_input_audit.json",audit); raise RuntimeError(f"Anchored baseline model score replay failed: max abs error {err:g}")
        for prior in ["broad","reference"]:
            mask=full2.prior.astype(str).str.lower().to_numpy()==prior
            sub=full2.loc[mask,["sample_id"]].copy(); ss=s[mask]; order=np.argsort(ss,kind="mergesort"); take=order[:5000]
            aset=set(sub.iloc[take].sample_id.astype(int)); cutoff=float(ss[take[-1]])
            if bool(sr.is_baseline):
                baseline_sets[prior]=aset
                frozen_set=set(full2.loc[(full2.prior.astype(str).str.lower()==prior)&(full2.accepted_corrected.astype(bool)),"sample_id"].astype(int))
                if aset!=frozen_set:
                    raise RuntimeError(f"Baseline accepted-set replay failed for {prior}: symmetric difference {len(aset^frozen_set)}")
            summary.append({"scheme_id":sid,"prior":prior,"split_ms":sr.split_ms,"primary_end_ms":sr.primary_end_ms,"late_end_ms":sr.late_end_ms,"is_baseline":bool(sr.is_baseline),"accepted_n":5000,"cutoff_score":cutoff,"final40_retained_n":len(aset&ids40) if prior=="broad" else np.nan,"dominant36_retained_n":len(aset&ids36) if prior=="broad" else np.nan,"native35_retained_n":len(aset&ids35) if prior=="broad" else np.nan,"primary_representative_retained":primary_id in aset if prior=="broad" else np.nan})
            if prior=="broad":
                for sid40 in sorted(ids40): retention.append({"scheme_id":sid,"sample_id":sid40,"retained":sid40 in aset,"final_regime":str(topo.set_index("sample_id").loc[sid40,"final_regime"])})
    summ=pd.DataFrame(summary)
    for prior in ["broad","reference"]:
        bset=baseline_sets[prior]
        for i,row in summ[summ.prior==prior].iterrows():
            sid=row.scheme_id; sr=schemes_df[schemes_df.scheme_id==sid].iloc[0]; s=score_matrix[sid]; mask=full2.prior.astype(str).str.lower().to_numpy()==prior; sub=full2.loc[mask,["sample_id"]]; order=np.argsort(s[mask],kind="mergesort")[:5000]; aset=set(sub.iloc[order].sample_id.astype(int)); summ.loc[i,"overlap_with_baseline_n"]=len(aset&bset); summ.loc[i,"jaccard_vs_baseline"]=jaccard(aset,bset)
    summ.to_csv(out/"05_model_window_sensitivity_by_prior.csv",index=False)
    pd.DataFrame(retention).to_csv(out/"06_final40_window_retention.csv",index=False)
    plot_window_sensitivity(summ,targets_df,out)

    broad=summ[summ.prior=="broad"]
    base_target=frozen_target
    decision={
        "status":"WINDOW_ROBUSTNESS_AND_GEOMETRY_COMPLETED",
        "n_window_schemes":int(len(schemes_df)),
        "baseline_scheme":scheme_id(*BASELINE_SCHEME),
        "broad_jaccard_min":float(broad.jaccard_vs_baseline.min()),
        "broad_jaccard_median":float(broad.jaccard_vs_baseline.median()),
        "final40_retained_min":int(broad.final40_retained_n.min()),
        "final40_retained_median":float(broad.final40_retained_n.median()),
        "dominant36_retained_min":int(broad.dominant36_retained_n.min()),
        "primary_representative_retained_all_schemes":bool(broad.primary_representative_retained.astype(bool).all()),
        "experimental_target_relative_ranges":{
            "E_over_H":[float(targets_df.target_E_over_H.min()/base_target["target_E_over_H"]),float(targets_df.target_E_over_H.max()/base_target["target_E_over_H"])],
            "L_over_H":[float(targets_df.target_L_over_H.min()/base_target["target_L_over_H"]),float(targets_df.target_L_over_H.max()/base_target["target_L_over_H"])],
            "Q_over_H":[float(targets_df.target_Q_over_H_ms.min()/base_target["target_Q_over_H_ms"]),float(targets_df.target_Q_over_H_ms.max()/base_target["target_Q_over_H_ms"])],
        },
        "interpretation_boundary":"Descriptive sensitivity analysis of prespecified operational windows. No new acceptance rule, response threshold, denominator cutoff, or Bayesian posterior is introduced.",
    }
    write_json(out/"10_scientific_summary.json",decision); audit["scientific_summary"]=decision; write_json(out/"00_input_audit.json",audit)

    readme=f"""# Step 31 results — final window robustness and parameter-space geometry\n\nStatus: **{decision['status']}**\n\nThe post-train start was held fixed at 140 ms. The operational boundaries were varied on a prespecified 3x3x3 grid: split = 180/200/220 ms, primary-window end = 230/250/270 ms, late/integral end = 450/500/550 ms. For every scheme, untreated SlowEPSC reference recordings were re-extracted from raw ABF data. The exact five frozen series-level baseline ratios were retained as anchors, while raw ABFs supplied only the within-series changes caused by moving the window boundaries. Experimental targets and floor-bounded robust tolerances were then recomputed, and all 100,000 Step-10 candidates were recalibrated. Model-window sensitivity uses a paired-ratio anchor: exact frozen Step-24 baseline observables are retained, while the prefix integrator supplies only the within-candidate change caused by moving a window boundary.\n\nThe baseline scheme 140-200 / 140-250 / 200-500 / 140-500 was required to replay all five frozen series-level ratios after anchoring and the frozen candidate calibration scores before any sensitivity result was accepted. Ro25 and memantine condition files are hard-excluded from calibration mapping.\n\nHeadline sensitivity quantities:\n\n- minimum broad-cohort Jaccard versus baseline: **{decision['broad_jaccard_min']:.3f}**\n- median broad-cohort Jaccard: **{decision['broad_jaccard_median']:.3f}**\n- minimum retained count among the frozen 40 experiment-scale solutions: **{decision['final40_retained_min']}/40**\n- minimum retained count from the dominant 36-solution network: **{decision['dominant36_retained_min']}/36**\n- primary typical representative retained under every scheme: **{decision['primary_representative_retained_all_schemes']}**\n\n`Fig31B_parameter_space_geometry` is the manuscript-facing geometry figure. Panel A shows all 5,000 control-compatible broad-prior models in a two-dimensional PCA projection of standardized log-rates, with the 40 experiment-scale solutions overlaid. Panel B adds only verified connectivity edges. Panel C shows calibration score against amplification magnitude for the 40 solutions. Panel D explicitly defines what a verified path means by plotting calibration score and amplification ratio along one path; this is a path through parameter space, not a receptor-state trajectory.\n\n`Fig31C_rate_filtering_marginals` shows sampled broad, control-compatible, and experiment-scale marginal rate distributions. These are filtering distributions from the search/calibration procedure and are **not** interpreted as Bayesian priors or posteriors.\n"""
    (out/"README_RESULTS.md").write_text(readme,encoding="utf-8")
    print(json.dumps(decision,indent=2),flush=True)
    print(f"Step31 complete: {out}",flush=True)


if __name__ == "__main__":
    main()
