#!/usr/bin/env python3
"""Reproduce the compact publication summaries from distributed source tables.

This script does not require the raw ABF recordings or the upstream Step10
parameter-sampling stage.  It verifies the frozen endpoint and topology-first
claims and regenerates simple publication-oriented diagnostic figures.
"""
from __future__ import annotations

from pathlib import Path
import argparse
import math
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def close(a, b, tol=5e-6):
    return math.isclose(float(a), float(b), rel_tol=0.0, abs_tol=tol)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[1])
    ap.add_argument("--output", type=Path, default=None)
    args = ap.parse_args()
    root = args.repo.resolve()
    out = (args.output or root / "reproduced").resolve()
    out.mkdir(parents=True, exist_ok=True)

    primary = pd.read_csv(root / "data/derived/plateau_minus40_primary.csv")
    targets = pd.read_csv(root / "data/derived/potentiating_targets_minus40.csv")
    topo = pd.read_csv(root / "data/model_outputs/base_topology_map.csv")
    amp = pd.read_csv(root / "data/model_outputs/bslow_paired_amplification.csv")
    summary = pd.read_csv(root / "provenance/publication_summary_table.csv")

    counts = primary.plateau_class_frozen.value_counts().to_dict()
    assert counts.get("plateau_potentiating", 0) == 5
    assert counts.get("plateau_suppressing", 0) == 5
    assert counts.get("plateau_neutral", 0) == 1
    assert close(targets.r_plateau.min(), 1.6359295439505144)
    assert close(targets.r_plateau.max(), 4.072558998570809)

    for prior, any_fraction, node_fraction in [
        ("broad", 0.3124, 1.0),
        ("reference", 0.9850, 0.9250),
    ]:
        s = summary[summary.prior == prior].iloc[0]
        assert close(s.Base_fraction_parameter_sets_any_Braess_somewhere, any_fraction)
        assert close(s.Base_fraction_forcing_nodes_with_any_Braess, node_fraction)

    report = pd.DataFrame([
        {"quantity":"potentiating_cells_minus40", "value":5},
        {"quantity":"suppressing_cells_minus40", "value":5},
        {"quantity":"neutral_cells_minus40", "value":1},
        {"quantity":"experimental_target_min", "value":targets.r_plateau.min()},
        {"quantity":"experimental_target_max", "value":targets.r_plateau.max()},
        {"quantity":"Base_broad_parameter_sets_any_Braess", "value":0.3124},
        {"quantity":"Base_reference_parameter_sets_any_Braess", "value":0.9850},
    ])
    report.to_csv(out / "verified_key_results.csv", index=False)

    # Experimental ratios.
    g = primary.sort_values("r_plateau")
    fig, ax = plt.subplots(figsize=(8, 4.8))
    ax.scatter(np.arange(len(g)), g.r_plateau)
    ax.axhline(1.0, linestyle="--", linewidth=1)
    ax.set_xticks(np.arange(len(g)))
    ax.set_xticklabels(g.cell_id, rotation=60, ha="right", fontsize=8)
    ax.set_ylabel("r_plateau")
    fig.tight_layout()
    fig.savefig(out / "experimental_plateau_ratios.png", dpi=220)
    plt.close(fig)

    # Base topology maps.
    for prior in ("broad", "reference"):
        x = topo[topo.prior == prior]
        p = x.pivot(index="tau_G_ms", columns="G_peak", values="braess_fraction").sort_index().sort_index(axis=1)
        fig, ax = plt.subplots(figsize=(7, 5))
        im = ax.imshow(p.values, origin="lower", aspect="auto")
        ax.set_xticks(np.arange(len(p.columns))); ax.set_xticklabels([f"{v:g}" for v in p.columns], rotation=45)
        ax.set_yticks(np.arange(len(p.index))); ax.set_yticklabels([f"{v:g}" for v in p.index])
        ax.set_xlabel("G_peak"); ax.set_ylabel("tau_G (ms)")
        fig.colorbar(im, ax=ax, label="fraction with r_plateau > 1")
        fig.tight_layout()
        fig.savefig(out / f"base_braess_map_{prior}.png", dpi=220)
        plt.close(fig)

    # Paired B-slow response.
    for prior in ("broad", "reference"):
        x = amp[amp.prior == prior].sort_values("b_timescale_scale")
        fig, ax = plt.subplots(figsize=(7, 4.5))
        ax.plot(x.b_timescale_scale, x.slowed_grid_context_braess_fraction_same_300, marker="o")
        ax.axhline(float(x.base_grid_context_braess_fraction_same_300.iloc[0]), linestyle="--", linewidth=1)
        ax.set_xscale("log")
        ax.set_xlabel("B-route timescale scale sB")
        ax.set_ylabel("fraction of matched contexts with r_plateau > 1")
        fig.tight_layout()
        fig.savefig(out / f"bslow_amplification_{prior}.png", dpi=220)
        plt.close(fig)

    print(f"PASS: publication tables verified; outputs written to {out}")


if __name__ == "__main__":
    main()
