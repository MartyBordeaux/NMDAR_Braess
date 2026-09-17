#!/usr/bin/env python3
from __future__ import annotations

import argparse
import io
import json
import math
import os
import zipfile
from pathlib import Path
from typing import Iterable

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

VERSION = '1.0.0'
STRONG_DEFAULT = 1.6359295439505144
ACCEPTED_DEFAULT = 5000
WITNESS_DEFAULT = 37318
THRESHOLDS_DEFAULT = [1.5, 1.6, STRONG_DEFAULT, 1.7, 1.8]
HERE = Path(__file__).resolve().parent


def log(msg: str) -> None:
    print(msg, flush=True)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description='Correct the Step10 model charge window from 0-600 ms to 140-500 ms and audit downstream consequences.')
    p.add_argument('--step10', default=None, help='CSV path or zip::member for matched_parameter_all_samples.csv')
    p.add_argument('--step19', default=None, help='Step19 results directory/zip; used as broad-only fallback')
    p.add_argument('--step18', default=None, help='Step18 results directory/zip')
    p.add_argument('--step20', default=None, help='Step20 results directory/zip')
    p.add_argument('--step21', default=None, help='Step21 results directory/zip')
    p.add_argument('--step23', default=None, help='Step23 results directory/zip')
    p.add_argument('--targets', default=str(HERE / 'inputs' / 'calibration_targets_step4.json'))
    p.add_argument('--out', default='results_step_24_calibration_window_correction')
    p.add_argument('--accepted-per-prior', type=int, default=ACCEPTED_DEFAULT)
    p.add_argument('--strong-threshold', type=float, default=STRONG_DEFAULT)
    p.add_argument('--witness-id', type=int, default=WITNESS_DEFAULT)
    p.add_argument('--allow-broad-only', action='store_true', help='Allow final completion without the full two-prior Step10 table.')
    return p.parse_args()


def candidate_roots() -> list[Path]:
    roots = []
    for x in [os.getcwd(), str(Path(os.getcwd()).parent), '/root/nmda2', '/root/nmda', '/home/vlad/nmda2', '/home/vlad/nmda']:
        p = Path(x)
        if p.exists():
            p = p.resolve()
            if p not in roots:
                roots.append(p)
    return roots


def read_csv_spec(spec: str, **kwargs) -> pd.DataFrame:
    if '::' in str(spec):
        z, member = str(spec).split('::', 1)
        with zipfile.ZipFile(z) as Z:
            b = Z.read(member)
        if member.endswith('.gz'):
            import gzip
            b = gzip.decompress(b)
        return pd.read_csv(io.BytesIO(b), **kwargs)
    return pd.read_csv(spec, compression='infer', **kwargs)


def read_json_spec(spec: str) -> dict:
    if '::' in str(spec):
        z, member = str(spec).split('::', 1)
        with zipfile.ZipFile(z) as Z:
            return json.loads(Z.read(member).decode('utf-8'))
    return json.loads(Path(spec).read_text())


def discover_step10(explicit: str | None, roots: Iterable[Path]) -> str | None:
    if explicit:
        return explicit
    preferred = [
        '/root/nmda/results_step10.zip::results/10_matched_rerouting/matched_parameter_all_samples.csv',
        '/root/nmda/results/10_matched_rerouting/matched_parameter_all_samples.csv',
        '/root/nmda2/results_step10.zip::results/10_matched_rerouting/matched_parameter_all_samples.csv',
    ]
    for s in preferred:
        p = Path(s.split('::', 1)[0])
        if p.exists():
            return s
    for root in roots:
        for z in root.rglob('*step10*.zip'):
            try:
                with zipfile.ZipFile(z) as Z:
                    for m in Z.namelist():
                        if m.endswith('matched_parameter_all_samples.csv'):
                            return f'{z}::{m}'
            except Exception:
                pass
        for f in root.rglob('matched_parameter_all_samples.csv'):
            return str(f)
    return None


def discover_results(explicit: str | None, step_num: int, stem: str, roots: Iterable[Path]) -> str | None:
    if explicit:
        return explicit
    preferred = [
        Path(f'/root/nmda2/step_{step_num}/results_step_{step_num}_{stem}'),
        Path(f'/root/nmda2/step_{step_num}/results_step_{step_num}_{stem}.zip'),
        Path(f'/root/nmda2/results_step_{step_num}_{stem}'),
        Path(f'/root/nmda2/results_step_{step_num}_{stem}.zip'),
    ]
    for p in preferred:
        if p.exists():
            return str(p)
    names = [f'results_step_{step_num}_{stem}', f'results_step_{step_num}_{stem}.zip']
    for root in roots:
        for name in names:
            hits = list(root.rglob(name))
            if hits:
                return str(hits[0])
    return None


def result_member(spec: str, filename: str) -> str:
    p = Path(spec)
    if p.is_dir():
        q = p / filename
        if not q.exists():
            raise FileNotFoundError(q)
        return str(q)
    if p.suffix == '.zip':
        with zipfile.ZipFile(p) as Z:
            matches = [m for m in Z.namelist() if m == filename or m.endswith('/' + filename)]
        if not matches:
            raise FileNotFoundError(f'{filename} not found in {p}')
        return f'{p}::{matches[0]}'
    raise ValueError(f'Unsupported results spec: {spec}')


def read_result_csv(spec: str, filename: str) -> pd.DataFrame:
    return read_csv_spec(result_member(spec, filename))


def read_result_json(spec: str, filename: str) -> dict:
    return read_json_spec(result_member(spec, filename))


def discover_step19_fallback(explicit: str | None, roots: Iterable[Path]) -> str | None:
    return discover_results(explicit, 19, 'calibration_frontier', roots)


def normalize_step10(df: pd.DataFrame) -> pd.DataFrame:
    d = df.copy()
    if 'accepted_control_calibration' not in d.columns and 'accepted' in d.columns:
        d['accepted_control_calibration'] = d['accepted'].astype(bool)
    if 'prior' not in d.columns:
        d['prior'] = 'broad'
    needed = [
        'prior', 'sample_id', 'calibration_score', 'accepted_control_calibration',
        'control_plateau_PO_140_250', 'control_early_PO_140_200',
        'control_late_PO_200_500', 'control_charge_PO_ms'
    ]
    missing = [c for c in needed if c not in d.columns]
    if missing:
        raise RuntimeError(f'Step10/fallback table missing columns: {missing}')
    for c in ['sample_id']:
        d[c] = pd.to_numeric(d[c], errors='raise').astype(int)
    for c in ['calibration_score', 'control_plateau_PO_140_250', 'control_early_PO_140_200', 'control_late_PO_200_500', 'control_charge_PO_ms']:
        d[c] = pd.to_numeric(d[c], errors='coerce')
    d['accepted_control_calibration'] = d['accepted_control_calibration'].astype(bool)
    if 'control_valid' not in d.columns:
        d['control_valid'] = 1.0
    d['control_valid'] = pd.to_numeric(d['control_valid'], errors='coerce').fillna(0.0)
    d['prior'] = d['prior'].astype(str).str.lower()
    return d


def load_step10_or_fallback(step10: str | None, step19: str | None, roots: list[Path]) -> tuple[pd.DataFrame, str, bool]:
    if step10:
        return normalize_step10(read_csv_spec(step10)), step10, False
    discovered = discover_step10(None, roots)
    if discovered:
        return normalize_step10(read_csv_spec(discovered)), discovered, False
    step19 = discover_step19_fallback(step19, roots)
    if not step19:
        raise RuntimeError('Could not discover Step10 all-samples table or Step19 fallback.')
    d = read_result_csv(step19, '00_merged_frontier_input.csv.gz')
    d['prior'] = 'broad'
    if 'accepted_control_calibration' not in d.columns:
        d['accepted_control_calibration'] = d['accepted']
    return normalize_step10(d), f'{step19}::00_merged_frontier_input.csv.gz', True


def score_from_observables(H: np.ndarray, E: np.ndarray, L: np.ndarray, J: np.ndarray, valid: np.ndarray, cal: dict) -> np.ndarray:
    H = np.asarray(H, float); E = np.asarray(E, float); L = np.asarray(L, float); J = np.asarray(J, float)
    valid = np.asarray(valid, bool)
    out = np.full(len(H), np.inf, float)
    ok = valid & np.isfinite(H) & np.isfinite(E) & np.isfinite(L) & np.isfinite(J) & (H > 1e-12) & (E > 0) & (L > 0) & (J > 0)
    if not np.any(ok):
        return out
    re = E[ok] / H[ok]
    rl = L[ok] / H[ok]
    rj = J[ok] / H[ok]
    out[ok] = (
        (np.log(re / float(cal['early_over_primary_median'])) / float(cal.get('early_over_primary_log_tolerance', 0.35))) ** 2
        + (np.log(rl / float(cal['late_over_primary_median'])) / float(cal.get('late_over_primary_log_tolerance', 0.35))) ** 2
        + 0.5 * (np.log(rj / float(cal['charge_over_primary_ms_median'])) / float(cal.get('charge_over_primary_log_tolerance', 0.50))) ** 2
    )
    return out


def midpoint_cutoff(scores: np.ndarray, accepted_mask: np.ndarray) -> tuple[float, float, float]:
    a = np.asarray(scores)[accepted_mask]
    r = np.asarray(scores)[~accepted_mask & np.isfinite(scores)]
    max_a = float(np.max(a)) if len(a) else float('nan')
    min_r = float(np.min(r)) if len(r) else float('nan')
    mid = (max_a + min_r) / 2 if np.isfinite(max_a) and np.isfinite(min_r) else max_a
    return max_a, min_r, mid


def assign_corrected(df: pd.DataFrame, cal: dict, accepted_per_prior: int) -> pd.DataFrame:
    d = df.copy()
    H = d['control_plateau_PO_140_250'].to_numpy(float)
    E = d['control_early_PO_140_200'].to_numpy(float)
    L = d['control_late_PO_200_500'].to_numpy(float)
    Jold = d['control_charge_PO_ms'].to_numpy(float)
    valid = d['control_valid'].to_numpy(float) > 0.5
    Jnew = 60.0 * E + 300.0 * L
    d['control_integrated_open_state_140_500_ms'] = Jnew
    d['charge_window_fraction_of_historical_full'] = np.divide(Jnew, Jold, out=np.full(len(d), np.nan), where=np.isfinite(Jold) & (np.abs(Jold) > 1e-12))
    d['calibration_score_replayed_historical'] = score_from_observables(H, E, L, Jold, valid, cal)
    d['calibration_score_corrected'] = score_from_observables(H, E, L, Jnew, valid, cal)
    d['accepted_corrected'] = False
    d['corrected_rank'] = np.nan
    d['historical_rank_recomputed'] = np.nan
    for prior, idx in d.groupby('prior').groups.items():
        ii = np.asarray(list(idx), int)
        finite = ii[np.isfinite(d.loc[ii, 'calibration_score_corrected'].to_numpy(float))]
        if len(finite) < accepted_per_prior:
            raise RuntimeError(f'{prior}: only {len(finite)} finite corrected scores for {accepted_per_prior} retained.')
        ord_new = finite[np.argsort(d.loc[finite, 'calibration_score_corrected'].to_numpy(float), kind='mergesort')]
        chosen = ord_new[:accepted_per_prior]
        d.loc[chosen, 'accepted_corrected'] = True
        rank_new = pd.Series(d.loc[ii, 'calibration_score_corrected'].to_numpy(float)).rank(method='first').to_numpy()
        rank_old = pd.Series(d.loc[ii, 'calibration_score'].to_numpy(float)).rank(method='first').to_numpy()
        d.loc[ii, 'corrected_rank'] = rank_new
        d.loc[ii, 'historical_rank_recomputed'] = rank_old
    return d


def acceptance_summary(d: pd.DataFrame, accepted_per_prior: int) -> pd.DataFrame:
    rows = []
    for prior, g in d.groupby('prior'):
        old = g['accepted_control_calibration'].to_numpy(bool)
        new = g['accepted_corrected'].to_numpy(bool)
        inter = int(np.sum(old & new)); union = int(np.sum(old | new)); sym = int(np.sum(old ^ new))
        old_max, old_min, old_mid = midpoint_cutoff(g['calibration_score'].to_numpy(float), old)
        new_max, new_min, new_mid = midpoint_cutoff(g['calibration_score_corrected'].to_numpy(float), new)
        corr = float(pd.Series(g['calibration_score']).corr(pd.Series(g['calibration_score_corrected']), method='spearman'))
        rows.append({
            'prior': prior, 'n_candidates': len(g), 'accepted_target_n': accepted_per_prior,
            'historical_accepted_n': int(old.sum()), 'corrected_accepted_n': int(new.sum()),
            'intersection_n': inter, 'symmetric_difference_n': sym,
            'jaccard': inter / union if union else float('nan'), 'score_spearman': corr,
            'historical_max_accepted_score': old_max, 'historical_min_rejected_score': old_min, 'historical_cutoff_midpoint': old_mid,
            'corrected_max_accepted_score': new_max, 'corrected_min_rejected_score': new_min, 'corrected_cutoff_midpoint': new_mid,
            'median_corrected_window_fraction_of_full_charge': float(np.nanmedian(g['charge_window_fraction_of_historical_full'])),
        })
    return pd.DataFrame(rows)


def reconstruct_threshold_strong(step18: str, step23: str | None, threshold: float) -> pd.DataFrame:
    grid = read_result_csv(step18, '03_candidate_grid_summary.csv').copy()
    grid['sample_id'] = pd.to_numeric(grid['sample_id']).astype(int)
    if abs(threshold - STRONG_DEFAULT) < 1e-10 and 'strong_anywhere_confirmed' in grid.columns:
        grid['strong_final'] = grid['strong_anywhere_confirmed'].astype(bool)
        return grid
    grid['strong_final'] = pd.to_numeric(grid['r_max_screen_dt0p2'], errors='coerce') >= threshold
    if step23:
        try:
            guard = read_result_csv(step23, '01A_threshold_exact_confirmations.csv.gz')
            gg = guard[np.isclose(pd.to_numeric(guard['threshold'], errors='coerce'), threshold)].copy()
            if len(gg):
                repl = dict(zip(pd.to_numeric(gg['sample_id']).astype(int), gg['strong_final'].astype(bool)))
                m = grid['sample_id'].isin(repl)
                grid.loc[m, 'strong_final'] = grid.loc[m, 'sample_id'].map(repl).astype(bool)
        except Exception as e:
            log(f'WARNING: Step23 exact-threshold guard unavailable for {threshold}: {e}')
    return grid


def corrected_score_for_point_table(df: pd.DataFrame, cal: dict) -> np.ndarray:
    H = pd.to_numeric(df['control_plateau_PO_140_250'], errors='coerce').to_numpy(float)
    E = pd.to_numeric(df['control_early_PO_140_200'], errors='coerce').to_numpy(float)
    L = pd.to_numeric(df['control_late_PO_200_500'], errors='coerce').to_numpy(float)
    J = 60.0 * E + 300.0 * L
    valid = np.ones(len(df), bool)
    return score_from_observables(H, E, L, J, valid, cal)


def plot_score_shift(broad: pd.DataFrame, old_cut: float, new_cut: float, outdir: Path) -> None:
    x = broad['calibration_score'].to_numpy(float)
    y = broad['calibration_score_corrected'].to_numpy(float)
    ok = np.isfinite(x) & np.isfinite(y)
    fig, ax = plt.subplots(figsize=(6.4, 5.2))
    ax.scatter(x[ok], y[ok], s=4, alpha=0.15)
    ax.axvline(old_cut, linewidth=1)
    ax.axhline(new_cut, linewidth=1)
    ax.set_xlabel('Historical calibration score')
    ax.set_ylabel('Corrected 140-500 ms calibration score')
    ax.set_title('Calibration-score shift after window correction')
    fig.tight_layout()
    fig.savefig(outdir / 'Fig24A_score_shift.pdf')
    fig.savefig(outdir / 'Fig24A_score_shift.png', dpi=220)
    plt.close(fig)


def plot_strong_shift(tab: pd.DataFrame, outdir: Path) -> None:
    fig, ax = plt.subplots(figsize=(6.4, 4.5))
    ax.plot(tab['threshold'], tab['historical_accepted_strong_n'], marker='o', label='historical retained')
    ax.plot(tab['threshold'], tab['corrected_accepted_strong_n'], marker='o', label='corrected retained')
    ax.set_xlabel('Strong-response threshold')
    ax.set_ylabel('Strong candidates among 5,000 retained broad sets')
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(outdir / 'Fig24B_strong_support_shift.pdf')
    fig.savefig(outdir / 'Fig24B_strong_support_shift.png', dpi=220)
    plt.close(fig)


def plot_local_shift(step20: pd.DataFrame | None, step21: pd.DataFrame | None, outdir: Path) -> None:
    if step20 is None and step21 is None:
        return
    fig, ax = plt.subplots(figsize=(6.7, 4.7))
    if step20 is not None and len(step20):
        ax.plot(step20['perturbation_pct'], step20['historical_strong_and_compatible_fraction'], marker='o', linestyle='--', label='Step20 historical')
        ax.plot(step20['perturbation_pct'], step20['corrected_strong_and_compatible_fraction'], marker='o', label='Step20 corrected')
    if step21 is not None and len(step21):
        # Map factor to a separate monotonic pseudo-percent axis only if Step20 absent; otherwise avoid mixing units.
        pass
    ax.set_xlabel('Simultaneous rate perturbation (%)')
    ax.set_ylabel('Strong + control-compatible fraction')
    ax.set_ylim(-0.02, 1.02)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(outdir / 'Fig24C_local_cloud_shift.pdf')
    fig.savefig(outdir / 'Fig24C_local_cloud_shift.png', dpi=220)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    outdir = Path(args.out).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    roots = candidate_roots()

    cal_payload = read_json_spec(args.targets)
    cal = cal_payload['calibration'] if 'calibration' in cal_payload else cal_payload

    step10_discovered = discover_step10(args.step10, roots)
    step19 = discover_step19_fallback(args.step19, roots)
    d, step10_source, broad_only = load_step10_or_fallback(step10_discovered, step19, roots)
    step18 = discover_results(args.step18, 18, 'prior_geometry', roots)
    step20 = discover_results(args.step20, 20, 'exact_witness_compatibility', roots)
    step21 = discover_results(args.step21, 21, 'basin_population', roots)
    step23 = discover_results(args.step23, 23, 'publication_consolidation', roots)

    audit = {
        'version': VERSION, 'step10_source': step10_source, 'broad_only_fallback': broad_only,
        'step18': step18, 'step20': step20, 'step21': step21, 'step23': step23,
        'targets': str(Path(args.targets).resolve()), 'accepted_per_prior': args.accepted_per_prior,
        'strong_threshold': args.strong_threshold, 'witness_id': args.witness_id,
        'priors_found': d['prior'].value_counts().to_dict(),
        'correction': 'model integrated occupancy window 0-600 ms -> 140-500 ms; all experimental targets and other score terms frozen',
        'corrected_integral_formula': 'J_140_500_ms = 60*E_140_200 + 300*L_200_500',
    }
    (outdir / '00_input_audit.json').write_text(json.dumps(audit, indent=2) + '\n')

    log('[1/7] Replaying historical score and computing corrected 140-500 ms score')
    d = assign_corrected(d, cal, args.accepted_per_prior)
    frozen = d['calibration_score'].to_numpy(float)
    replay = d['calibration_score_replayed_historical'].to_numpy(float)
    finite = np.isfinite(frozen) & np.isfinite(replay)
    maxerr = float(np.max(np.abs(frozen[finite] - replay[finite]))) if np.any(finite) else float('inf')
    n_nonfinite_mismatch = int(np.sum(np.isfinite(frozen) != np.isfinite(replay)))
    replay_status = 'PASS' if maxerr < 1e-8 and n_nonfinite_mismatch == 0 else 'FAIL'
    replay_json = {
        'status': replay_status, 'max_abs_score_error': maxerr, 'nonfinite_pattern_mismatches': n_nonfinite_mismatch,
        'criterion': 'max abs score error < 1e-8 and identical finite/nonfinite pattern',
        'target_values': cal,
    }
    (outdir / '01_historical_score_replay.json').write_text(json.dumps(replay_json, indent=2) + '\n')
    if replay_status != 'PASS':
        raise RuntimeError('Historical calibration-score replay failed; refusing correction analysis.')

    keep_cols = [
        'prior','sample_id','calibration_score','calibration_score_replayed_historical','calibration_score_corrected',
        'historical_rank_recomputed','corrected_rank','accepted_control_calibration','accepted_corrected',
        'control_plateau_PO_140_250','control_early_PO_140_200','control_late_PO_200_500','control_charge_PO_ms',
        'control_integrated_open_state_140_500_ms','charge_window_fraction_of_historical_full'
    ]
    extra = [c for c in ['pulse_amplitude_au','glutamate_tau_ms'] if c in d.columns]
    d[keep_cols + extra].to_csv(outdir / '02_corrected_candidate_scores.csv.gz', index=False, compression='gzip')

    log('[2/7] Comparing historical and corrected retained sets')
    acc = acceptance_summary(d, args.accepted_per_prior)
    acc.to_csv(outdir / '03_acceptance_shift_by_prior.csv', index=False)
    broad = d[d['prior'].eq('broad')].copy()
    if len(broad) == 0:
        raise RuntimeError('No broad prior found; downstream manuscript audit requires broad candidates.')
    broad_acc_row = acc[acc['prior'].eq('broad')].iloc[0]
    old_cut = float(broad_acc_row['historical_cutoff_midpoint'])
    new_cut = float(broad_acc_row['corrected_cutoff_midpoint'])

    witness = broad[broad['sample_id'].eq(args.witness_id)]
    witness_info = None
    if len(witness) == 1:
        w = witness.iloc[0]
        witness_info = {
            'sample_id': args.witness_id, 'historical_score': float(w['calibration_score']), 'corrected_score': float(w['calibration_score_corrected']),
            'historical_accepted': bool(w['accepted_control_calibration']), 'corrected_accepted': bool(w['accepted_corrected']),
            'historical_rank': float(w['historical_rank_recomputed']), 'corrected_rank': float(w['corrected_rank'])
        }
    (outdir / '03A_witness_shift.json').write_text(json.dumps(witness_info, indent=2) + '\n')
    plot_score_shift(broad, old_cut, new_cut, outdir)

    strong_primary = None
    threshold_tab = None
    if step18:
        log('[3/7] Recomputing global strong-response support under corrected retention')
        rows = []
        for th in THRESHOLDS_DEFAULT:
            grid = reconstruct_threshold_strong(step18, step23, th)
            m = broad[['sample_id','accepted_control_calibration','accepted_corrected']].merge(grid[['sample_id','strong_final']], on='sample_id', how='inner', validate='one_to_one')
            st = m['strong_final'].astype(bool).to_numpy(); old = m['accepted_control_calibration'].to_numpy(bool); new = m['accepted_corrected'].to_numpy(bool)
            rows.append({
                'threshold': th, 'precal_strong_n': int(st.sum()), 'precal_strong_fraction': float(st.mean()),
                'historical_accepted_strong_n': int(np.sum(st & old)), 'historical_accepted_strong_fraction': float(np.sum(st & old) / max(1, old.sum())),
                'corrected_accepted_strong_n': int(np.sum(st & new)), 'corrected_accepted_strong_fraction': float(np.sum(st & new) / max(1, new.sum())),
            })
        threshold_tab = pd.DataFrame(rows)
        threshold_tab.to_csv(outdir / '04A_threshold_strong_support_shift.csv', index=False)
        strong_primary = threshold_tab.iloc[np.argmin(np.abs(threshold_tab['threshold'].to_numpy(float) - args.strong_threshold))].to_dict()
        pd.DataFrame([strong_primary]).to_csv(outdir / '04_primary_strong_support_shift.csv', index=False)
        plot_strong_shift(threshold_tab, outdir)

        grid0 = read_result_csv(step18, '03_candidate_grid_summary.csv')
        top = broad[broad['accepted_corrected']][['sample_id','calibration_score_corrected']].merge(grid0[['sample_id','r_max_confirmed_or_screen','G_peak_at_max','tau_G_ms_at_max']], on='sample_id', how='left')
        top = top.sort_values('r_max_confirmed_or_screen', ascending=False).head(100)
        top.to_csv(outdir / '04B_corrected_accepted_top_response_candidates.csv', index=False)

    step20_summary = None
    if step20:
        log('[4/7] Re-scoring Step20 local clouds')
        c = read_result_csv(step20, '04A_cloud_points.csv.gz').copy()
        c['calibration_score_corrected'] = corrected_score_for_point_table(c, cal)
        c['control_compatible_corrected'] = c['calibration_score_corrected'] <= new_cut
        c['strong_and_compatible_corrected'] = c['strong'].astype(bool) & c['control_compatible_corrected'].astype(bool)
        rows = []
        for pct, g in c.groupby('perturbation_pct'):
            rows.append({
                'perturbation_pct': pct, 'n': len(g),
                'historical_compatible_fraction': float(g['control_compatible_midpoint'].astype(bool).mean()),
                'corrected_compatible_fraction': float(g['control_compatible_corrected'].mean()),
                'historical_strong_and_compatible_fraction': float(g['strong_and_compatible'].astype(bool).mean()),
                'corrected_strong_and_compatible_fraction': float(g['strong_and_compatible_corrected'].mean()),
                'delta_joint_fraction': float(g['strong_and_compatible_corrected'].mean() - g['strong_and_compatible'].astype(bool).mean()),
                'median_historical_score': float(np.median(g['calibration_score_exact'])),
                'median_corrected_score': float(np.median(g['calibration_score_corrected'])),
            })
        step20_summary = pd.DataFrame(rows).sort_values('perturbation_pct')
        step20_summary.to_csv(outdir / '05_step20_local_cloud_corrected_summary.csv', index=False)
        c[['perturbation_pct','point_id','strong','calibration_score_exact','calibration_score_corrected','control_compatible_midpoint','control_compatible_corrected','strong_and_compatible','strong_and_compatible_corrected']].to_csv(outdir / '05A_step20_local_cloud_corrected_points.csv.gz', index=False, compression='gzip')

    step21_summary = None
    multistart_summary = None
    if step21:
        log('[5/7] Re-scoring Step21 expanded basin and multistart points')
        b = read_result_csv(step21, '01A_expanded_basin_points.csv.gz').copy()
        b['calibration_score_corrected'] = corrected_score_for_point_table(b, cal)
        b['control_compatible_corrected'] = b['calibration_score_corrected'] <= new_cut
        b['strong_and_compatible_fixed_corrected'] = b['strong_fixed'].astype(bool) & b['control_compatible_corrected'].astype(bool)
        b['strong_and_compatible_highG_corrected'] = b['strong_highG'].astype(bool) & b['control_compatible_corrected'].astype(bool)
        rows = []
        for fac, g in b.groupby('factor'):
            rows.append({
                'factor': fac, 'n': len(g),
                'historical_compatible_fraction': float(g['control_compatible'].astype(bool).mean()),
                'corrected_compatible_fraction': float(g['control_compatible_corrected'].mean()),
                'historical_strong_and_compatible_fixed_fraction': float(g['strong_and_compatible_fixed'].astype(bool).mean()),
                'corrected_strong_and_compatible_fixed_fraction': float(g['strong_and_compatible_fixed_corrected'].mean()),
                'historical_strong_and_compatible_highG_fraction': float(g['strong_and_compatible_highG'].astype(bool).mean()),
                'corrected_strong_and_compatible_highG_fraction': float(g['strong_and_compatible_highG_corrected'].mean()),
                'delta_fixed_joint_fraction': float(g['strong_and_compatible_fixed_corrected'].mean() - g['strong_and_compatible_fixed'].astype(bool).mean()),
                'median_historical_score': float(np.median(g['calibration_score_exact'])),
                'median_corrected_score': float(np.median(g['calibration_score_corrected'])),
            })
        step21_summary = pd.DataFrame(rows).sort_values('factor')
        step21_summary.to_csv(outdir / '06_step21_expanded_basin_corrected_summary.csv', index=False)
        b[['factor','point_id','strong_fixed','strong_highG','calibration_score_exact','calibration_score_corrected','control_compatible','control_compatible_corrected','strong_and_compatible_fixed','strong_and_compatible_fixed_corrected']].to_csv(outdir / '06B_step21_expanded_basin_corrected_points.csv.gz', index=False, compression='gzip')

        try:
            ms = read_result_csv(step21, '03A_multistart_search_points.csv.gz').copy()
            ms['calibration_score_corrected'] = corrected_score_for_point_table(ms, cal)
            ms['control_compatible_corrected'] = ms['calibration_score_corrected'] <= new_cut
            ms['strong_and_compatible_corrected'] = ms['strong'].astype(bool) & ms['control_compatible_corrected'].astype(bool)
            ok = ms['strong_and_compatible_corrected'].astype(bool)
            multistart_summary = {
                'n_points': int(len(ms)), 'historical_strong_and_compatible_n': int(ms['strong_and_compatible'].astype(bool).sum()),
                'corrected_strong_and_compatible_n': int(ok.sum()),
                'corrected_max_r_search': float(ms.loc[ok, 'r_search_max'].max()) if ok.any() else None,
                'corrected_min_score_among_strong': float(ms.loc[ms['strong'].astype(bool), 'calibration_score_corrected'].min()) if ms['strong'].astype(bool).any() else None,
            }
            (outdir / '06A_step21_multistart_corrected_summary.json').write_text(json.dumps(multistart_summary, indent=2) + '\n')
        except Exception as e:
            (outdir / '06A_step21_multistart_corrected_summary.json').write_text(json.dumps({'status':'unavailable','error':repr(e)}, indent=2) + '\n')

    plot_local_shift(step20_summary, step21_summary, outdir)

    log('[6/7] Freezing headline comparison and scientific decision')
    headline = []
    headline.append({'quantity':'broad_retained_set_jaccard','historical':1.0,'corrected':float(broad_acc_row['jaccard']),'delta':float(broad_acc_row['jaccard'])-1.0})
    if witness_info:
        headline.append({'quantity':'witness_37318_calibration_score','historical':witness_info['historical_score'],'corrected':witness_info['corrected_score'],'delta':witness_info['corrected_score']-witness_info['historical_score']})
    if strong_primary:
        headline.append({'quantity':'retained_broad_strong_n_primary_threshold','historical':strong_primary['historical_accepted_strong_n'],'corrected':strong_primary['corrected_accepted_strong_n'],'delta':strong_primary['corrected_accepted_strong_n']-strong_primary['historical_accepted_strong_n']})
    if step20_summary is not None and np.any(np.isclose(step20_summary['perturbation_pct'],5)):
        r = step20_summary[np.isclose(step20_summary['perturbation_pct'],5)].iloc[0]
        headline.append({'quantity':'step20_plusminus5_joint_fraction','historical':r['historical_strong_and_compatible_fraction'],'corrected':r['corrected_strong_and_compatible_fraction'],'delta':r['delta_joint_fraction']})
    if step21_summary is not None and np.any(np.isclose(step21_summary['factor'],1.5)):
        r = step21_summary[np.isclose(step21_summary['factor'],1.5)].iloc[0]
        headline.append({'quantity':'step21_factor1p5_joint_fraction','historical':r['historical_strong_and_compatible_fixed_fraction'],'corrected':r['corrected_strong_and_compatible_fixed_fraction'],'delta':r['delta_fixed_joint_fraction']})
    headline_df = pd.DataFrame(headline)
    headline_df.to_csv(outdir / '07_headline_comparison.csv', index=False)

    criteria = {}
    criteria['broad_retained_jaccard_below_0p90'] = bool(float(broad_acc_row['jaccard']) < 0.90)
    criteria['primary_strong_retained_count_changed'] = bool(strong_primary and int(strong_primary['historical_accepted_strong_n']) != int(strong_primary['corrected_accepted_strong_n']))
    criteria['witness_loses_compatibility'] = bool(witness_info is not None and witness_info['historical_accepted'] and not witness_info['corrected_accepted'])
    if step20_summary is not None and np.any(np.isclose(step20_summary['perturbation_pct'],5)):
        r = step20_summary[np.isclose(step20_summary['perturbation_pct'],5)].iloc[0]
        criteria['step20_plusminus5_joint_shift_gt_0p05'] = bool(abs(float(r['delta_joint_fraction'])) > 0.05)
    else:
        criteria['step20_plusminus5_joint_shift_gt_0p05'] = None
    if step21_summary is not None and np.any(np.isclose(step21_summary['factor'],1.5)):
        r = step21_summary[np.isclose(step21_summary['factor'],1.5)].iloc[0]
        criteria['step21_factor1p5_joint_shift_gt_0p05'] = bool(abs(float(r['delta_fixed_joint_fraction'])) > 0.05)
    else:
        criteria['step21_factor1p5_joint_shift_gt_0p05'] = None
    material = any(v is True for v in criteria.values())
    complete_full = not broad_only and {'broad','reference'}.issubset(set(d['prior'].unique()))
    decision = {
        'status': 'MATERIAL_CHANGE' if material else 'ROBUST_TO_WINDOW_CORRECTION',
        'full_two_prior_input': complete_full,
        'publication_ready_decision': bool(complete_full),
        'broad_only_run_allowed_by_flag': bool(broad_only and args.allow_broad_only),
        'criteria': criteria,
        'corrected_broad_cutoff_midpoint': new_cut,
        'historical_broad_cutoff_midpoint': old_cut,
        'required_follow_up': [
            'Use corrected 140-500 ms model integral in the calibration score.',
            'Regenerate Step19 calibration-frontier reporting from corrected scores.',
            'Regenerate any Step20/Step21 compatibility-boundary results that depend on the calibration cutoff.',
            'Rerun Step22 component likelihood/model comparison if compatible support changes materially; the count-based experimental occupancy update itself is unchanged.',
            'Regenerate Step23 publication consolidation and revise the manuscript before submission.'
        ] if material else [
            'Update Methods/Eq.7 definitions and provenance to use the corrected 140-500 ms model integral.',
            'Archive this Step24 audit as a sensitivity/correction check.'
        ],
        'note': 'The pipeline does not alter experimental labels or raw-data measurements.'
    }
    (outdir / '08_scientific_decision.json').write_text(json.dumps(decision, indent=2) + '\n')

    lines = [
        '# Step24 results', '',
        f'- Historical score replay: **{replay_status}**; max absolute error {maxerr:.3g}.',
        f'- Broad retained-set Jaccard after correcting the model integration window: **{float(broad_acc_row["jaccard"]):.4f}**.',
        f'- Historical broad cutoff: {old_cut:.6g}; corrected broad cutoff: {new_cut:.6g}.',
    ]
    if witness_info:
        lines.append(f'- Witness {args.witness_id}: score {witness_info["historical_score"]:.6g} -> {witness_info["corrected_score"]:.6g}; corrected retained = {witness_info["corrected_accepted"]}.')
    if strong_primary:
        lines.append(f'- Strong candidates among retained broad sets at r >= {args.strong_threshold:.6g}: {int(strong_primary["historical_accepted_strong_n"])} -> **{int(strong_primary["corrected_accepted_strong_n"])}**.')
    if broad_only:
        lines.append('- WARNING: full Step10 two-prior table was not found; this run is a broad-only audit and should not be used as the final repository freeze unless explicitly accepted.')
    lines.extend(['', f'**Scientific decision: {decision["status"]}.**', '', 'See `08_scientific_decision.json` and the detailed CSV outputs.'])
    (outdir / 'README_RESULTS.md').write_text('\n'.join(lines) + '\n')

    log('[7/7] Done')
    print(json.dumps({'status':'PASS','scientific_decision':decision['status'],'outdir':str(outdir),'broad_jaccard':float(broad_acc_row['jaccard']),'corrected_broad_cutoff':new_cut}, indent=2))


if __name__ == '__main__':
    main()
