"""Evaluate frozen held-out predictions; see docs/zachary_evaluation.md."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

MODELS = ('constant', 'ridge', 'random_forest', 'gradient_boosting')


def metric(y, p):
    mae = float(np.mean(np.abs(y - p)))
    rho = float(spearmanr(y, p).statistic) if len(np.unique(y)) > 1 and len(np.unique(p)) > 1 else np.nan
    return mae, rho


def interval(values):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if len(values) < 100:
        return (np.nan, np.nan)
    return tuple(np.percentile(values, [2.5, 97.5]))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--predictions', required=True)
    ap.add_argument('--frozen', required=True, help='Pre-test frozen JSON decision record')
    ap.add_argument('--out', required=True)
    ap.add_argument('--bootstrap', type=int, default=2000)
    ap.add_argument('--seed', type=int, default=42)
    ap.add_argument('--scores', help='CSV gene,score produced using the confirmed course guide')
    args = ap.parse_args()
    if args.bootstrap < 100:
        ap.error('at least 100 bootstrap samples are required')
    frozen = json.loads(Path(args.frozen).read_text(encoding='utf-8'))
    primary = frozen['primary_model']
    if primary not in MODELS[1:] or frozen['primary_metric'] != 'overall_spearman':
        raise ValueError('primary must be a trained model and metric overall_spearman')
    if not frozen.get('frozen_before_test') or not frozen.get('tolerance_justification'):
        raise ValueError('record pre-test freeze and a scientific tolerance justification')
    target, tolerance = float(frozen['target']), float(frozen['tolerance'])
    if not np.isfinite(target) or not np.isfinite(tolerance) or tolerance < 0:
        raise ValueError('target and nonnegative tolerance must be finite')
    df = pd.read_csv(args.predictions)
    required = ['gene', 'split', 'block', 'is_surface', 'dependency'] + ['pred_' + m for m in MODELS] + ['pred_primary_no_amplification']
    missing = set(required) - set(df.columns)
    if missing:
        raise ValueError(f'Missing columns: {sorted(missing)}')
    if df.gene.isna().any() or df.gene.duplicated().any():
        raise ValueError('gene IDs must be present and unique')
    if not df['split'].eq('test').all():
        raise ValueError('input must contain test genes only')
    if df.block.isna().any() or not df.is_surface.isin([0, 1]).all():
        raise ValueError('block and binary is_surface must be present')
    numeric = ['dependency'] + ['pred_' + m for m in MODELS] + ['pred_primary_no_amplification']
    df[numeric] = df[numeric].apply(pd.to_numeric, errors='raise')
    if not np.isfinite(df[numeric].to_numpy()).all():
        raise ValueError('all labels and predictions must be finite on the same test genes')
    if df.is_surface.sum() < 20:
        raise ValueError('at least 20 held-out surface candidates required')
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    rows, boot = [], {}
    for subset, part in [('overall', df), ('surface', df[df.is_surface.eq(1)])]:
        part = part.reset_index(drop=True)
        groups = [g.index.to_numpy() for _, g in part.groupby('block', sort=True)]
        if len(groups) < 2:
            raise ValueError(f'{subset} has fewer than two independent blocks')
        samples = [np.concatenate([groups[i] for i in rng.integers(0, len(groups), len(groups))]) for _ in range(args.bootstrap)]
        y = part.dependency.to_numpy()
        for name in (*MODELS, 'primary_no_amplification'):
            p = part['pred_' + name].to_numpy()
            mae, rho = metric(y, p)
            draws = np.array([metric(y[idx], p[idx]) for idx in samples])
            boot[subset, name] = draws
            ml, mh = interval(draws[:, 0]); rl, rh = interval(draws[:, 1])
            rows.append(dict(subset=subset, model=name, n=len(part), blocks=len(groups), mae=mae, mae_ci_low=ml, mae_ci_high=mh, spearman=rho, spearman_ci_low=rl, spearman_ci_high=rh))
    summary = pd.DataFrame(rows)
    summary.to_csv(out / 'evaluation_metrics.csv', index=False)
    comparisons = []
    for subset in ('overall', 'surface'):
        for other in [m for m in MODELS[1:] if m != primary] + ['primary_no_amplification']:
            a, b = boot[subset, primary], boot[subset, other]
            x = summary.set_index(['subset', 'model'])
            # Positive MAE and Spearman differences favor the primary model.
            delta_mae = x.loc[(subset, other), 'mae'] - x.loc[(subset, primary), 'mae']
            delta_rho = x.loc[(subset, primary), 'spearman'] - x.loc[(subset, other), 'spearman']
            ml, mh = interval(b[:, 0] - a[:, 0]); rl, rh = interval(a[:, 1] - b[:, 1])
            comparisons.append(dict(subset=subset, primary=primary, comparator=other, delta_mae=delta_mae, mae_ci_low=ml, mae_ci_high=mh, delta_spearman=delta_rho, spearman_ci_low=rl, spearman_ci_high=rh))
    pd.DataFrame(comparisons).to_csv(out / 'paired_comparisons.csv', index=False)
    overall = summary[summary.subset.eq('overall')].set_index('model')
    candidates = overall.loc[[m for m in MODELS[1:] if m != primary], 'spearman']
    if not np.isfinite(candidates).all() or not np.isfinite(overall.loc[primary, 'spearman']):
        raise ValueError('trained model Spearman is undefined; comparability cannot be assessed')
    strongest = candidates.idxmax()
    value = float(overall.loc[primary, 'spearman'])
    gap = float(candidates.max() - value)
    ci_low = float(overall.loc[primary, 'spearman_ci_low'])
    verdict = dict(primary_model=primary, primary_metric='overall_spearman', primary_value=value, strongest_comparator=strongest, comparator_value=float(candidates.max()), gap=gap, tolerance=tolerance, comparable=gap <= tolerance, target=target, target_met=bool(value >= target and np.isfinite(ci_low) and ci_low > 0), bootstrap_unit='genomic block', bootstrap_samples=args.bootstrap, bootstrap_seed=args.seed)
    (out / 'comparability_verdict.json').write_text(json.dumps(verdict, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    if args.scores:
        if not frozen.get('ranking_formula_source'):
            raise ValueError('ranking requires the confirmed guide formula source in frozen JSON')
        scores = pd.read_csv(args.scores)
        if not {'gene', 'score'}.issubset(scores.columns) or scores.gene.isna().any() or scores.gene.duplicated().any():
            raise ValueError('scores must have unique gene and score columns')
        ranked = df[df.is_surface.eq(1)].merge(scores[['gene', 'score']], on='gene', how='left', validate='one_to_one')
        ranked['score'] = pd.to_numeric(ranked.score, errors='raise')
        if not np.isfinite(ranked.score).all():
            raise ValueError('every eligible gene requires a finite guide score')
        ranked = ranked.sort_values(['score', 'gene'], ascending=[False, True])
        ranked.insert(0, 'rank', np.arange(1, len(ranked) + 1))
        ranked.to_csv(out / 'surface_candidates.csv', index=False)


if __name__ == '__main__':
    main()
