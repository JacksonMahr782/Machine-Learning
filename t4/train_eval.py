"""Train, tune on validation, freeze, then evaluate once on held-out test genes.

Models (METHODS_AND_EQUATIONS A4): constant, Ridge, random forest, gradient boosting.
Feature sets (A6): base, interactions (improvement), no_amp (amplification ablation).
Every trained model gets the same budget: 12 validation-scored configurations per feature set.

    python -m t4.train_eval --config configs/t4_cp2.yaml
"""
import argparse, itertools, time, platform, datetime
import numpy as np
import pandas as pd
import joblib
import sklearn
from scipy.stats import spearmanr
from sklearn.dummy import DummyRegressor
from sklearn.linear_model import Ridge
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from .common import load_config, write_json, sha256, FEATURE_SETS

MODELS = ["constant", "ridge", "random_forest", "gradient_boosting"]
TRAINED = MODELS[1:]


def build(name, params, seed):
    if name == "constant":
        return DummyRegressor(strategy="mean")
    if name == "ridge":  # scaler is fit on training genes only (A3)
        return make_pipeline(StandardScaler(), Ridge(**params))
    if name == "random_forest":
        return RandomForestRegressor(random_state=seed, n_jobs=-1, **params)
    return GradientBoostingRegressor(random_state=seed, loss="squared_error", **params)


def grid(cfg, name):
    if name == "constant":
        return [{}]
    spec = cfg["grids"][name]
    keys = list(spec)
    return [dict(zip(keys, v)) for v in itertools.product(*(spec[k] for k in keys))]


def metrics(y, p, surf):
    out = {}
    for tag, m in (("all", np.ones_like(surf, bool)), ("surface", surf == 1)):
        yt, pt = y[m], p[m]
        out[f"mae_{tag}"] = float(np.mean(np.abs(yt - pt)))
        out[f"mse_{tag}"] = float(np.mean((yt - pt) ** 2))
        out[f"spearman_{tag}"] = float(spearmanr(yt, pt).statistic) if np.std(pt) > 0 else 0.0
        out[f"n_{tag}"] = int(m.sum())
    return out


def bootstrap(y, preds, surf, primary, B, seed):
    """Percentile 95% CIs per model plus paired differences vs the primary model."""
    rng = np.random.default_rng(seed)
    n = len(y)
    stats = {k: {"mae": [], "sp": []} for k in preds}
    for _ in range(B):
        i = rng.integers(0, n, n)
        for k, p in preds.items():
            stats[k]["mae"].append(np.mean(np.abs(y[i] - p[i])))
            stats[k]["sp"].append(spearmanr(y[i], p[i]).statistic if np.std(p[i]) > 0 else 0.0)
    rows = []
    for k in preds:
        mae, sp = np.array(stats[k]["mae"]), np.nan_to_num(np.array(stats[k]["sp"]))
        d_mae = mae - np.array(stats[primary]["mae"])        # >0: primary has lower error
        d_sp = np.nan_to_num(np.array(stats[primary]["sp"])) - sp  # >0: primary ranks better
        rows.append({"model": k,
                     "mae_ci_lo": np.percentile(mae, 2.5), "mae_ci_hi": np.percentile(mae, 97.5),
                     "spearman_ci_lo": np.percentile(sp, 2.5), "spearman_ci_hi": np.percentile(sp, 97.5),
                     "primary_minus_model_mae_ci": [np.percentile(-d_mae, 2.5), np.percentile(-d_mae, 97.5)],
                     "primary_minus_model_spearman_ci": [np.percentile(d_sp, 2.5), np.percentile(d_sp, 97.5)],
                     "p_primary_not_better_spearman": float(np.mean(d_sp <= 0)) if k != primary else None})
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="configs/t4_cp2.yaml")
    cfg = load_config(ap.parse_args().config)
    res, mdir = cfg["results_dir"], cfg["models_dir"]
    t_start = time.time()

    g = pd.read_csv(cfg["processed_dir"] / "gene_table.csv")
    man = pd.read_csv(cfg["splits_dir"] / "split_manifest.csv")[["gene", "split"]]
    g = g.merge(man, on="gene", how="inner")
    S = {s: g[g["split"] == s].reset_index(drop=True) for s in ("train", "val", "test")}
    y = {s: S[s]["dependency"].to_numpy() for s in S}
    surf = {s: S[s]["is_surface"].to_numpy() for s in S}

    # ---------- 1. validation: hyper-parameters and feature set -------------------
    rows = []
    for fs, cols in FEATURE_SETS.items():
        for name in MODELS:
            for params in grid(cfg, name):
                t0 = time.time()
                m = build(name, params, cfg["seeds"][0]).fit(S["train"][cols], y["train"])
                fit_s = time.time() - t0
                rows.append({"feature_set": fs, "model": name, "params": params, "fit_seconds": fit_s,
                             **metrics(y["val"], m.predict(S["val"][cols]), surf["val"])})
                print(f"[val] {fs:12s} {name:18s} {params} spearman={rows[-1]['spearman_all']:.3f}")
    val = pd.DataFrame(rows)
    val.to_csv(res / "validation_grid.csv", index=False)

    best = (val.sort_values("mse_all").groupby(["feature_set", "model"], as_index=False).first())
    cand = best[best["model"].isin(TRAINED) & best["feature_set"].isin(["base", "interactions"])]
    prim = cand.sort_values("spearman_all", ascending=False).iloc[0]
    frozen = {
        "frozen_at_utc": datetime.datetime.utcnow().isoformat(timespec="seconds"),
        "rule": "hyper-parameters: lowest validation MSE per (model, feature set); "
                "primary: highest validation Spearman among trained models on base/interactions",
        "primary": {"model": prim["model"], "feature_set": prim["feature_set"]},
        "configs": {f"{r.model}|{r.feature_set}": r.params for r in best.itertuples()},
        "evaluation": cfg["evaluation"],
        "ranking_weights": cfg["ranking"]["weights"],
        "split_manifest_sha256": sha256(cfg["splits_dir"] / "split_manifest.csv"),
        "gene_table_sha256": sha256(cfg["processed_dir"] / "gene_table.csv"),
    }
    write_json(frozen, res / "frozen_choices.json")   # written BEFORE any test prediction
    print("FROZEN:", frozen["primary"])
    pkey = f"{prim['model']}|{prim['feature_set']}"

    # ---------- 2. test: evaluated once, every model on the same test genes -------
    rows, preds0, pred_table = [], {}, S["test"][["gene", "is_surface", "dependency"]].copy()
    for r in best.itertuples():
        cols, key = FEATURE_SETS[r.feature_set], f"{r.model}|{r.feature_set}"
        seeds = cfg["seeds"] if r.model in ("random_forest", "gradient_boosting") else cfg["seeds"][:1]
        for sd in seeds:
            t0 = time.time()
            m = build(r.model, r.params, sd).fit(S["train"][cols], y["train"])
            fit_s = time.time() - t0
            p = m.predict(S["test"][cols])
            rows.append({"model": r.model, "feature_set": r.feature_set, "seed": sd,
                         "fit_seconds": fit_s, **metrics(y["test"], p, surf["test"])})
            if sd == cfg["seeds"][0]:
                preds0[key] = p
                pred_table[key] = p
                joblib.dump({"model": m, "features": cols, "params": r.params},
                            mdir / f"{r.model}__{r.feature_set}.joblib")
    test = pd.DataFrame(rows)
    test.to_csv(res / "test_metrics_by_seed.csv", index=False)
    pred_table.to_csv(res / "test_predictions.csv", index=False)

    num = [c for c in test.columns if c.startswith(("mae", "mse", "spearman", "fit"))]
    summ = test.groupby(["model", "feature_set"])[num].agg(["mean", "std"])
    summ.columns = [f"{a}_{b}" for a, b in summ.columns]
    summ = summ.reset_index()

    # bootstrap CIs on the primary feature set (comparison table) and on every set
    ci = bootstrap(y["test"], preds0, surf["test"], pkey, cfg["bootstrap"], cfg["split"]["seed"])
    ci[["model", "feature_set"]] = ci["model"].str.split("|", expand=True)
    summ = summ.merge(ci, on=["model", "feature_set"])
    summ["is_primary"] = (summ["model"] == prim["model"]) & (summ["feature_set"] == prim["feature_set"])
    summ.to_csv(res / "test_summary.csv", index=False)

    # amplification ablation and improvement, per model (A6)
    piv = summ.pivot(index="model", columns="feature_set", values=["mae_all_mean", "spearman_all_mean"])
    abl = pd.DataFrame({
        "delta_mae_no_amp_minus_base": piv[("mae_all_mean", "no_amp")] - piv[("mae_all_mean", "base")],
        "delta_spearman_base_minus_no_amp": piv[("spearman_all_mean", "base")] - piv[("spearman_all_mean", "no_amp")],
        "delta_mae_base_minus_interactions": piv[("mae_all_mean", "base")] - piv[("mae_all_mean", "interactions")],
        "delta_spearman_interactions_minus_base": piv[("spearman_all_mean", "interactions")] - piv[("spearman_all_mean", "base")],
    }).reset_index()
    abl.to_csv(res / "ablation_improvement.csv", index=False)

    # comparability standard and success target (declared in config, frozen above)
    ev = cfg["evaluation"]
    pm = summ.set_index(["model", "feature_set"])
    prim_score = pm.loc[(prim["model"], prim["feature_set"]), ev["primary_metric"] + "_mean"]
    comps = summ[summ["model"].isin(TRAINED) & (summ["feature_set"] == prim["feature_set"])
                 & (summ["model"] != prim["model"])]
    strongest = comps.sort_values(ev["primary_metric"] + "_mean", ascending=False).iloc[0]
    gap = prim_score - strongest[ev["primary_metric"] + "_mean"]
    lo = pm.loc[(prim["model"], prim["feature_set"]), "spearman_ci_lo"]
    verdict = {"primary": pkey, "primary_metric": ev["primary_metric"], "primary_value": prim_score,
               "success_target": ev["success_target"], "target_met": bool(prim_score >= ev["success_target"] and lo > 0),
               "primary_ci95": [lo, pm.loc[(prim["model"], prim["feature_set"]), "spearman_ci_hi"]],
               "strongest_comparator": f"{strongest['model']}|{strongest['feature_set']}",
               "comparator_value": strongest[ev["primary_metric"] + "_mean"], "gap": gap,
               "tolerance": ev["tolerance"], "comparability_met": bool(gap >= -ev["tolerance"]),
               "paired_bootstrap_gap_ci95": strongest["primary_minus_model_spearman_ci"]}
    write_json(verdict, res / "comparability_verdict.json")
    print("\nVERDICT", verdict)

    write_json({"python": platform.python_version(), "sklearn": sklearn.__version__,
                "platform": platform.platform(), "processor": platform.processor(),
                "total_seconds": round(time.time() - t_start, 1),
                "n_train": len(S["train"]), "n_val": len(S["val"]), "n_test": len(S["test"])},
               res / "run_info.json")

    show = summ[summ["feature_set"] == prim["feature_set"]][
        ["model", "mae_all_mean", "spearman_all_mean", "spearman_ci_lo", "spearman_ci_hi",
         "mae_surface_mean", "spearman_surface_mean"]]
    print("\nTEST (feature set = %s, primary = %s)\n" % (prim["feature_set"], prim["model"]),
          show.round(3).to_string(index=False))
    print("\nABLATION / IMPROVEMENT\n", abl.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
