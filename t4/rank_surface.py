"""Rank held-out (test) surface candidates (METHODS_AND_EQUATIONS A7).

Eligibility: HPA plasma-membrane location (is_surface = 1) and in the test split.
Score: weighted percentile ranks fixed in the config before test scoring.
Predicted dependency from the frozen primary model is supporting evidence only;
a nonessential surface antigen can still rank high on expression evidence.

    python -m t4.rank_surface --config configs/t4_cp2.yaml
"""
import argparse, json
import joblib
import pandas as pd
from .common import load_config


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="configs/t4_cp2.yaml")
    cfg = load_config(ap.parse_args().config)
    res = cfg["results_dir"]
    frozen = json.load(open(res / "frozen_choices.json"))
    p = frozen["primary"]
    bundle = joblib.load(cfg["models_dir"] / f"{p['model']}__{p['feature_set']}.joblib")

    g = pd.read_csv(cfg["processed_dir"] / "gene_table.csv")
    man = pd.read_csv(cfg["splits_dir"] / "split_manifest.csv")[["gene", "split"]]
    t = g.merge(man, on="gene")
    t = t[(t["split"] == "test") & (t["is_surface"] == 1)].copy()
    t["predicted_dependency"] = bundle["model"].predict(t[bundle["features"]])

    w = frozen["ranking_weights"]
    t["priority_score"] = sum(wt * t[k].rank(pct=True) for k, wt in w.items())
    t = t.sort_values("priority_score", ascending=False)
    t["rank"] = range(1, len(t) + 1)
    t["essential_like"] = (t["predicted_dependency"] >= 0.5).map({True: "yes", False: "no"})
    cols = ["rank", "gene", "priority_score", "expr_enrichment", "expr_mean_tumor", "gistic_amp_freq",
            "predicted_dependency", "dependency", "essential_like", "hpa_reliability"]
    t[cols].to_csv(res / "surface_candidates_test.csv", index=False, float_format="%.4g")
    print(f"{len(t)} eligible held-out surface genes; top {cfg['ranking']['top_n']}:")
    print(t[cols].head(cfg["ranking"]["top_n"]).round(3).to_string(index=False))
    print("Note: 'dependency' (measured) is shown for evaluation only; it is not used in the score.")


if __name__ == "__main__":
    main()
