"""Figures for the report and slides, from files in results/.

    python -m t4.make_figures --config configs/t4_cp2.yaml
"""
import argparse, json
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import spearmanr
from .common import load_config

LABEL = {"constant": "Constant", "ridge": "Ridge", "random_forest": "Random forest",
         "gradient_boosting": "Gradient boosting"}
ORDER = ["constant", "ridge", "random_forest", "gradient_boosting"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="configs/t4_cp2.yaml")
    cfg = load_config(ap.parse_args().config)
    res = cfg["results_dir"]
    fig_dir = res / "figures"
    fig_dir.mkdir(exist_ok=True)
    frozen = json.load(open(res / "frozen_choices.json"))
    fs, pm = frozen["primary"]["feature_set"], frozen["primary"]["model"]
    summ = pd.read_csv(res / "test_summary.csv")
    s = summ[summ["feature_set"] == fs].set_index("model").loc[ORDER]

    # 1. test Spearman and MAE with 95% bootstrap CIs, all models, primary feature set
    fig, ax = plt.subplots(1, 2, figsize=(8, 3))
    x = range(len(ORDER))
    colors = ["#999999" if m == "constant" else ("#1f77b4" if m == pm else "#9ecae1") for m in ORDER]
    for a, met in zip(ax, ["spearman", "mae"]):
        mid = s[f"{met}_all_mean"]
        err = [mid - s[f"{met}_ci_lo"], s[f"{met}_ci_hi"] - mid]
        a.bar(x, mid, yerr=err, color=colors, capsize=3)
        a.set_xticks(list(x), [LABEL[m] for m in ORDER], rotation=20, ha="right", fontsize=8)
        a.set_title(f"Test {'Spearman ρ' if met == 'spearman' else 'MAE'} (all genes)", fontsize=9)
    ax[0].axhline(cfg["evaluation"]["success_target"], ls="--", c="k", lw=0.8)
    fig.tight_layout()
    fig.savefig(fig_dir / "test_comparison.png", dpi=200)

    # 2. ablation / improvement: Spearman by feature set per trained model
    fig, a = plt.subplots(figsize=(5, 3))
    piv = summ.pivot(index="model", columns="feature_set", values="spearman_all_mean").loc[ORDER[1:]]
    piv[["no_amp", "base", "interactions"]].plot.bar(ax=a, rot=0, color=["#fdae6b", "#9ecae1", "#1f77b4"])
    a.set_xticklabels([LABEL[m] for m in piv.index], fontsize=8)
    a.set_ylabel("Test Spearman ρ"); a.set_xlabel("")
    a.legend(["No amplification", "Base (CP1)", "+ Interactions"], fontsize=7)
    fig.tight_layout()
    fig.savefig(fig_dir / "ablation.png", dpi=200)

    # 3. predicted vs measured dependency for the primary model
    pr = pd.read_csv(res / "test_predictions.csv")
    col = f"{pm}|{fs}"
    fig, a = plt.subplots(figsize=(4, 3.5))
    for surf, c, lab in ((0, "#c6dbef", "non-surface"), (1, "#d62728", "surface")):
        d = pr[pr["is_surface"] == surf]
        a.scatter(d[col], d["dependency"], s=4, c=c, label=lab, alpha=0.7)
    rho = spearmanr(pr[col], pr["dependency"]).statistic
    a.set_xlabel("Predicted dependency"); a.set_ylabel("Measured (−median Chronos)")
    a.set_title(f"{LABEL[pm]}, test genes (ρ = {rho:.2f})", fontsize=9)
    a.legend(fontsize=7, markerscale=3)
    fig.tight_layout()
    fig.savefig(fig_dir / "pred_vs_true.png", dpi=200)
    print("figures written to", fig_dir)


if __name__ == "__main__":
    main()
