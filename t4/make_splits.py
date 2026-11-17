"""Gene-level 70/15/15 split, stratified by surface status, grouped by genomic block.

Change from CP1: genes are kept together in 10 Mb genomic blocks. Neighbouring genes
share GISTIC amplification calls (same amplicon), so a random gene split would put
near-identical feature rows in train and test. Grouping by block prevents that leak.

    python -m t4.make_splits --config configs/t4_cp2.yaml
"""
import argparse
import pandas as pd
from sklearn.model_selection import StratifiedGroupKFold
from .common import load_config, write_json


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="configs/t4_cp2.yaml")
    cfg = load_config(ap.parse_args().config)
    sc = cfg["split"]

    g = pd.read_csv(cfg["processed_dir"] / "gene_table.csv")
    g = g[g["dependency"].notna()].copy()
    g["block"] = g["chrom"].astype(str) + ":" + (g["start"] // (sc["block_mb"] * 1_000_000)).astype(int).astype(str)

    # 20 folds of 5% each -> 14 train, 3 validation, 3 test.
    k = 20
    n_tr, n_va = round(sc["fractions"][0] * k), round(sc["fractions"][1] * k)
    sgkf = StratifiedGroupKFold(n_splits=k, shuffle=True, random_state=sc["seed"])
    fold = pd.Series(-1, index=g.index)
    for i, (_, idx) in enumerate(sgkf.split(g, g["is_surface"], groups=g["block"])):
        fold.iloc[idx] = i
    g["split"] = pd.cut(fold, [-1, n_tr - 1, n_tr + n_va - 1, k], labels=["train", "val", "test"]).astype(str)

    assert g.groupby("block")["split"].nunique().max() == 1, "a block crosses splits"
    counts = (g.groupby("split").agg(genes=("gene", "size"), surface_genes=("is_surface", "sum"),
                                     blocks=("block", "nunique"), label_mean=("dependency", "mean"))
                .reindex(["train", "val", "test"]))
    counts["surface_share"] = counts["surface_genes"] / counts["genes"]
    assert counts.loc["test", "surface_genes"] >= sc["min_test_surface"], "too few test surface genes"

    g[["gene", "split", "block", "chrom", "start", "is_surface"]].to_csv(
        cfg["splits_dir"] / "split_manifest.csv", index=False)
    write_json({"seed": sc["seed"], "unit": "gene (grouped by %d Mb block)" % sc["block_mb"],
                "counts": counts.round(4).to_dict(orient="index")},
               cfg["splits_dir"] / "split_counts.json")
    print(counts.round(3).to_string())


if __name__ == "__main__":
    main()
