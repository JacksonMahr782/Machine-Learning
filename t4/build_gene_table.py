"""Build the one-row-per-gene table (METHODS_AND_EQUATIONS A1, A2).

Features come only from Xena expression and cBioPortal GISTIC copy number.
DepMap supplies only the label. HPA supplies only the surface flag (not a feature).

    python -m t4.build_gene_table --config configs/t4_cp2.yaml
"""
import argparse
import numpy as np
import pandas as pd
from .common import load_config, raw, write_json, sha256


def tcga_sample_type(barcode):
    """TCGA barcode chars 14-15: 01-09 tumour, 10-19 normal."""
    code = int(barcode.split("-")[3][:2])
    return "tumor" if code < 10 else ("normal" if code < 20 else "other")


def expression_features(cfg, log):
    pm = pd.read_csv(raw(cfg, "xena_probemap"), sep="\t")
    pm = pm.rename(columns={"id": "ensembl_id", "gene": "gene", "chromStart": "start"})
    tpm = pd.read_csv(raw(cfg, "xena_tpm"), sep="\t", index_col=0, dtype={0: str})
    tpm = tpm.astype(np.float32)  # values are already log2(TPM + 1) on the Xena GDC hub
    types = pd.Series({s: tcga_sample_type(s) for s in tpm.columns})
    T, N = types.index[types == "tumor"], types.index[types == "normal"]
    log["xena_genes_raw"], log["xena_samples_tumor"], log["xena_samples_normal"] = len(tpm), len(T), len(N)
    log["xena_patients_tumor"] = int(pd.Index(T).str[:12].nunique())
    log["xena_patients_normal"] = int(pd.Index(N).str[:12].nunique())
    f = pd.DataFrame({
        "expr_mean_tumor": tpm[T].mean(axis=1),
        "expr_mean_normal": tpm[N].mean(axis=1),
        "expr_sd_tumor": tpm[T].std(axis=1, ddof=1),
    })
    f["expr_enrichment"] = f["expr_mean_tumor"] - f["expr_mean_normal"]
    f = f.join(pm.set_index("ensembl_id")[["gene", "chrom", "start"]], how="inner")
    f = f[~f.index.str.endswith("_PAR_Y")]
    # Several Ensembl IDs can share a symbol: keep the most expressed one.
    f = (f.sort_values("expr_mean_tumor", ascending=False)
          .reset_index(names="ensembl_id").drop_duplicates("gene").set_index("gene"))
    log["xena_genes_unique_symbol"] = len(f)
    return f


def cna_features(cfg, log):
    cna = pd.read_csv(raw(cfg, "cbio_cna"), sep="\t", low_memory=False)
    cna = cna.dropna(subset=["Hugo_Symbol"]).drop_duplicates("Hugo_Symbol").set_index("Hugo_Symbol")
    calls = cna.drop(columns=[c for c in ("Entrez_Gene_Id", "Cytoband") if c in cna.columns])
    calls = calls.apply(pd.to_numeric, errors="coerce")
    n = calls.notna().sum(axis=1)
    log["cna_genes"], log["cna_samples"] = len(calls), calls.shape[1]
    return pd.DataFrame({
        "gistic_amp_freq": (calls >= 1).sum(axis=1) / n,
        "gistic_highamp_freq": (calls == 2).sum(axis=1) / n,
    })


def dependency_label(cfg, log):
    model = pd.read_csv(raw(cfg, "depmap_model"), low_memory=False)
    luad = set(model.loc[model["OncotreeCode"] == cfg["label"]["oncotree_code"], "ModelID"])
    eff = pd.read_csv(raw(cfg, "depmap_effect"), index_col=0)
    eff = eff.loc[eff.index.isin(luad)]
    eff.columns = [c.split(" (")[0] for c in eff.columns]
    eff = eff.loc[:, ~eff.columns.duplicated()]
    n_lines = eff.notna().sum(axis=0)
    y = -eff.median(axis=0, skipna=True)
    y[n_lines < cfg["label"]["min_lines"]] = np.nan
    log["depmap_luad_models_in_model_csv"], log["depmap_luad_lines_with_effects"] = len(luad), len(eff)
    log["depmap_genes"], log["depmap_genes_labeled"] = eff.shape[1], int(y.notna().sum())
    return pd.DataFrame({"dependency": y, "n_luad_lines": n_lines})


def surface_flag(cfg, log):
    path = raw(cfg, "hpa_location")
    if not path.exists() and path.with_suffix(".tsv.zip").exists():
        path = path.with_suffix(".tsv.zip")
    hpa = pd.read_csv(path, sep="\t")
    loc = hpa["Main location"].fillna("") + ";" + hpa["Additional location"].fillna("")
    hpa["is_surface"] = loc.str.contains("Plasma membrane", regex=False).astype(int)
    hpa = hpa.drop_duplicates("Gene name").set_index("Gene name")
    log["hpa_genes"], log["hpa_surface_genes"] = len(hpa), int(hpa["is_surface"].sum())
    return hpa[["is_surface", "Reliability"]].rename(columns={"Reliability": "hpa_reliability"})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", default="configs/t4_cp2.yaml")
    cfg = load_config(ap.parse_args().config)
    log = {"inputs": {k: {"file": v, "sha256": sha256(cfg["raw_dir"] / v)}
                      for k, v in cfg["files"].items() if (cfg["raw_dir"] / v).exists()}}

    expr = expression_features(cfg, log)
    cna = cna_features(cfg, log)
    lab = dependency_label(cfg, log)
    hpa = surface_flag(cfg, log)

    g = expr.join(cna, how="inner")              # genes with every feature
    log["genes_with_all_features"] = len(g)
    g = g.join(lab, how="left").join(hpa, how="left")
    g["hpa_annotated"] = g["is_surface"].notna().astype(int)
    g["is_surface"] = g["is_surface"].fillna(0).astype(int)
    g["mean_x_amp"] = g["expr_mean_tumor"] * g["gistic_amp_freq"]
    g["enrich_x_amp"] = g["expr_enrichment"] * g["gistic_amp_freq"]
    g.index.name = "gene"

    lab_mask = g["dependency"].notna()
    log["genes_labeled"] = int(lab_mask.sum())
    log["genes_labeled_surface"] = int((lab_mask & (g["is_surface"] == 1)).sum())
    log["label_summary"] = g.loc[lab_mask, "dependency"].describe().round(4).to_dict()

    out = cfg["processed_dir"] / "gene_table.csv"
    g.reset_index().to_csv(out, index=False, float_format="%.6g")
    write_json(log, cfg["processed_dir"] / "build_log.json")
    print(f"wrote {out}: {len(g)} genes, {log['genes_labeled']} labeled, "
          f"{log['genes_labeled_surface']} labeled surface")


if __name__ == "__main__":
    main()
