# T4: Tumor surface target prioritization (CPSC 4430/6430, Fall 2026)

Group T4: Jackson Mahr, Dante Mallozzi, Veer Patel, Zachary Schoon.

**Repository:** https://github.com/JacksonMahr782/Machine-Learning
**Submitted revision (Checkpoint 2):** tag `cp2-submission`, commit `<FILL IN: git rev-parse cp2-submission>`

This repo holds two separate pieces of work:

1. **Paper reproduction** (TROP2, Vaquero-Siguero et al., Nature 2026): TACSTD2/LGR5 expression by cell state and donor cell-state proportions from processed E-MTAB-16585 epithelial data, using the authors' `figures/Figure_01.ipynb`. See `paper_repro/` (owner: Dante).
2. **LUAD gene-dependency predictor** (this README's commands): one row per gene, public expression and copy-number features, label = −median DepMap Chronos effect over ≥3 LUAD lines. Compares a constant control with Ridge, random forest and gradient boosting, then ranks held-out surface candidates.

Dependency predictions support prioritization only. They do not establish antigen validity or therapeutic efficacy.

## Folder layout

```
configs/t4_cp2.yaml        all paths, split rules, hyper-parameter grids, success target, ranking weights
t4/build_gene_table.py     raw files -> data/processed/gene_table.csv (+ build_log.json with counts, file hashes)
t4/make_splits.py          gene split -> data/splits/split_manifest.csv (+ split_counts.json)
t4/train_eval.py           validation grid, freeze, single test evaluation, bootstrap CIs, ablation
t4/rank_surface.py         held-out surface-candidate ranking
t4/make_figures.py         results/figures/*.png
run_all.sh                 runs the five steps in order
slurm/t4_cp2.sbatch        Palmetto job that runs run_all.sh
scripts/package_zip.py     builds the Canvas ZIP from the committed revision
scripts/make_data_table.py CP1 row/column counter for the raw files
data/raw/                  downloaded inputs (not in GitHub, see data/raw/README.md)
data/dictionaries/         CP1 source, gene/label and split dictionaries, download log
docs/                      CP1 METHODS_AND_EQUATIONS and DATA_TABLE
results/                   metrics, predictions, frozen choices, figures (committed after the run)
models/                    trained checkpoints, *.joblib (git-ignored; included in the ZIP)
```

## Setup

Python 3.11. Either conda (recommended on Palmetto):

```bash
conda env create -f environment.yml
conda activate t4cp2
```

or pip: `pip install -r requirements.txt`.

## Data access

All inputs are public. Download them into `data/raw/` with the exact names listed in `data/raw/README.md`
(URLs, releases and download dates are in `data/dictionaries/download_log.csv`). Total about 1.6 GB compressed.
Our copies are also on the group Box folder. The processed `data/processed/gene_table.csv` (with `build_log.json`) and `data/splits/` are committed, so training and ranking run without the raw files; rebuild them from `data/raw/` to reproduce from scratch.

| Input | Release | Used for |
|---|---|---|
| UCSC Xena GDC TCGA-LUAD `star_tpm` + gencode v36 probemap | version 05-10-2024 (GDC DR40) | expression features (tumour vs normal by TCGA barcode) |
| cBioPortal `luad_tcga_pan_can_atlas_2018` `data_cna.txt` | PanCancer Atlas 2018 | GISTIC amplification features |
| DepMap `CRISPRGeneEffect.csv`, `Model.csv` | Public 26Q3 | label only (OncotreeCode = LUAD) |
| HPA `subcellular_location.tsv` | v25 (Ensembl 109); main site is now v26, use the v25 archive | surface eligibility and split stratification only |

## Commands

From the repo root:

```bash
bash run_all.sh                      # whole pipeline
# or on Palmetto
sbatch slurm/t4_cp2.sbatch           # log in logs/t4_cp2_<jobid>.out
```

Individual steps (same order):

```bash
python -m t4.build_gene_table --config configs/t4_cp2.yaml
python -m t4.make_splits      --config configs/t4_cp2.yaml
python -m t4.train_eval       --config configs/t4_cp2.yaml
python -m t4.rank_surface     --config configs/t4_cp2.yaml
python -m t4.make_figures     --config configs/t4_cp2.yaml
```

Demo (uses the saved checkpoints, no retraining): `python -m t4.rank_surface --config configs/t4_cp2.yaml`.

## What each step does

* **Features** (`t4/build_gene_table.py`): from Xena log2(TPM+1), mean tumour expression, tumour − normal enrichment, tumour SD; from GISTIC, amplification (≥+1) and high-level amplification (=+2) frequency. Improvement adds two interaction terms (mean × amp, enrichment × amp). Symbols with several Ensembl IDs keep the most expressed one.
* **Label:** −median Chronos effect over LUAD lines with a measured value; genes with < 3 lines stay unlabeled. Nothing derived from DepMap is a feature.
* **Split** (`t4/make_splits.py`): labeled genes only, 70/15/15 train/validation/test, seed 42, stratified on HPA surface status. **Change from CP1:** genes are grouped into 10 Mb genomic blocks and a whole block goes to one split, because neighbouring genes share amplicons and would otherwise leak near-identical GISTIC features across splits. The script fails if any block crosses splits or the test set has < 20 surface genes.
* **Training** (`t4/train_eval.py`): MSE loss for every model. Ridge standardises features with training-gene statistics only. Each trained model gets the same budget, 12 hyper-parameter configurations per feature set, scored on validation genes. Hyper-parameters are chosen by validation MSE, and the primary model is the trained model with the best validation Spearman. These choices, the success target and the ranking weights are written to `results/frozen_choices.json` **before** any test prediction. Tree models are then refit with 5 seeds for variability, and 95% CIs come from 2,000 bootstrap resamples of test genes (paired across models).
* **Success target** (declared in the config before testing): primary metric = test Spearman over all held-out genes; target ≥ 0.20 with a 95% CI excluding 0; comparability tolerance = the primary may trail the strongest trained comparator by at most 0.02.
* **Ranking** (`t4/rank_surface.py`): eligible = HPA plasma-membrane genes in the test split. Score = 0.4·pct(enrichment) + 0.3·pct(tumour expression) + 0.1·pct(amplification) + 0.2·pct(predicted dependency). Measured dependency is printed for evaluation only and is not in the score.

## Expected outputs

| File | Contents |
|---|---|
| `data/processed/build_log.json` | sample counts (tumour/normal), genes per source, labeled genes, input SHA-256 |
| `data/splits/split_manifest.csv`, `split_counts.json` | gene → split, block; counts and surface share per split |
| `results/validation_grid.csv` | every model × feature set × configuration on validation |
| `results/frozen_choices.json` | frozen configs, primary model, target, tolerance, ranking weights |
| `results/test_metrics_by_seed.csv`, `test_summary.csv` | test MAE/MSE/Spearman, all and surface, mean ± SD over seeds, 95% CIs |
| `results/comparability_verdict.json` | primary vs strongest comparator: gap, tolerance, target met |
| `results/ablation_improvement.csv` | no-amplification and interaction deltas per model |
| `results/surface_candidates_test.csv` | ranked held-out surface candidates with evidence |
| `results/figures/*.png` | test comparison, ablation, predicted vs measured |
| `models/*.joblib` | fitted model, feature list and parameters (seed 0) |

## Hardware and run time

Runs on CPU only; no GPU needed. Peak memory is about 4 GB (loading the Xena matrix).
Measured on Palmetto, 8 cores, 32 GB: `<FILL IN from logs/t4_cp2_<jobid>.out>` (expected under 30 minutes).

## Packaging the submission

```bash
git tag cp2-submission && git push origin main --tags
python scripts/package_zip.py --names Mahr Mallozzi Patel Schoon
```

Put `report.pdf`, `slides.pptx`, slide videos and the execution `.mp4` in `deliverables/` first.

## Known limitations

* Genes are joined across sources by HGNC symbol, so genes renamed between gencode v36, cBioPortal and DepMap 26Q3 drop out.
* The label is cell-line dependency, not tumour dependency or antigen quality; surface genes are often nonessential, so low predicted dependency must not exclude a candidate.
* TCGA normal lung samples are few (about 60), so the enrichment feature is noisy.
* Ranking weights are our own fixed choice; the course public-data guide's locked score should replace them if provided.
* HPA covers about 13,000 genes; unannotated genes are treated as non-surface.
