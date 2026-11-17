# CP2 report: Data section (draft, owner Veer)

Numbers are from the real run on 2026-11-17 (`data/processed/build_log.json`, `data/splits/split_counts.json`).

---

## Data, features and splits

**Sources.** We use four public sources, all downloaded Oct 15–19, 2026 (`data/dictionaries/download_log.csv`). DepMap Public 26Q3 (`CRISPRGeneEffect.csv`, `Model.csv`) supplies the label only. The UCSC Xena GDC TCGA-LUAD STAR-TPM matrix (version 05-10-2024, GDC DR40, log2(TPM+1), gencode v36) and cBioPortal `luad_tcga_pan_can_atlas_2018` GISTIC calls (`data_cna.txt`) supply the features. HPA v25 `subcellular_location.tsv` decides surface eligibility only. Genes are joined by HGNC symbol.

**Label.** For gene *g*, y_g = −median of Chronos effect over the LUAD lines (OncotreeCode = LUAD) where *g* was measured. Higher values mean stronger dependency. Genes measured in fewer than 3 lines stay unlabeled and are not in any split. Model.csv lists 98 LUAD models, of which 53 have CRISPR effects. Each labeled gene was measured in 18–53 lines (median 53). 16,756 genes in the final table are labeled. No predictor is derived from effect data.

**Features.** All features are computed per gene from 530 tumour samples (516 patients) and 59 normal samples (59 patients) in Xena, and from 511 tumour samples in GISTIC. With x = log2(TPM+1):

| Feature | Definition |
|---|---|
| `expr_mean_tumor` | mean of x over tumour samples |
| `expr_enrichment` | mean_tumor(x) − mean_normal(x) |
| `expr_sd_tumor` | SD of x over tumour samples |
| `gistic_amp_freq` | fraction of samples with GISTIC call ≥ +1 |
| `gistic_highamp_freq` | fraction of samples with GISTIC call = +2 |
| `mean_x_amp`, `enrich_x_amp` | `expr_mean_tumor` × `gistic_amp_freq`, `expr_enrichment` × `gistic_amp_freq` (improvement set only) |

**Preprocessing.** Xena values are already log2(TPM+1), with aliquots of the same sample averaged by Xena. Samples are labelled tumour (TCGA sample code 01–09) or normal (10–19). Ensembl IDs are mapped to symbols with the gencode v36 probemap, and when several IDs share a symbol we keep the most expressed one. PAR_Y duplicates are removed. GISTIC frequencies use each gene's non-missing calls as the denominator. Genes missing any feature are dropped (0 genes on the real data). The table has 22,646 genes with all features. None of these steps learns anything from the data. The only learned preprocessing is Ridge's standardisation, which is fit on training genes. We do no data-driven feature selection. The only feature choice is among three fixed sets (base, +interactions, −amplification), and it is made on validation genes.

**Surface rule.** A gene is a surface candidate if its HPA main or additional location includes "Plasma membrane" (the CP1 definition; 2,263 of 13,597 HPA v25 genes). Cell Junctions alone does not qualify, and all reliability levels are kept. Genes not in HPA count as non-surface. The flag is used only to stratify the split and to choose which genes get ranked. It is never a feature or the label.

**Split.** The unit is the gene. Neighbouring genes share amplicons, which gives them near-identical GISTIC features, so genes in the same 10 Mb genomic block always go to the same split (this is a change from CP1). We used StratifiedGroupKFold (20 folds, seed 42, stratified on `is_surface`) and assigned 14 folds to train, 3 to validation and 3 to test. The script checks that no block crosses splits and that the test set has at least 20 surface genes.

| Split | Genes | Surface | Blocks | Mean y |
|---|---|---|---|---|
| Train | 11,749 (70.1%) | 1,435 (12.2%) | 212 | 0.148 |
| Val | 2,504 (14.9%) | 302 (12.1%) | 48 | 0.130 |
| Test | 2,503 (14.9%) | 304 (12.1%) | 46 | 0.155 |
| Total | 16,756 | 2,041 | 306 | 0.146 |

**Leakage controls.** DepMap supplies only y and HPA supplies only eligibility. The Ridge scaler is fit inside the pipeline on training genes. Hyperparameters, the feature set and the primary-model rule are tuned on validation genes and written to `frozen_choices.json` before any test prediction. The test set is scored once.

**Limitations.** An amplicon that crosses a block boundary can still be split across sets, and paralogs are not grouped. Renamed genes drop out of the symbol join. With only 59 normal samples, enrichment is noisy. 4,783 labeled genes have no HPA annotation and count as non-surface. The label is right-skewed (median 0.03, max 5.3 for core essentials such as RPL15), which inflates MSE/MAE more than Spearman. The label measures cell-line dependency, not tumour dependency or antigen quality.

---

## Slides (2-minute slot: dataset, split units and counts, preprocessing, features)

**Slide 1: Data and label (~40 s)**
- A flow diagram of the four sources and what each supplies: DepMap → y, Xena + GISTIC → features, HPA → surface flag.
- y = −median Chronos over ≥ 3 LUAD lines. Show 53 lines and 16,756 labeled genes.
- Line to say: "Effect data never enters a feature, and HPA never enters the model."

**Slide 2: Features (~40 s)**
- The five base features plus two interactions, each with a one-line formula (use the table above).
- Sample counts: 530 tumour / 59 normal / 511 GISTIC.
- A histogram of y with the surface genes overlaid (optional; take it from `make_figures.py` if it exists).

**Slide 3: Split and leakage (~40 s)**
- A cartoon of chromosome blocks: whole 10 Mb blocks are coloured train, val or test.
- The count table (train/val/test × genes/surface).
- Checklist: no block crosses splits ✓, ≥ 20 test surface ✓, scaler fit on train ✓, choices frozen before test ✓.
