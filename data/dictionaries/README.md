# T4 Checkpoint 1 - Data start evidence

Group T4, Tumor surface target prioritization (CPSC 4430). Data start owner: Jackson Mahr.

**GitHub repository (code + dictionaries):** https://github.com/JacksonMahr782/Machine-Learning

Status: all 6 public sources accessed and downloaded (Oct 15-19, 2026). All 12 processed E-MTAB-16585 .txt files confirmed present (Oct 21).

| File | What it is |
|---|---|
| `DATA_TABLE.pdf` | Table of every downloaded file (source, version, size, date, what a row/column is) and the layout of our gene table |
| `METHODS_AND_EQUATIONS.pdf` | The equation behind every feature, the label, each model (constant, Ridge, random forest, gradient boosting), the loss, the metrics (MAE, Spearman), the improvement/ablation, and the paper's statistical tests |
| `download_log.csv` | Every source: URL, files, release/version, date, size, status, proof screenshot |
| `data_files_table.csv` | Table 1 of DATA_TABLE.pdf as a CSV |
| `gene_table_example.csv` | Column layout of the one-row-per-gene table (example values, not real data) |
| `source_dictionary.csv` | What each source contains and how the sources join together |
| `gene_label_dictionary.csv` | Planned columns of the gene table: features, label, excluded columns |
| `split_dictionary.csv` | Gene-level 70/15/15 split rules, stratified by surface status |
| `make_data_table.py` | Script that counts rows/columns and prints the first rows of each downloaded file |
| `T4_CP1_data_start.xlsx` | All of the above tables in one workbook |
| `proof/` | Screenshots of source pages, downloads list and downloaded files |

## Key rules
- Leakage: nothing derived from DepMap gene effects is used as a feature. DepMap only supplies the label.
- Scope: the TROP2 paper reproduction (E-MTAB-16585) and the LUAD dependency model are separate tasks.

## CP2 changes to these dictionaries
The CSVs now use the column names the code writes (`t4/build_gene_table.py`, `t4/make_splits.py`). `T4_CP1_data_start.xlsx` is kept as submitted at CP1.
- `gene_symbol` -> `gene`, `label_dependency` -> `dependency`, split label `validation` -> `val`, `splits.csv` -> `data/splits/split_manifest.csv` + `split_counts.json`.
- New feature columns: `mean_x_amp`, `enrich_x_amp` (interaction improvement). New QC/grouping columns: `chrom`, `start`, `block`, `hpa_annotated`, `hpa_reliability`.
- Split unit is still the gene, but genes are grouped by 10 Mb genomic block so a shared amplicon cannot leak across splits.
- Expression CV, mentioned in the CP1 report, is not used; SD is the variability feature.
- Surface rule is set in `configs/t4_cp2.yaml` (`surface:`). Releases recorded: Xena star_tpm 05-10-2024, HPA v25.

## Where things live
- GitHub (`data/dictionaries/`, `docs/`): CSVs, equations and code. Small text files only.
- Box (`data/raw/`): the downloaded data files (several GB, too large for GitHub).

## Next action (owner: Jackson)
Load the 12 E-MTAB-16585 files with the authors' Figure_01.ipynb on Palmetto, and run `make_data_table.py` on the raw files to record exact row and column counts, before Checkpoint 2.
