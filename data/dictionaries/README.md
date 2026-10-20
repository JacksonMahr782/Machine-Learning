# T4 Checkpoint 1 - Data start evidence

Owner: Jackson Mahr. Status: all 6 public sources accessed and downloaded (Oct 15-19, 2026).

| File | What it is |
|---|---|
| `download_log.csv` | Every source: URL, files, release/version, date, size, status, proof screenshot |
| `source_dictionary.csv` | What each source contains and how the sources join together |
| `gene_label_dictionary.csv` | Planned columns of the one-row-per-gene table: features, label, excluded columns |
| `split_dictionary.csv` | Gene-level 70/15/15 split rules, stratified by surface status |
| `proof/` | Screenshots of source pages, downloads list and downloaded files |

## Key rules
- Leakage: nothing derived from DepMap gene effects is used as a feature. DepMap only supplies the label.
- Scope: the TROP2 paper reproduction (E-MTAB-16585) and the LUAD dependency model are separate tasks.

## Where things live
- GitHub (`data/dictionaries/`): these CSVs and code. Small text files only.
- Box (`data/raw/`): the downloaded data files (several GB, too large for GitHub).

## Next action (owner: Jackson)
Open the three E-MTAB-16585 zips, confirm all 12 processed .txt files are present, and load them with the authors' Figure_01.ipynb before Checkpoint 2.
