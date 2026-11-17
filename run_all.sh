#!/usr/bin/env bash
# Full T4 Checkpoint 2 dependency pipeline. Run from the repo root.
set -euo pipefail
CONFIG=${1:-configs/t4_cp2.yaml}
python -m t4.build_gene_table --config "$CONFIG"   # data/processed/gene_table.csv, build_log.json
python -m t4.make_splits      --config "$CONFIG"   # data/splits/split_manifest.csv, split_counts.json
python -m t4.train_eval       --config "$CONFIG"   # results/*.csv|json, models/*.joblib
python -m t4.rank_surface     --config "$CONFIG"   # results/surface_candidates_test.csv
python -m t4.make_figures     --config "$CONFIG"   # results/figures/*.png
