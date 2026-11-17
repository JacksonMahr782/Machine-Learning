"""Shared helpers: config loading, paths, feature-set definitions."""
from pathlib import Path
import json, hashlib
import yaml

ROOT = Path(__file__).resolve().parents[1]

BASE_FEATURES = ["expr_mean_tumor", "expr_enrichment", "expr_sd_tumor",
                 "gistic_amp_freq", "gistic_highamp_freq"]
INTERACTIONS = ["mean_x_amp", "enrich_x_amp"]
AMP_FEATURES = ["gistic_amp_freq", "gistic_highamp_freq"] + INTERACTIONS

FEATURE_SETS = {
    "base": BASE_FEATURES,                                   # CP1 plan, A1
    "interactions": BASE_FEATURES + INTERACTIONS,            # improvement, A6
    "no_amp": [f for f in BASE_FEATURES if f not in AMP_FEATURES],  # ablation, A6
}


def load_config(path):
    with open(path) as fh:
        cfg = yaml.safe_load(fh)
    for k in ("raw_dir", "processed_dir", "splits_dir", "results_dir", "models_dir"):
        cfg[k] = (ROOT / cfg[k]).resolve()
        if k != "raw_dir":
            cfg[k].mkdir(parents=True, exist_ok=True)
    return cfg


def raw(cfg, key):
    return cfg["raw_dir"] / cfg["files"][key]


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_json(obj, path):
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=2, default=str)
