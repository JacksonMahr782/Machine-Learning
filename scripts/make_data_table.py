"""Count rows and columns in each downloaded data file and write a summary table.

Usage (from the folder that holds the downloads):
    python make_data_table.py <data_folder>

Writes data_table_counts.csv and data_table_counts.md next to this script.
Uses only the Python standard library. The E-MTAB-16585 zips are large, so
counting their lines can take several minutes.
"""
import csv
import gzip
import io
import sys
import tarfile
import zipfile
from pathlib import Path

SAMPLE_ROWS = 3
SAMPLE_COLS = 6


def summarize(name, lines, delim):
    """Return rows/columns/header/sample for an iterator of text lines."""
    reader = csv.reader(lines, delimiter=delim)
    header = next(reader, [])
    sample, n_rows = [], 0
    for row in reader:
        if n_rows < SAMPLE_ROWS:
            sample.append(row[:SAMPLE_COLS])
        n_rows += 1
    return {"file": name, "rows": n_rows, "columns": len(header),
            "header": header[:SAMPLE_COLS], "sample": sample}


def text(binary):
    return io.TextIOWrapper(binary, encoding="utf-8", errors="replace")


def delim_for(name):
    return "," if name.lower().endswith(".csv") else "\t"


def inspect(path):
    name = path.name
    low = name.lower()
    if low.endswith((".tar.gz", ".tgz")):
        with tarfile.open(path, "r:gz") as tar:
            for member in tar:
                if member.name.endswith("data_cna.txt"):
                    yield summarize(f"{name} :: data_cna.txt", text(tar.extractfile(member)), "\t")
    elif low.endswith(".zip"):
        with zipfile.ZipFile(path) as z:
            for inner in z.namelist():
                if inner.lower().endswith((".txt", ".tsv", ".csv")):
                    print(f"  reading {inner}", flush=True)
                    with z.open(inner) as f:
                        yield summarize(f"{name} :: {inner}", text(f), delim_for(inner))
    elif low.endswith(".gz"):
        with gzip.open(path, "rb") as f:
            yield summarize(name, text(f), delim_for(name[:-3]))
    elif low.endswith((".csv", ".tsv", ".txt")):
        with open(path, "rb") as f:
            yield summarize(name, text(f), delim_for(name))


def main():
    folder = Path(sys.argv[1] if len(sys.argv) > 1 else ".")
    out_dir = Path(__file__).resolve().parent
    results = []
    for path in sorted(folder.rglob("*")):
        if path.is_file():
            print(f"Checking {path.name}", flush=True)
            results.extend(inspect(path))

    with open(out_dir / "data_table_counts.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["file", "data_rows", "columns", "first_columns"])
        for r in results:
            w.writerow([r["file"], r["rows"], r["columns"], " | ".join(r["header"])])

    with open(out_dir / "data_table_counts.md", "w") as f:
        f.write("# Downloaded data: size of each table\n\n")
        f.write("| File | Data rows | Columns |\n|---|---:|---:|\n")
        for r in results:
            f.write(f"| {r['file']} | {r['rows']:,} | {r['columns']:,} |\n")
        f.write("\n## First rows of each table (first 6 columns)\n")
        for r in results:
            f.write(f"\n### {r['file']}\n\n")
            f.write("| " + " | ".join(r["header"]) + " |\n")
            f.write("|" + "---|" * len(r["header"]) + "\n")
            for row in r["sample"]:
                f.write("| " + " | ".join(c[:30] for c in row) + " |\n")
    print(f"Done: {len(results)} tables written to {out_dir / 'data_table_counts.md'}")


if __name__ == "__main__":
    main()
