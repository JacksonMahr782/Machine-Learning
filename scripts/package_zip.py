"""Build the Canvas ZIP from the committed revision plus deliverables.

Put report.pdf, slides.pptx (or slides.pdf), any slide videos, and the execution
MP4 in deliverables/ first, commit and tag, then run from the repo root:

    python scripts/package_zip.py --names Mahr Mallozzi Patel Schoon

The code inside the ZIP comes from `git archive HEAD`, so it matches the submitted
commit exactly. Trained checkpoints (models/*.joblib, git-ignored) are added separately.
"""
import argparse, subprocess, zipfile, io
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REQUIRED = ["report.pdf", "README.md"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--names", nargs="+", required=True, help="member last names")
    ap.add_argument("--rev", default="HEAD")
    a = ap.parse_args()
    out = ROOT / f"T4_Checkpoint2_{'_'.join(a.names)}.zip"
    commit = subprocess.check_output(["git", "rev-parse", a.rev], cwd=ROOT, text=True).strip()
    if subprocess.check_output(["git", "status", "--porcelain", "--untracked-files=no"], cwd=ROOT, text=True).strip():
        print("WARNING: uncommitted changes; the ZIP uses the committed version only.")

    deliv = ROOT / "deliverables"
    have = {p.name for p in deliv.glob("*")} if deliv.exists() else set()
    missing = [f for f in ["report.pdf"] if f not in have]
    if not any(n.startswith("slides.") for n in have):
        missing.append("slides.pptx or slides.pdf")
    if not any(n.endswith(".mp4") for n in have):
        missing.append("execution video (.mp4)")
    if missing:
        raise SystemExit(f"missing in deliverables/: {missing}")

    code = subprocess.check_output(["git", "archive", "--format=zip", commit], cwd=ROOT)
    with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as z:
        with zipfile.ZipFile(io.BytesIO(code)) as src:
            for info in src.infolist():
                if not info.is_dir():
                    z.writestr(f"code/{info.filename}", src.read(info))
        for p in sorted(deliv.glob("*")):
            z.write(p, p.name)                       # report.pdf, slides, videos at top level
        z.write(ROOT / "README.md", "README.md")
        for p in sorted((ROOT / "models").glob("*.joblib")):
            z.write(p, f"code/models/{p.name}")
        z.writestr("COMMIT.txt", f"https://github.com/JacksonMahr782/Machine-Learning\ncommit {commit}\n")
    print(f"wrote {out.name} ({out.stat().st_size / 1e6:.1f} MB) from commit {commit[:10]}")


if __name__ == "__main__":
    main()
