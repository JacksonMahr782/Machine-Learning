@echo off
REM Windows version of run_all.sh. Run from the repo root: run_all.bat
set CONFIG=%1
if "%CONFIG%"=="" set CONFIG=configs/t4_cp2.yaml
for %%S in (build_gene_table make_splits train_eval rank_surface make_figures) do (
  echo === %%S ===
  python -m t4.%%S --config %CONFIG% || exit /b 1
)
