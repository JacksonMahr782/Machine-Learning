# Checkpoint 2 evaluation handoff

`scripts/zachary_evaluate.py` analyzes predictions made by the team's separate training pipeline. It does not train, choose a model, validate the split, verify tuning budgets, or calculate the course ranking formula. Freeze choices before opening test predictions. Keep the frozen record under version control with a timestamp or commit predating test evaluation. A Boolean in the file is an attestation, not proof of timing.

The README proposes overall test Spearman as primary metric, a target of at least 0.20 with a 95% confidence interval excluding zero, and a tolerance of 0.02. The tolerance requires a scientific justification agreed before testing; the README does not provide one. The proposed primary selection rule is best validation Spearman after each model's hyperparameters are selected on validation MSE. Document the rule and selection result before testing, and explain the choice in the report. Do not use this script to choose a primary model based on test results.

Example decision record (replace placeholders before testing; do not treat this as an approved tolerance justification):

```json
{
  "frozen_before_test": true,
  "primary_model": "ridge",
  "primary_metric": "overall_spearman",
  "target": 0.20,
  "tolerance": 0.02,
  "tolerance_justification": "REPLACE with agreed scientific justification",
  "ranking_formula_source": "REPLACE with exact guide and formula location"
}
```

Predictions CSV: one row per held-out gene, columns `gene`, `split` (test), `block`, `is_surface` (0/1), `dependency`, `pred_constant`, `pred_ridge`, `pred_random_forest`, `pred_gradient_boosting`, `pred_primary_no_amplification`. Include additional evidence columns for the candidate table. The primary ablation must be retrained without both GISTIC columns and all amplification interactions, using the same split and validation-only selection. Verify predictions correspond to the frozen model and same test manifest; this script cannot independently establish provenance or ensure 2,000 labeled genes across all splits.

```bash
python scripts/zachary_evaluate.py --predictions results/test_predictions.csv --frozen results/frozen_choices.json --out results/zachary
```

Outputs: overall and surface MAE and Spearman with 95% intervals, paired model and ablation differences, and the target/comparability verdict. Positive differences favor the primary. Percentile intervals use 2,000 paired resamples of genomic blocks, reflecting the split's 10 Mb grouping, not independent genes; they do not represent retraining variability. Surface intervals resample blocks represented in that subset. Undefined Spearman values are omitted from bootstrap quantiles; fewer than 100 valid replicates makes the interval undefined. The constant's Spearman is undefined. Interpret small block counts and correlated genes cautiously. Report repeated-run variability separately when available. Model comparability also requires checking tuning budgets, features, split, and metrics; the JSON verdict covers only the numerical gap.

The exact course prioritization formula is not in this repository. Do not use the provisional README weights as the required score. After confirming and freezing the formula and eligibility rules, independently compute a CSV with `gene,score` for **all** eligible test genes and run with `--scores path/to/scores.csv`. The script joins and sorts those scores, retaining evidence columns, but cannot verify their formula or whether they improperly use measured dependency. Confirm that the score uses predicted dependency as supporting evidence and that nonessential antigens remain eligible. Without the score file, no ranking is generated. Record formula provenance and code alongside the results.

Check outputs before including them in the report: split counts and >=2,000 labeled genes; >=20 held-out surface genes; common test genes; target, primary gap and justified tolerance; paired intervals; ablation; and errors or failures. Neither predicted dependency nor rank establishes antigen validity or therapeutic efficacy.
