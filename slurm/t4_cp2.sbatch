#!/bin/bash
# Palmetto (Slurm) job for the dependency pipeline. Submit from the repo root:
#   sbatch slurm/t4_cp2.sbatch
# Edit the partition/module lines to match your Palmetto account if needed.
#SBATCH --job-name=t4_cp2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=logs/t4_cp2_%j.out

set -euo pipefail
mkdir -p logs
module load anaconda3 2>/dev/null || module load miniforge3 2>/dev/null || true
source activate t4cp2 2>/dev/null || conda activate t4cp2

echo "host: $(hostname)  cpus: $SLURM_CPUS_PER_TASK  start: $(date)"
echo "commit: $(git rev-parse HEAD)"
python -c "import sklearn, pandas, numpy; print('sklearn', sklearn.__version__, 'pandas', pandas.__version__, 'numpy', numpy.__version__)"
/usr/bin/time -v bash run_all.sh configs/t4_cp2.yaml
echo "end: $(date)"
