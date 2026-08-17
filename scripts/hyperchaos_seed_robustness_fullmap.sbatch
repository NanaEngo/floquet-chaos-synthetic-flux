#!/bin/bash
#SBATCH --job-name=qom_seed_full
#SBATCH --partition=production
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=logs/slurm_seed_fullmap_%j.log
#SBATCH --error=logs/slurm_seed_fullmap_%j.err

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
mkdir -p logs results

export OMP_NUM_THREADS=8
export OPENBLAS_NUM_THREADS=8
export MKL_NUM_THREADS=8

PY=/home/nanaengo/miniforge3/envs/qom/bin/python3
exec "$PY" scripts/hyperchaos_seed_robustness.py \
  --output results/hyperchaos_seed_robustness_fullmap.json \
  --drives 0.2,1.0,2.0,4.0,8.0 \
  --points 12 \
  --n-seeds 3 \
  --g1 0.3 --g2 0.27 \
  --n-steps 1000000 \
  --transient-steps 100000 \
  --workers 6 \
  --seed-base 20260817
