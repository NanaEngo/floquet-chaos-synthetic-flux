#!/bin/bash
#SBATCH --job-name=qom_chaos_fisher
#SBATCH --partition=production
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=logs/slurm_chaos_fisher_%j.log
#SBATCH --error=logs/slurm_chaos_fisher_%j.err

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
mkdir -p logs results

export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4
export MKL_NUM_THREADS=4

PY=/home/nanaengo/miniforge3/envs/qom/bin/python3
exec "$PY" scripts/chaotic_attractor_fisher.py \
  --output results/chaotic_attractor_fisher.json \
  --n-ens 1024 \
  --n-steps 120000 \
  --transient-steps 3000 \
  --seed 42
