#!/bin/bash
#SBATCH --job-name=qom_fisher_uq60
#SBATCH --partition=production
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --output=logs/slurm_fisher_uq60_%j.log
#SBATCH --error=logs/slurm_fisher_uq60_%j.err

set -euo pipefail
cd /home/nanaengo/Floquet_COM_Synthetic_Magnetism_Analysis/Q1_Papers/manuscripts/track_D_chaos
mkdir -p logs results

export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4
export MKL_NUM_THREADS=4

PY=/home/nanaengo/miniforge3/envs/qom/bin/python3
exec "$PY" scripts/matched_fisher_uq.py \
  --output results/matched_fisher_uq.json \
  --n-samples 60 \
  --seed 42 \
  --theta 3.141592653589793 \
  --hopping-range 1e-3 1e-1 \
  --n-th-range 0.01 1.0 \
  --drive-range 0.18 0.22 \
  --detector-variance 0.01
