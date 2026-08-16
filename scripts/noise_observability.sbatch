#!/bin/bash
#SBATCH --job-name=qom_noise_obs
#SBATCH --partition=production
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/slurm_noise_obs_%j.log
#SBATCH --error=logs/slurm_noise_obs_%j.err

set -euo pipefail
cd /home/nanaengo/Floquet_COM_Synthetic_Magnetism_Analysis/Q1_Papers/manuscripts/track_D_chaos
mkdir -p logs results

export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4
export MKL_NUM_THREADS=4

PY=/home/nanaengo/miniforge3/envs/qom/bin/python3
exec "$PY" scripts/noise_observability.py \
  --output results/noise_observability.json \
  --n-ens 32 \
  --n-steps 50000 \
  --transient-steps 5000 \
  --dt 0.01 \
  --detector-variance 0.01 \
  --n-th 0.1 0.1 \
  --seed 42
