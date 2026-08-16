#!/bin/bash
#SBATCH --job-name=qom_basin_v2
#SBATCH --partition=production
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=logs/slurm_basin_v2_%j.log
#SBATCH --error=logs/slurm_basin_v2_%j.err

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
mkdir -p logs results
PY=/home/nanaengo/miniforge3/bin/python3
exec "$PY" scripts/fast_basin_capture.py \
  --screen results/physical_parameter_reduced_bad_cavity_screen_robust_v4.json \
  --output results/fast_basin_capture_validation_v2.json \
  --samples 5 10 15 \
  --initial-condition-replicates 5 \
  --max-periods 5000 \
  --substeps 320 \
  --capture-tol 1e-5 \
  --residual-tol 1e-7 \
  --residual-check-interval 25
