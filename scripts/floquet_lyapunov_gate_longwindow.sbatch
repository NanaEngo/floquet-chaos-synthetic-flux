#!/bin/bash
#SBATCH --job-name=qom_floquet_longwin
#SBATCH --partition=production
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=logs/slurm_floquet_longwin_%j.log
#SBATCH --error=logs/slurm_floquet_longwin_%j.err

# Finding K follow-up: the default 1e6-step gate run gives
# max Floquet-vs-QR difference 5.45e-4 > 5e-4 (FAIL), while the lost 2026-08-13
# pilot artifact reportedly passed at 2.68e-4 with identical code (hashes match)
# and identical orbit residual (8.485e-9). Test whether a longer QR window
# converges to the documented PASS. Writes to a scratch dir; promoted to
# results/floquet_lyapunov.json only after the numbers are checked.

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
mkdir -p logs results/pilot_longwindow

export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=4
export MKL_NUM_THREADS=4

PY=/home/nanaengo/miniforge3/envs/qom/bin/python3
"$PY" scripts/floquet_lyapunov_gate.py \
  --output-dir results/pilot_longwindow \
  --n-steps 5000000 \
  --transient-steps 500000
