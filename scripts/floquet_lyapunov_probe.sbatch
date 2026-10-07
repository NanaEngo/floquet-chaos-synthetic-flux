#!/bin/bash
#SBATCH --job-name=qom_flq_probe
#SBATCH --partition=production
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/slurm_flq_probe_%j.log
#SBATCH --error=logs/slurm_flq_probe_%j.err

# Finding K follow-up: the default-window floquet_lyapunov_gate run
# (n_steps 1e6) FAILs its own 5e-4 rate/spectrum criterion (5.45e-4) while
# the lost 2026-08-13 artifact passed with 2.68e-4. Probe longer QR windows
# (usage: sbatch scripts/floquet_lyapunov_probe.sbatch <n_steps> [transient])
# into a scratch subdir; the best PASS is adopted into results/ afterwards.

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
mkdir -p logs "results/floquet_probe_${1}_${2:-100000}"

export OMP_NUM_THREADS=2
export OPENBLAS_NUM_THREADS=2
export MKL_NUM_THREADS=2

PY=/home/nanaengo/miniforge3/envs/qom/bin/python3
exec "$PY" scripts/floquet_lyapunov_gate.py \
  --output-dir "results/floquet_probe_${1}_${2:-100000}" \
  --n-steps "${1}" \
  --transient-steps "${2:-100000}"
