#!/bin/bash
#SBATCH --job-name=qom_fstm
#SBATCH --partition=production
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/slurm_fstm_%j.log
#SBATCH --error=logs/slurm_fstm_%j.err

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
mkdir -p logs results
PY=/home/nanaengo/miniforge3/envs/qom/bin/python3
exec "$PY" scripts/full_spectrum_transition_map.py \
  --output results/full_spectrum_transition_map.json \
  --drives 0.2,0.5,1.0,2.0,4.0,8.0 \
  --points 9 \
  --n-steps 1000000 \
  --transient-steps 100000 \
  --workers 6 \
  --seed 20260815
