#!/bin/bash
#SBATCH --job-name=qom_hmap
#SBATCH --partition=production
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/slurm_hmap_%j.log
#SBATCH --error=logs/slurm_hmap_%j.err

set -euo pipefail
cd "${SLURM_SUBMIT_DIR}"
mkdir -p logs results
PY=/home/nanaengo/miniforge3/envs/qom/bin/python3
exec "$PY" scripts/full_spectrum_hyperchaos_map.py \
  --output results/full_spectrum_hyperchaos_map_g03.json \
  --drives 0.2,1.0,2.0,4.0,8.0 \
  --points 12 \
  --g1 0.3 --g2 0.27 \
  --n-steps 1000000 \
  --transient-steps 100000 \
  --workers 6 \
  --seed 20260815
