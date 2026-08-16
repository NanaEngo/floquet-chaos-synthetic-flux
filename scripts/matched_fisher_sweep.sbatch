#!/bin/bash
#SBATCH --job-name=qom_fisher_sweep
#SBATCH --partition=production
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=logs/slurm_fisher_sweep_%j.log
#SBATCH --error=logs/slurm_fisher_sweep_%j.err

set -euo pipefail
cd /home/nanaengo/Floquet_COM_Synthetic_Magnetism_Analysis/Q1_Papers/manuscripts/track_D_chaos
mkdir -p logs results

PY=/home/nanaengo/miniforge3/envs/qom/bin/python3
exec "$PY" scripts/matched_fisher_sweep.py \
  --output results/matched_fisher_sweep.json \
  --drives 0.1,0.2,0.5,1.0 \
  --n-ths 0.0,0.1,0.5,1.0
