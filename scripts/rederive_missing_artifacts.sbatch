#!/bin/bash
#SBATCH --job-name=qom_rederive_K
#SBATCH --partition=production
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/slurm_rederive_K_%j.log
#SBATCH --error=logs/slurm_rederive_K_%j.err

# Finding K (report 2026-10-05): regenerate artifacts lost from results/ that
# back manuscript claims: feasibility frontier + 3-seed onset, optimized
# readout reference, pilot Floquet-Lyapunov checkpoint. Stored producer
# defaults only; each step's exit code is recorded so a late failure does not
# hide earlier successes.

set -uo pipefail
cd "${SLURM_SUBMIT_DIR}"
mkdir -p logs results

export OMP_NUM_THREADS=8
export OPENBLAS_NUM_THREADS=8
export MKL_NUM_THREADS=8

PY=/home/nanaengo/miniforge3/envs/qom/bin/python3
rc_all=0

echo "== [1/4] feasibility_frontier (5 drives x 8 g1, 6 workers) =="
"$PY" scripts/feasibility_frontier.py \
  --output results/feasibility_frontier.json \
  --workers 6 || rc_all=$?

echo "== [2/4] feasibility_frontier_onset (3 seeds x 8 g1, sequential) =="
"$PY" scripts/feasibility_frontier_onset.py \
  --output results/feasibility_frontier_onset.json || rc_all=$?

echo "== [3/4] optimized_linear_sensor (360 angles, 5 thetas) =="
"$PY" scripts/optimized_linear_sensor.py \
  --output results/optimized_linear_sensor.json || rc_all=$?

echo "== [4/4] floquet_lyapunov_gate (pilot checkpoint, 1e6 steps) =="
"$PY" scripts/floquet_lyapunov_gate.py \
  --output-dir results || rc_all=$?

echo "== done, last-failure rc=$rc_all =="
exit "$rc_all"
