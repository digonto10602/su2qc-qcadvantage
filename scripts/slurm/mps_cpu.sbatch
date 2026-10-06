#!/bin/bash
# MPS run on CPU nodes. Fill in the placeholders before submitting:
#   <PARTITION> <ACCOUNT> <ENV> (conda/venv activation command) <CONFIG> (configs/*.json)
#SBATCH --job-name=su2torus-mps
#SBATCH --partition=<PARTITION>
#SBATCH --account=<ACCOUNT>
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=200G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x-%j.out
set -euo pipefail
<ENV>
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
cd "$SLURM_SUBMIT_DIR"
python3 scripts/run_classical.py --config <CONFIG>
