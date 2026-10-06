#!/bin/bash
# MPS run on one GPU (config "backend": "cupy" or "torch"). Fill in the placeholders:
#   <PARTITION> <ACCOUNT> <ENV> (environment activation command) <CONFIG> (configs/*.json)
#SBATCH --job-name=su2torus-mps-gpu
#SBATCH --partition=<PARTITION>
#SBATCH --account=<ACCOUNT>
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-task=1
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x-%j.out
set -euo pipefail
<ENV>
cd "$SLURM_SUBMIT_DIR"
python3 scripts/run_classical.py --config <CONFIG>
