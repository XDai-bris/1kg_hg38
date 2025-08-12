#!/bin/bash
# SLURM submission script for BC4 (call: sbatch qsub.sh)

#SBATCH --job-name=run_pipeline_arr
#SBATCH --partition=veryshort
#SBATCH --array=0-114%84          # 115 tasks total, max 56 running at once (~2 nodes worth)
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G                   # per task memory; bump to 6G only if needed
#SBATCH --time=04:30:00            # under 6h improves backfill
#SBATCH --account=smed001801
#SBATCH --output=/user/home/xd14188/repo/1kg_hg38/logs/arr-%A_%a.out
#SBATCH --error=/user/home/xd14188/repo/1kg_hg38/logs/arr-%A_%a.err

set -euo pipefail

# Jump to repo directory (adjust if needed)
cd /user/home/xd14188/repo/1kg_hg38

# Run the pipeline (modules/threads are handled inside)
bash ./run_chr_pop_pipeline_imputedVarID.sh
