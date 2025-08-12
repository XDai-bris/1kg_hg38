#!/bin/bash
# SLURM launcher for multithreaded population merges (5 pops × 8 CPUs each)

#SBATCH --job-name=merge_1kg_8cpu
#SBATCH --partition=veryshort        # use 'cpu' if your merge might exceed 6h
#SBATCH --ntasks=5                   # one task per population
#SBATCH --cpus-per-task=8            # plink2 --threads 8
#SBATCH --mem-per-cpu=4G             # 8 CPUs → 32G per task
#SBATCH --time=04:00:00
#SBATCH --account=smed001801
#SBATCH --output=/user/home/xd14188/repo/1kg_hg38/merged_genome/slurm-%j.out
#SBATCH --error=/user/home/xd14188/repo/1kg_hg38/merged_genome/slurm-%j.err

set -euo pipefail
cd /user/home/xd14188/repo/1kg_hg38

# mergeChrBed.sh should launch one srun per population and use --threads ${SLURM_CPUS_PER_TASK}
bash ./mergeChrBed.sh