#!/bin/bash

#SBATCH --job-name=download_job
#SBATCH --partition=mrcieu,cpu,test
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=1
#SBATCH --time=48:0:0
#SBATCH --mem=24G
#SBATCH --account=smed001801

# Load required modules
module load bcftools/1.19-openblas-2333
module load plink2/2.00a4.3-netlib-lapack-qbk6
# Optional: check module environment
which bcftools
bcftools --version
which parallel
which plink2
plink2 --version
# Run
chmod +x /user/home/xd14188/repo/1kg_hg38/scripts/process_1kg_pipeline.sh
bash /user/home/xd14188/repo/1kg_hg38/scripts/process_1kg_pipeline.sh
# ============================================================================