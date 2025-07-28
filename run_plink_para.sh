#!/bin/bash

#SBATCH --job-name=download_job
#SBATCH --partition=mrcieu,cpu,test
#SBATCH --nodes=5
#SBATCH --ntasks-per-node=1
#SBATCH --time=2:0:0
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
# bash plink_parallel.sh 
# or for the final script
# chmod +x plink_paralel.sh
# or
# bash plink_para_2_bcftools-R.sh
# or 
# bash merge_vcf.sh
bash plink_para_2_bcftools-R.sh