#!/bin/bash

#SBATCH --job-name=download_job
#SBATCH --partition=mrcieu,cpu,test
#SBATCH --nodes=6
#SBATCH --ntasks-per-node=6
#SBATCH --time=1:0:0
#SBATCH --mem=24G
#SBATCH --account=smed001801

# Load required modules
module load bcftools/1.19-openblas-2333
module load plink2/2.00a4.3-netlib-lapack-qbk6
module load languages/R/4.4.3 
module load plink/1.9-beta6.27-openblas-ibxp
# Optional: check module environment
which bcftools
bcftools --version
which parallel
which plink2
plink2 --version
which plink
plink --version
R --version



# bash run_plink_para.sh
# bash mergeChrBed.sh
# bash mergeSolv.sh
 bash run_chrX_pipeline.sh