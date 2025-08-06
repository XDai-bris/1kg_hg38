#!/bin/bash

#SBATCH --job-name=run_pipeline
#SBATCH --partition=mrcieu,cpu
#SBATCH --nodes=12
#SBATCH --ntasks-per-node=1
#SBATCH --time=26:0:0
#SBATCH --mem=60G
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



bash run_chr_pop_pipeline_imputedVarID.sh
bash mergeChrBed.sh
