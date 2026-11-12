#!/bin/sh
#SBATCH --job-name=spatial-mediation
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1

# Submit from the repository root: sbatch Simulations/submit.sh
set -eu
srun Rscript Simulations/Spatial_Mediation_MLE_Gaussian_Primary.R "$@"
