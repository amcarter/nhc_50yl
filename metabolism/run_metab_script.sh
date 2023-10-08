#!/bin/bash -l
#SBATCH --job-name NHC_metab
#SBATCH --time=48:00:00
#SBATCH --mem=30G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=alice.carter@flbs.umt.edu
#SBATCH --account=modelscape
#SBATCH --output=slurmlogs/SV_run.out

echo "SLURM_JOB_ID:" $SLURM_JOB_ID
echo "SLURM_JOB_NAME:" $SLURM_JOB_NAME

module load gcc/12.2.0 arcc/1.0 r/4.2.2

cd /project/modelscape/users/acarte26/NHC_metab

Rscript --vanilla run_NHC_streamMetabolizer_raymond_K600.R

