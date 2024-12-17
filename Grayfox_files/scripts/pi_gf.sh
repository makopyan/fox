#!/bin/sh
#SBATCH --account=jazlynmo_738
#SBATCH --job-name=pi_gf
#SBATCH --output=/scratch1/marjanak/pifst/logs/pi_gf.out
#SBATCH --error=/scratch1/marjanak/pifst/logs/pi_gf.err
#SBATCH --partition=qcb 
#SBATCH --time=6:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 
#SBATCH --array=1-32

module purge

module load gcc/11.3.0 intel/19.0.4 vcftools/0.1.14

vcftools --gzvcf /scratch1/marjanak/gvcf6/gf/gf.west.AN12.gvcf.gz --window-pi 50000 --chr chr${SLURM_ARRAY_TASK_ID} --out gf_west_pi50kb_chr${SLURM_ARRAY_TASK_ID}

vcftools --gzvcf /scratch1/marjanak/gvcf6/gf/gf.east.AN12.gvcf.gz --window-pi 50000 --chr chr${SLURM_ARRAY_TASK_ID} --out gf_east_pi50kb_chr${SLURM_ARRAY_TASK_ID}