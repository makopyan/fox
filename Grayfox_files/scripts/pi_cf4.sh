#!/bin/sh
#SBATCH --account=jazlynmo_738
#SBATCH --job-name=pi_cf4
#SBATCH --output=/scratch1/marjanak/pi/logs/pi_cf4.out
#SBATCH --error=/scratch1/marjanak/pi/logs/pi_cf4.err
#SBATCH --time=6:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 
#SBATCH --array=1-38

module purge

module load gcc/11.3.0 intel/19.0.4 vcftools/0.1.14

vcftools --gzvcf /scratch1/marjanak/gvcf6/cf4/cf4.west.AN12.gvcf.gz --window-pi 50000 --chr chr${SLURM_ARRAY_TASK_ID} --out cf4/cf4_west_pi50kb_chr${SLURM_ARRAY_TASK_ID}

vcftools --gzvcf /scratch1/marjanak/gvcf6/cf4/cf4.east.AN12.gvcf.gz --window-pi 50000 --chr chr${SLURM_ARRAY_TASK_ID} --out cf4/cf4_east_pi50kb_chr${SLURM_ARRAY_TASK_ID}