#!/bin/sh
#SBATCH --job-name=phyloP_lift_cf4
#SBATCH --output=/scratch1/marjanak/phyloP/phyloP_lift_cf4.out
#SBATCH --error=/scratch1/marjanak/phyloP/phyloP_lift_cf4.err
#SBATCH --partition=qcb
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=20G




module purge

eval "$(conda shell.bash hook)"

conda activate /home1/marjanak/.conda/envs/ucsctools

bigWigAverageOverBed zoonomia-cf3.1-lifted-to-cf4.liftover.phylop.20210708.bw cf4_gf_fst_outlier.bed cf4_outlier_phylop.out