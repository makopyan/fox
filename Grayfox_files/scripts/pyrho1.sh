#!/bin/sh
#SBATCH --job-name=pyrho1
#SBATCH --output=/scratch1/genchev/pyrho1.out  
#SBATCH --error=/scratch1/genchev/pyrho1.err
#SBATCH --time=04:00:00
#SBATCH --account=jazlynmo_738
#SBATCH --partition=main #partition to run the job on 
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3 #one cpu per task
#SBATCH --mem-per-cpu=5000MB #this is equivalent to 5G
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=genchev@usc.edu

module load gcc/8.3.0
module load python/3.6.8

mamba activate my-pyrho-env


# Arctic fox, East
pyrho make_table --smcpp_file split_vcf_n6/east.af.vcf.gz -n 12 -m 4.5e-9 -o arcticfox_rmap/arcticfox_east.hdf

pyrho hyperparam --samplesize 12 --tablefile arcticfox_rmap/arcticfox_east.hdf --ploidy 2 --mu 4.5e-9 --split_vcf_n6/east.af.vcf.gz -o arcticfox_rmap/arcticfox_east_hyperparam


# Arctic fox, West
pyrho make_table --smcpp_file split_vcf_n6/west.af.vcf.gz -n 12 -m 4.5e-9 -o arcticfox_rmap/arcticfox_west.hdf

pyrho hyperparam --samplesize 12 --tablefile arcticfox_rmap/arcticfox_west.hdf --ploidy 2 --mu 4.5e-9 --split_vcf_n6/west.af.vcf.gz -o arcticfox_rmap/arcticfox_west_hyperparam


# Gray fox, East
pyrho make_table --smcpp_file split_vcf_n6/east.gf.vcf.gz -n 12 -m 4.5e-9 -o grayfox_rmap/grayfox_east.hdf

pyrho hyperparam --samplesize 12 --tablefile grayfox_rmap/grayfox_east.hdf --ploidy 2 --mu 4.5e-9 --split_vcf_n6/east.gf.vcf.gz -o grayfox_rmap/grayfox_east_hyperparam


# Gray fox, West
pyrho make_table --smcpp_file split_vcf_n6/west.gf.vcf.gz -n 12 -m 4.5e-9 -o grayfox_rmap/grayfox_west.hdf

pyrho hyperparam --samplesize 12 --tablefile grayfox_rmap/grayfox_west.hdf --ploidy 2 --mu 4.5e-9 --split_vcf_n6/west.gf.vcf.gz -o grayfox_rmap/grayfox_west_hyperparam


# Canfam4, East
pyrho make_table --smcpp_file split_vcf_n6/east.cf4.vcf.gz -n 12 -m 4.5e-9 -o canfam4_rmap/canfam4_east.hdf

pyrho hyperparam --samplesize 12 --tablefile canfam4_rmap/canfam4_east.hdf --ploidy 2 --mu 4.5e-9 --split_vcf_n6/east.cf4.vcf.gz -o canfam4_rmap/canfam4_east_hyperparam


#Canfam4, West
pyrho make_table --smcpp_file split_vcf_n6/west.cf4.vcf.gz -n 12 -m 4.5e-9 -o canfam4_rmap/canfam4_west.hdf

pyrho hyperparam --samplesize 12 --tablefile canfam4_rmap/canfam4_west.hdf --ploidy 2 --mu 4.5e-9 --split_vcf_n6/west.cf4.vcf.gz -o canfam4_rmap/canfam4_west_hyperparam


# Canfam3, East
pyrho make_table --smcpp_file split_vcf_n6/east.cf3.vcf.gz -n 12 -m 4.5e-9 -o canfam3_rmap/canfam3_east.hdf

pyrho hyperparam --samplesize 12 --tablefile canfam3_rmap/canfam3_east.hdf --ploidy 2 --mu 4.5e-9 --split_vcf_n6/east.cf3.vcf.gz -o canfam3_rmap/canfam3_east_hyperparam


# Canfam3, West
pyrho make_table --smcpp_file split_vcf_n6/west.cf3.vcf.gz -n 12 -m 4.5e-9 -o canfam3_rmap/canfam3_west.hdf

pyrho hyperparam --samplesize 12 --tablefile canfam3_rmap/canfam3_west.hdf --ploidy 2 --mu 4.5e-9 --split_vcf_n6/west.cf3.vcf.gz -o canfam3_rmap/canfam3_west_hyperparam
