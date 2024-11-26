#!/bin/sh
#SBATCH --job-name=msmc2_size_e
#SBATCH --output=/scratch1/marjanak/gvcf6/gf/logs/msmc2_size_e.out
#SBATCH --error=/scratch1/marjanak/gvcf6/gf/logs/msmc2_size_e.err
#SBATCH --partition=qcb
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=100G

# Purge existing modules
module purge

# Activate the Conda environment for MSMC2
eval "$(conda shell.bash hook)"
conda activate /home1/marjanak/.conda/envs/mymsmc2

# Define directories
INPUTDIR=/scratch1/marjanak/gvcf6/gf/multihet    # Path to multihetsep files
OUTDIR=/scratch1/marjanak/gvcf6/gf/results      # Path to save MSMC2 output

# MSMC2 command for East population on all chromosomes
msmc2_Linux -t 10 \
    -o $OUTDIR/east.msmc2 \
    -I 0-1,2-3,4-5,6-7,8-9,10-11 \
    $INPUTDIR/gf_east_chr1.multihetsep.txt \
    $INPUTDIR/gf_east_chr2.multihetsep.txt \
    $INPUTDIR/gf_east_chr3.multihetsep.txt \
    $INPUTDIR/gf_east_chr4.multihetsep.txt \
    $INPUTDIR/gf_east_chr5.multihetsep.txt \
    $INPUTDIR/gf_east_chr6.multihetsep.txt \
    $INPUTDIR/gf_east_chr7.multihetsep.txt \
    $INPUTDIR/gf_east_chr8.multihetsep.txt \
    $INPUTDIR/gf_east_chr9.multihetsep.txt \
    $INPUTDIR/gf_east_chr10.multihetsep.txt \
    $INPUTDIR/gf_east_chr11.multihetsep.txt \
    $INPUTDIR/gf_east_chr12.multihetsep.txt \
    $INPUTDIR/gf_east_chr13.multihetsep.txt \
    $INPUTDIR/gf_east_chr14.multihetsep.txt \
    $INPUTDIR/gf_east_chr15.multihetsep.txt \
    $INPUTDIR/gf_east_chr16.multihetsep.txt \
    $INPUTDIR/gf_east_chr17.multihetsep.txt \
    $INPUTDIR/gf_east_chr18.multihetsep.txt \
    $INPUTDIR/gf_east_chr19.multihetsep.txt \
    $INPUTDIR/gf_east_chr20.multihetsep.txt \
    $INPUTDIR/gf_east_chr21.multihetsep.txt \
    $INPUTDIR/gf_east_chr22.multihetsep.txt \
    $INPUTDIR/gf_east_chr23.multihetsep.txt \
    $INPUTDIR/gf_east_chr24.multihetsep.txt \
    $INPUTDIR/gf_east_chr25.multihetsep.txt \
    $INPUTDIR/gf_east_chr26.multihetsep.txt \
    $INPUTDIR/gf_east_chr27.multihetsep.txt \
    $INPUTDIR/gf_east_chr28.multihetsep.txt \
    $INPUTDIR/gf_east_chr29.multihetsep.txt \
    $INPUTDIR/gf_east_chr30.multihetsep.txt \
    $INPUTDIR/gf_east_chr31.multihetsep.txt \
    $INPUTDIR/gf_east_chr32.multihetsep.txt

