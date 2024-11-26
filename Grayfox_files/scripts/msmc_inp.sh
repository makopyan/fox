#!/bin/sh
#SBATCH --job-name=msmc_gf_inp
#SBATCH --output=/scratch1/marjanak/gvcf6/gf/logs/msmc_inp_%A_%a.out
#SBATCH --error=/scratch1/marjanak/gvcf6/gf/logs/msmc_inp_%A_%a.err
#SBATCH --partition=qcb
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8000MB
#SBATCH --array=1-32  # Job array for 32 chromosomes

# Purge existing modules
module purge

# Activate the Conda environment for MSMC2
eval "$(conda shell.bash hook)"
conda activate /home1/marjanak/.conda/envs/mymsmc2

chrom="$SLURM_ARRAY_TASK_ID"

generate_multihetsep.py \
--mask=/scratch1/marjanak/gvcf6/gf/gf_west_bed/gf.west.mask.chr${chrom}.bed \
--mask=/scratch1/marjanak/gvcf6/gf/gf_west_bed/gf.west.mask.chr${chrom}.bed \
--mask=/scratch1/marjanak/gvcf6/gf/gf_west_bed/gf.west.mask.chr${chrom}.bed \
--mask=/scratch1/marjanak/gvcf6/gf/gf_west_bed/gf.west.mask.chr${chrom}.bed \
--mask=/scratch1/marjanak/gvcf6/gf/gf_west_bed/gf.west.mask.chr${chrom}.bed \
--mask=/scratch1/marjanak/gvcf6/gf/gf_west_bed/gf.west.mask.chr${chrom}.bed \
--negative_mask=/scratch1/marjanak/gvcf6/gf/neg_masks/chr${chrom}_mask.bed.gz \
/scratch1/marjanak/gvcf6/gf/split_gvcf/gf_west_SRR24465284_chr${chrom}.gvcf.gz \
/scratch1/marjanak/gvcf6/gf/split_gvcf/gf_west_SRR24465287_chr${chrom}.gvcf.gz \
/scratch1/marjanak/gvcf6/gf/split_gvcf/gf_west_SRR24465288_chr${chrom}.gvcf.gz \
/scratch1/marjanak/gvcf6/gf/split_gvcf/gf_west_SRR24465290_chr${chrom}.gvcf.gz \
/scratch1/marjanak/gvcf6/gf/split_gvcf/gf_west_SRR24465292_chr${chrom}.gvcf.gz \
/scratch1/marjanak/gvcf6/gf/split_gvcf/gf_west_SRR24465296_chr${chrom}.gvcf.gz > multihet/gf_west_chr${chrom}.multihetsep.txt

generate_multihetsep.py \
--mask=/scratch1/marjanak/gvcf6/gf/gf_east_bed/gf.east.mask.chr${chrom}.bed \
--mask=/scratch1/marjanak/gvcf6/gf/gf_east_bed/gf.east.mask.chr${chrom}.bed \
--mask=/scratch1/marjanak/gvcf6/gf/gf_east_bed/gf.east.mask.chr${chrom}.bed \
--mask=/scratch1/marjanak/gvcf6/gf/gf_east_bed/gf.east.mask.chr${chrom}.bed \
--mask=/scratch1/marjanak/gvcf6/gf/gf_east_bed/gf.east.mask.chr${chrom}.bed \
--mask=/scratch1/marjanak/gvcf6/gf/gf_east_bed/gf.east.mask.chr${chrom}.bed \
--negative_mask=/scratch1/marjanak/gvcf6/gf/neg_masks/chr${chrom}_mask.bed.gz \
/scratch1/marjanak/gvcf6/gf/split_gvcf/gf_east_SRR24465286_chr${chrom}.gvcf.gz \
/scratch1/marjanak/gvcf6/gf/split_gvcf/gf_east_SRR24465275_chr${chrom}.gvcf.gz \
/scratch1/marjanak/gvcf6/gf/split_gvcf/gf_east_SRR24465272_chr${chrom}.gvcf.gz \
/scratch1/marjanak/gvcf6/gf/split_gvcf/gf_east_SRR24465271_chr${chrom}.gvcf.gz \
/scratch1/marjanak/gvcf6/gf/split_gvcf/gf_east_SRR24465305_chr${chrom}.gvcf.gz \
/scratch1/marjanak/gvcf6/gf/split_gvcf/gf_east_SRR24465304_chr${chrom}.gvcf.gz > multihet/gf_east_chr${chrom}.multihetsep.txt
