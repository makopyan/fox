#!/bin/sh
#SBATCH --job-name=split_vcfs
#SBATCH --output=/scratch1/marjanak/gvcf6/logs/split_gvcf_%A_%a.out
#SBATCH --error=/scratch1/marjanak/gvcf6/logs/split_gvcf_%A_%a.err
#SBATCH --partition=qcb
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8000MB
#SBATCH --array=1-32  # Job array for 32 chromosomes

module purge
module load gcc/11.3.0 openblas/0.3.20 bcftools

# Define chromosome list (32 chromosomes)
chromosomes=({1..32})

# Read sample lists from alleast.txt and allwest.txt
samples_east=($(cat /scratch1/marjanak/gvcf6/gf_east_n6.txt))
samples_west=($(cat /scratch1/marjanak/gvcf6/gf_west_n6.txt))

# Get the chromosome for the current task
chrom=${chromosomes[$SLURM_ARRAY_TASK_ID-1]}

# Define input GVCF files and prefixes
east_gvcf="gf.east.AN12.gvcf.gz"
west_gvcf="gf.west.AN12.gvcf.gz"
east_prefix="gf_east"
west_prefix="gf_west"

# Loop through each sample in the east and split by chromosome
for sample in "${samples_east[@]}"; do
  # Create the split GVCF file for the east samples and chromosome
  bcftools view -r chr${chrom} -s ${sample} $east_gvcf -Oz -o ${east_prefix}_${sample}_chr${chrom}.gvcf.gz
  
  # Index the resulting east GVCF file
  bcftools index -t ${east_prefix}_${sample}_chr${chrom}.gvcf.gz
done

# Loop through each sample in the west and split by chromosome
for sample in "${samples_west[@]}"; do
  # Create the split GVCF file for the west samples and chromosome
  bcftools view -r chr${chrom} -s ${sample} $west_gvcf -Oz -o ${west_prefix}_${sample}_chr${chrom}.gvcf.gz
  
  # Index the resulting west GVCF file
  bcftools index -t ${west_prefix}_${sample}_chr${chrom}.gvcf.gz
done
