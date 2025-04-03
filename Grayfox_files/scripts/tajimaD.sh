#!/bin/sh
#SBATCH --account=jazlynmo_738
#SBATCH --partition=qcb
#SBATCH --job-name=tajD
#SBATCH --output=/scratch1/marjanak/tajD/logs/tajD.out
#SBATCH --error=/scratch1/marjanak/tajD/logs/tajD.err
#SBATCH --time=6:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 
#SBATCH --array=1-38  # Set to the max chromosome count (38 from cf4)

module purge
module load gcc/13.3.0 vcftools/0.1.16 


# Current chromosome number
chr=${SLURM_ARRAY_TASK_ID}

# Process each population-genome combination
declare -A genome_chr_counts=(
  ["af"]=24
  ["cf4"]=38
  ["gf"]=32
)

for pop in "east" "west"; do
  for genome in "af" "cf4" "gf"; do
    # Skip if current chromosome is beyond this genome's chromosome count
    max_chr=${genome_chr_counts[$genome]}
    if [ $chr -gt $max_chr ]; then
      continue
    fi
    
    input_file="${pop}.${genome}.vcf.gz"
    
    # Make sure output directory exists
    mkdir -p ${genome}
    
    # Run vcftools
    vcftools --gzvcf /project/jazlynmo_738/share_data/refbias/split_vcf_n6/${input_file} \
      --TajimaD 50000 \
      --max-missing 1 \
      --chr chr${chr} \
      --out ${genome}/${pop}_${genome}_tajD_50kb_chr${chr}
  done
done