#!/bin/bash
#SBATCH --job-name=sfs
#SBATCH --output=/home1/marjanak/n6/sfs_%A_%a.out
#SBATCH --error=/home1/marjanak/n6/sfs_%A_%a.err
#SBATCH --time=3:00:00
#SBATCH --partition=qcb
#SBATCH --array=0-5 # 6 tasks: 2 populations x 3 VCFs
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8000MB

module purge
module load gcc/11.3.0
module load openblas/0.3.20
module load bcftools
module load vcftools

# Define arrays for VCFs, populations, and their corresponding genome identifiers
vcfs=(
    "/project/jazlynmo_738/DataRepository/Canids/Variants/GrayFox/Mainland/arcticfox_filtered.renameChroms.Mainland.drop295.QUAL30_DPgr205lt500_ANgr59.vcf.gz"
    "/project/jazlynmo_738/DataRepository/Canids/Variants/GrayFox/Mainland/grayfox_filtered.renameChroms.Mainland.drop295.QUAL30_DPgr205lt500_ANgr59.vcf.gz"
    "/project/jazlynmo_738/DataRepository/Canids/Variants/GrayFox/Mainland/Canfam4_filtered.renameChroms.Mainland.drop295.QUAL30_DPgr205lt500_ANgr59.vcf.gz"
)
pops=("gf_east_n6.txt" "gf_west_n6.txt")
pop_prefixes=("east" "west")
genome_prefixes=("af" "gf" "cf4")

# Calculate which VCF and population this task should run based on the array index
vcf_idx=$((SLURM_ARRAY_TASK_ID / 2))
pop_idx=$((SLURM_ARRAY_TASK_ID % 2))

vcf="${vcfs[$vcf_idx]}"
pop_file="${pops[$pop_idx]}"
pop_prefix="${pop_prefixes[$pop_idx]}"
genome_prefix="${genome_prefixes[$vcf_idx]}"
out_prefix="${pop_prefix}_${genome_prefix}"

# Step 1: Filter using bcftools based on the population file
bcftools view -Oz -S "$pop_file" "$vcf" > ${out_prefix}.vcf.gz

# Step 2: Index the VCF file
bcftools index -t ${out_prefix}.vcf.gz

# Step 3: Run vcftools counts
vcftools --gzvcf ${out_prefix}.vcf.gz --counts --out ${out_prefix} >& ${out_prefix}.log

# Step 4: Convert the count file to a tab-delimited format
sed -e 's/:/\t/g' ${out_prefix}.frq.count > ${out_prefix}.count

# Step 5: Create summary for each count file
fixed=0
singleton=0
others=0

awk 'NR>1 { 
    count1=$6+0; count2=$8+0; 
    if (count1 == 1 || count2 == 1) 
        singleton++;
    else if (count1 == 0 || count2 == 0) 
        fixed++;
    else 
        others++;
} 
END {print fixed, singleton, others}' ${out_prefix}.count > ${out_prefix}_snp_summary.txt

# Step 6: Calculate the Site Frequency Spectrum (SFS)
output_dir="sfs_output"
mkdir -p $output_dir

# Number of genotypes in each population
n_chr=12

# Initialize an associative array to store SFS bin counts
declare -A sfs

# Read the count file and process counts
awk -v n_chr=$n_chr '
BEGIN { FS="\t"; OFS="\t" }
NR > 1 && $4 == n_chr {
    # Bin by the minimum and maximum counts between count1 and count2
    min_count = ($6 < $8) ? $6 : $8
    max_count = ($6 >= $8) ? $6 : $8
    bin = min_count "," max_count
    sfs[bin]++
}
END {
    for (bin in sfs) {
        print bin, sfs[bin]
    }
}' ${out_prefix}.count > $output_dir/${out_prefix}_sfs.txt

# Notify user of completion
echo "SFS calculations completed for ${out_prefix}. Results saved in $output_dir/${out_prefix}_sfs.txt"

# Unset the associative array to reset for next population
unset sfs

echo "Analysis completed for ${out_prefix}"
