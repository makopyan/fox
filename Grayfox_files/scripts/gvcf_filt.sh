#!/bin/bash
#SBATCH --job-name=gvcf_filt
#SBATCH --output=/scratch1/marjanak/gvcf6/gvcf_filt_%A_%a.out
#SBATCH --error=/scratch1/marjanak/gvcf6/gvcf_filt_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --partition=qcb
#SBATCH --array=0-7 # 8 tasks, one for each combination
#SBATCH --cpus-per-task=8

module purge
module load gcc/11.3.0
module load openblas/0.3.20
module load bcftools

# Define your arrays for prefixes and GVCF paths
prefixes=("gf" "af" "cf3" "cf4")
gvcfs=(
    "/project/jazlynmo_738/DataRepository/Canids/Invariant/GrayFox/Mainland/grayfox_filtered.renameChroms.Mainland.drop295.QUAL30_DPgr205lt500_ANgr59.gvcf.gz" \
    "/project/jazlynmo_738/DataRepository/Canids/Invariant/GrayFox/Mainland/arcticfox_filtered.renameChroms.Mainland.drop295.QUAL30_DPgr205lt500_ANgr59.gvcf.gz" \
    "/project/jazlynmo_738/DataRepository/Canids/Invariant/GrayFox/Mainland/Canfam3.1_filtered.renameChroms.Mainland.drop295.QUAL30_DPgr205lt500_ANgr59.gvcf.gz" \
    "/project/jazlynmo_738/DataRepository/Canids/Invariant/GrayFox/Mainland/Canfam4_filtered.renameChroms.Mainland.drop295.QUAL30_DPgr205lt500_ANgr59.gvcf.gz"
)
populations=("east" "west")

# Calculate which combination this task should run based on array index
prefix_idx=$((SLURM_ARRAY_TASK_ID / 2))
pop_idx=$((SLURM_ARRAY_TASK_ID % 2))

prefix="${prefixes[$prefix_idx]}"
gvcf="${gvcfs[$prefix_idx]}"
population="${populations[$pop_idx]}"

# Define sample file based on population
if [ "$population" == "east" ]; then
    sample_file="gf_east_n6.txt"
else
    sample_file="gf_west_n6.txt"
fi

# Set AN threshold to 12
an=12

# Step 1: Filter samples and apply AN filter
bcftools view -S "$sample_file" "$gvcf" | bcftools filter -i "AN >= $an" -Oz > "${prefix}.${population}.AN${an}.gvcf.gz"

# Step 2: Create BED file
bcftools query -f '%CHROM\t%POS\n' "${prefix}.${population}.AN${an}.gvcf.gz" > "${prefix}.${population}.mask.bed"

# Step 3: Index the resulting filtered GVCF
bcftools index --threads 8 -t "${prefix}.${population}.AN${an}.gvcf.gz"

echo "Processing and indexing completed for ${prefix} ${population}"

