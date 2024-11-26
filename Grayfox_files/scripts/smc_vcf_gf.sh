#!/bin/sh
#SBATCH --job-name=smc_vcf_gf
#SBATCH --output=smc_vcf_gf_%a.out
#SBATCH --error=smc_vcf_gf_%a.err
#SBATCH --array=1-32

module purge

eval "$(conda shell.bash hook)"

conda activate /home1/marjanak/.conda/envs/smcpp_env

export OMP_NUM_THREADS=2

# Define VCF and map files
VCF_FILE="grayfox_filtered.renameChroms.Mainland.drop295.QUAL30_DPgr205lt500_ANgr59.vcf.gz"
BED_FILE="grayfox.rep.map.genic.bed.gz"

# Define arrays for WES and EAS samples
WES_SAMPLES=("SRR24465296" "SRR24465287" "SRR24465288" "SRR24465284" "SRR24465292" "SRR24465290")
EAS_SAMPLES=("SRR24465286" "SRR24465275" "SRR24465272" "SRR24465271" "SRR24465305" "SRR24465304")

# Loop through WES samples
for SAMPLE in "${WES_SAMPLES[@]}"; do
    # Dynamically create a list of WES samples excluding the current SAMPLE
    WES_OTHER_SAMPLES=("${WES_SAMPLES[@]/$SAMPLE}")  # Remove $SAMPLE from the array
    smc++ vcf2smc \
    -d $SAMPLE $SAMPLE -m $MAP_FILE $VCF_FILE wes_${SAMPLE}_chr${SLURM_ARRAY_TASK_ID}.gf.smc.gz \
    chr${SLURM_ARRAY_TASK_ID} WES:${WES_OTHER_SAMPLES[*]}  # Use the updated list
done

# Loop through EAS samples
for SAMPLE in "${EAS_SAMPLES[@]}"; do
    # Dynamically create a list of EAS samples excluding the current SAMPLE
    EAS_OTHER_SAMPLES=("${EAS_SAMPLES[@]/$SAMPLE}")  # Remove $SAMPLE from the array
    smc++ vcf2smc \
    -d $SAMPLE $SAMPLE -m $MAP_FILE $VCF_FILE eas_${SAMPLE}_chr${SLURM_ARRAY_TASK_ID}.gf.smc.gz \
    chr${SLURM_ARRAY_TASK_ID} EAS:${EAS_OTHER_SAMPLES[*]}  # Use the updated list
done