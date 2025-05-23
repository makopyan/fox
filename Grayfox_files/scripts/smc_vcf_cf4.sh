#!/bin/sh
#SBATCH --job-name=smc_vcf_cf4
#SBATCH --output=smc_vcf_cf4_%a.out
#SBATCH --error=smc_vcf_cf4_%a.err
#SBATCH --array=1-38

module purge

eval "$(conda shell.bash hook)"

conda activate /home1/marjanak/.conda/envs/smcpp_env

export OMP_NUM_THREADS=2

# Define VCF and map files
VCF_FILE="Canfam4_filtered.renameChroms.Mainland.drop295.QUAL30_DPgr205lt500_ANgr59.vcf.gz"
BED_FILE="canfam4.rep.map.genic.bed.gz"

# Define arrays for WES and EAS samples
WES_SAMPLES=("SRR24465296" "SRR24465287" "SRR24465288" "SRR24465284" "SRR24465292" "SRR24465290")
EAS_SAMPLES=("SRR24465286" "SRR24465275" "SRR24465272" "SRR24465271" "SRR24465305" "SRR24465304")

# Loop through WES samples
for SAMPLE in "${WES_SAMPLES[@]}"; do
    WES_OTHER_SAMPLES=("${WES_SAMPLES[@]}")  
    smc++ vcf2smc \
    -d $SAMPLE $SAMPLE -m $MAP_FILE $VCF_FILE wes_${SAMPLE}_chr${SLURM_ARRAY_TASK_ID}.cf4.smc.gz \
    chr${SLURM_ARRAY_TASK_ID} WES:${WES_OTHER_SAMPLES[*]}
done

# Loop through EAS samples
for SAMPLE in "${EAS_SAMPLES[@]}"; do
    EAS_OTHER_SAMPLES=("${EAS_SAMPLES[@]}")  
    smc++ vcf2smc \
    -d $SAMPLE $SAMPLE -m $MAP_FILE $VCF_FILE eas_${SAMPLE}_chr${SLURM_ARRAY_TASK_ID}.cf4.smc.gz \
    chr${SLURM_ARRAY_TASK_ID} EAS:${EAS_OTHER_SAMPLES[*]}
done