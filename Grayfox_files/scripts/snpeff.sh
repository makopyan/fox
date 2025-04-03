#!/bin/sh
#SBATCH --job-name=snpeff
#SBATCH --output=/scratch1/marjanak/mash/snpeff_%A_%a.out
#SBATCH --error=/scratch1/marjanak/mash/snpeff_%A_%a.err
#SBATCH --partition=qcb
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=30G
#SBATCH --array=0-2  # Array for three reference genomes

module purge
eval "$(conda shell.bash hook)"
conda activate /home1/marjanak/.conda/envs/snpeff_env

# Define arrays for references, annotations, and VCFs
GENOMES=("arcticfox" "canfam4" "grayfox")
FASTA_FILES=(
    "/scratch1/marjanak/mash/arcticfox_genome.fna"
    "/scratch1/marjanak/mash/canfam4_genome.fna"
    "/scratch1/marjanak/mash/grayfox_genome.fna"
)
ANNOTATIONS=(
    "/scratch1/marjanak/mash/arcticfox.chr.gtf"
    "/scratch1/marjanak/mash/canfam4.chr.gtf"
    "/scratch1/marjanak/mash/grayfox.chr.gtf"
)
VCFS=(
    "/scratch1/marjanak/mash/arcticfox.n12.vcf.gz"
    "/scratch1/marjanak/mash/canfam4.n12.vcf.gz"
    "/scratch1/marjanak/mash/grayfox.n12.vcf.gz"
)

# Select genome, annotation, and VCF based on SLURM array index
GENOME=${GENOMES[$SLURM_ARRAY_TASK_ID]}
FASTA=${FASTA_FILES[$SLURM_ARRAY_TASK_ID]}
ANNOTATION=${ANNOTATIONS[$SLURM_ARRAY_TASK_ID]}
VCF=${VCFS[$SLURM_ARRAY_TASK_ID]}

# SnpEff setup
SNP_EFF_DIR="/scratch1/marjanak/mash/snpeff_data"
mkdir -p $SNP_EFF_DIR/data/$GENOME

# Copy genome and annotation for SnpEff database
cp $FASTA $SNP_EFF_DIR/data/$GENOME/sequences.fa
cp $ANNOTATION $SNP_EFF_DIR/data/$GENOME/genes.gtf

# Add genome to SnpEff config (prevent duplicates)
CONFIG_FILE="$SNP_EFF_DIR/snpEff.config"
grep -q "^$GENOME.genome" $CONFIG_FILE || echo "$GENOME.genome : $GENOME" >> $CONFIG_FILE

# Build the SnpEff database
cd $SNP_EFF_DIR
snpEff -Xmx30G build -gtf22 -noCheckProtein -noCheckCds -v $GENOME

# Annotate the VCF
OUT_VCF="/scratch1/marjanak/mash/${GENOME}_annotated.vcf"
snpEff -v $GENOME $VCF > $OUT_VCF

# Extract synonymous and non-synonymous mutations
grep -E "synonymous_variant" $OUT_VCF > "/scratch1/marjanak/mash/${GENOME}_synonymous.txt"
grep -E "missense_variant" $OUT_VCF > "/scratch1/marjanak/mash/${GENOME}_nonsynonymous.txt"

echo "SnpEff annotation completed for $GENOME"
