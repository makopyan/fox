### Install paf2chain
#git clone https://github.com/AndreaGuarracino/paf2chain
#cd paf2chain
#cargo install --force --path


### Install crossmap
#conda install bioconda::crossmap


## Get pairwise alignment of dog (Canfam4) and Arctic fox genomes to grayfox genome from Ellie

cp canfam4-grayfox.paf cf4-gf-chr.paf

cp arcticfox-grayfox.paf af-gf-chr.paf


## Make pairwise alignment between Canfam4 and Arctic fox using Minimap2

minimap2 -x asm20 canfam4_genome.fna arcticfox_genome.fna -t 8 > aln_af_cf4_asm20.paf


## Rename scaffolds to chromosomes

bash cf4.gf.rename.sed cf4-gf-chr.paf

bash af.gf.rename.sed af-gf-chr.paf


## Convert renamed paf to chain

paf2chain -i cf4-gf-chr.paf > cf4-gf-chr.chain

paf2chain -i af-gf-chr.paf > af-gf-chr.chain

paf2chain -i aln_af_cf4_asm20.paf > aln_af_cf4_asm20.chain




## Get SNPs per pop and genome


module load gcc/11.3.0  
module load openblas/0.3.20
module load bcftools


## Filter vcf for 6 west samples, with no missing data

bcftools view -Oz -S gf_west_n6.txt /project/jazlynmo_738/DataRepository/Canids/Variants/GrayFox/Mainland/grayfox_filtered.renameChroms.Mainland.drop295.QUAL30_DPgr205lt500_ANgr59.vcf.gz > west6.vcf.gz

bcftools index -t west6.vcf.gz

bcftools view --min-ac=1 --max-ac 11 west6.vcf.gz | bcftools view -i 'F_MISSING<0.1' -Oz -o snpwest6.vcf.gz

bcftools index -t snpwest6.vcf.gz


## Get snp positions and covert to bed format

bcftools query -f '%CHROM %POS\n' snpwest6.vcf.gz  > snpwest.txt

awk 'OFS="\t"{print $1, $2-1, $2}' snpwest.txt > snpwest.bed


## Filter vcf for 6 east samples, with no missing data

bcftools view -Oz -S gf_east_n6.txt /project/jazlynmo_738/DataRepository/Canids/Variants/GrayFox/Mainland/grayfox_filtered.renameChroms.Mainland.drop295.QUAL30_DPgr205lt500_ANgr59.vcf.gz > east6.vcf.gz

bcftools index -t east6.vcf.gz

bcftools view --min-ac=1 --max-ac 11 east6.vcf.gz | bcftools view -i 'F_MISSING<0.1' -Oz -o snpeast6.vcf.gz

bcftools index -t snpeast6.vcf.gz


## Get snp positions and covert to bed format

bcftools query -f '%CHROM %POS\n' snpeast6.vcf.gz  > snpeast.txt

awk 'OFS="\t"{print $1, $2-1, $2}' snpeast.txt > snpeast.bed




eval "$(conda shell.bash hook)"

conda activate crossmap



## Convert snp positions

CrossMap bed af-gf-chr.chain snpeast.bed | tee af_gf_east.txt

CrossMap bed af-gf-chr.chain snpwest.bed | tee af_gf_west.txt

CrossMap bed cf4-gf-chr.chain snpeast.bed | tee cf4_gf_east.txt

CrossMap bed cf4-gf-chr.chain snpwest.bed | tee cf4_gf_west.txt


## Liftover 50kb FST regions to grayfox

CrossMap region cf4-gf-chr.chain cf4_fst_50kb.txt | tee cf4_gf_fst.txt

CrossMap region af-gf-chr.chain af_fst_50kb.txt | tee af_gf_fst.txt


## Liftover 50kb FST regions between heterospecifics

CrossMap region /scratch1/marjanak/mash/aln_af_cf4_asm20.chain af_fst_50kb.txt  | tee af_cf4_fst.txt



## Liftover 50kb recombination regions to grayfox

CrossMap region cf4-gf-chr.chain cf4_recomb50kb.txt | tee cf4_gf_recomb.txt

CrossMap region af-gf-chr.chain af_recomb50kb.txt | tee af_gf_recomb.txt
