##Download, trim, and align

1. Run dlfa.sh to download raw fastq from SRA
2. Run trim.sh to trim adapter reads using bbduk.sh from bbmap
3. Download grayfox genome and run bwa index
4. Align trimmed reads to indexed genome


5. Merge aligned bams and remove dupicates

##Genotype (GATK)

6. Run HaplotypeCaller in GATK to generate per-sample GVCFs

7. Run GenomicsDBImport to input GVCFs for joint genotyping

8. Run GenotypeGVCFs for joint genotyping of all samples


##Filter (bcftools)

9. Filter out repeat regions and low mappability regions

10. Rename chromosomes and output all biallelic sites

