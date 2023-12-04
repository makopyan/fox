variants
================

start with two vcf files, one for gray foxes and one for island foxes
with variants called using the gray fox reference genome

Copy files to my scratch (/scratch1/marjanak/) from OG directory:

`/project/jazlynmo_738/DataRepository/Canids/Variants/GrayFox/Mainland/grayfox_filtered.renameChroms.Mainland.ACgr61_DPgr165lt500.vcf.gz`

and

`/project/jazlynmo_738/DataRepository/Canids/Variants/GrayFox/ChannelIsland/grayfox_filtered.renameChroms.ACgr25_DPgr165lt500.vcf.gz`

Run *vcf-stats* module from vcftools with shell script `vcftools.sh`:  
<details>
<summary>
Show code
</summary>
<p>

``` bash
#!/bin/sh
#SBATCH --job-name=vcfmainland
#SBATCH --output=/scratch1/marjanak/vcfACgr61.out
#SBATCH --error=/scratch1/marjanak/vcfACgr61.err
#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8000MB
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=marjanak@usc.edu

module load vcftools

vcf-stats grayfox_filtered.renameChroms.Mainland.ACgr61_DPgr165lt500.vcf.gz
```

</p>
</details>

Command: `sbatch ./vcfstats.sh`

#### Mainland Fox

Includes 42 samples and 602,838 variants

#### Island Fox

Includes 16 samples and 4,160,435 variants
