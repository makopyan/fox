#!/bin/sh
#SBATCH --job-name=gf_east_pyrho2
#SBATCH --output=/scratch1/genchev/logfile/gf_east_pyrho2_%a.out  
#SBATCH --error=/scratch1/genchev/logfile/gf_east_pyrho2_%a.err
#SBATCH --time=03:00:00
#SBATCH --account=jazlynmo_738
#SBATCH --partition=main #partition to run the job on 
#SBATCH --ntasks=1
#SBATCH --array=1-32
#SBATCH --cpus-per-task=5 #one cpu per task
#SBATCH --mem-per-cpu=5000MB #this is equivalent to 5G
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=genchev@usc.edu


#load all the modules and start up enviornment
module purge
eval "$(conda shell.bash hook)"
module load conda
module load gcc/8.3.0
module load python/3.6.8

mamba activate my-pyrho-env #activate your environment change as needed




# UPDATE !!!! SUPER IMPORTANT !!!!!! UPDATE THE NAME!!!!!!!!
NAME="grayfox_east"
CHROM=$SLURM_ARRAY_TASK_ID
#CHROM=32 #used this to test

# n (sample size)
declare -i nw=12
declare -i ne=12


# Determined through analysis of hyperparameter optimizations
# Block penalty
declare -i BP=25

# Window size
declare -i WS=50




pyrho optimize --tablefile "/home1/genchev/grayfox_rmap/"$NAME".hdf" \
	--vcffile "/home1/genchev/split_vcf_n6/chr"$CHROM".east.gf.vcf.gz"  \
	--outfile "/home1/genchev/grayfox_rmap/"$NAME"_"$CHROM".rmap" \
	--blockpenalty $BP --windowsize $WS \
	--ploidy 2 --numthreads 5

echo "optimize done!"


# pyrho also has a utility to compute the distribution of r^2, a measure
# of linkage-disequillibrium from the lookup tables. This can be of interest
# to see if inferred demographies or recombination maps fit the data well.
# Note that the distribution of r^2 depends heavily on the demographic
# history.

pyrho compute_r2 --quantiles .25,.5,.75 --compute_mean --samplesize $ne \
	--tablefile "/home1/genchev/grayfox_rmap/"$NAME".hdf" \
	--outfile "/home1/genchev/grayfox_rmap/"$NAME"_"$CHROM"_r2.tsv"

echo "compute_r2 done!"

sleep 180 #sleep for 3 minutes no matter what
