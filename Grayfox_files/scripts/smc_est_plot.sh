#!/bin/sh
#SBATCH --job-name=smc_est
#SBATCH --output=smc_est.out
#SBATCH --error=smc_est.err

module purge

eval "$(conda shell.bash hook)"

conda activate /home1/marjanak/.conda/envs/smcpp_env

export OMP_NUM_THREADS=2


smc++ estimate 45e-10 wes_*.gf.smc.gz --base gf_west --timepoints 2.5e3 5e5 -v --cores 20 --thinning 1792 --knots 12

smc++ estimate 45e-10 eas_*.gf.smc.gz --base gf_east --timepoints 2.5e3 5e5 -v --cores 20 --thinning 1792 --knots 12

smc++ estimate 45e-10 wes_*.af.smc.gz --base af_west --timepoints 2.5e3 5e5 -v --cores 20 --thinning 1792 --knots 12

smc++ estimate 45e-10 eas_*.af.smc.gz --base af_east --timepoints 2.5e3 5e5 -v --cores 20 --thinning 1792 --knots 12

smc++ estimate 45e-10 wes_*.cf4.smc.gz --base cf4_west --timepoints 2.5e3 5e5 -v --cores 20 --thinning 1792 --knots 12

smc++ estimate 45e-10 eas_*.cf4.smc.gz --base cf4_east --timepoints 2.5e3 5e5 -v --cores 20 --thinning 1792 --knots 12


smc++ plot -c gf_east_plot.pdf gf_east.final.json -g 2

smc++ plot -c gf_west_plot.pdf gf_west.final.json -g 2

smc++ plot -c af_east_plot.pdf af_east.final.json -g 2

smc++ plot -c af_west_plot.pdf af_west.final.json -g 2

smc++ plot -c cf4_east_plot.pdf cf4_east.final.json -g 2

smc++ plot -c cf4_west_plot.pdf cf4_west.final.json -g 2