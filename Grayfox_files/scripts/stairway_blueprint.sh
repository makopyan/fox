#!/bin/bash
#SBATCH --job-name=sw_blueprint
#SBATCH --output=sw_%A_%a.out
#SBATCH --error=sw_%A_%a.err
#SBATCH --time=20:00:00
#SBATCH --partition=qcb
#SBATCH --array=0-5

##blueprint files are located in the input directory

# Load the required module
module load jdk/17.0.5

# Define arrays for blueprint files and scripts
blueprints=(
    "af_east_single.blueprint"
    "cf4_east_single.blueprint"
    "gf_east_single.blueprint"
    "af_west_single.blueprint"
    "cf4_west_single.blueprint"
    "gf_west_single.blueprint"
)

scripts=(
    "af_east_single.blueprint.sh"
    "cf4_east_single.blueprint.sh"
    "gf_east_single.blueprint.sh"
    "af_west_single.blueprint.sh"
    "cf4_west_single.blueprint.sh"
    "gf_west_single.blueprint.sh"
)

# Get the current blueprint file and script based on the array index
blueprint=${blueprints[$SLURM_ARRAY_TASK_ID]}
script=${scripts[$SLURM_ARRAY_TASK_ID]}

# Step 1: Run Java Stairbuilder for the blueprint file
java -cp stairway_plot_es Stairbuilder "$blueprint"

# Step 2: Ensure the script is executable
chmod +x "$script"

# Step 3: Run the corresponding script
"$script"