#!/bin/sh
#SBATCH --output=dispatcher.out
#SBATCH --error=dispatcher.err
#SBATCH --mem=1G

module purge
module load singularity/3.5

pixi run ./dispatcher.sh -r -m -i F11023