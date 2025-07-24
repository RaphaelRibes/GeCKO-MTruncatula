#!/bin/sh
#SBATCH --output=dispatcher.out
#SBATCH --error=dispatcher.err

pixi run ./dispatcher.sh -r -m -i F11023