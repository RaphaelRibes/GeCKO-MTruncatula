#!/bin/bash
#SBATCH --gres=gpu:a100:1

module () {
    eval `/usr/bin/modulecmd bash $*`
}

pbrun haplotypecaller --ref "$1" --in-bam "$2" --out-variants "$3" # --haplotypecaller-options "$4"
