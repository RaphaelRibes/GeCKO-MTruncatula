#!/bin/bash
#SBATCH --gres=gpu:a100:1

module () {
    eval `/usr/bin/modulecmd bash $*`
}

pbrun indexgvcf --input "$1"
