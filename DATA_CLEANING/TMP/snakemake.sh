#!/bin/bash 
#
#SBATCH -J DEMULT
#SBATCH -o DEMULT."%j".out
#SBATCH -e DEMULT."%j".err 

#SBATCH -p agap_normal

module purge
module load snakemake/5.13.0 
module load python/Anaconda/3-5.1.0
mkdir -p logs
snakePipe=$1 



#snakemake --snakefile data_cleaning.smk --dryrun --dag --cluster-config cluster_config.yml  --configfile config.yml --jobs 200 --printshellcmds --use-envmodules  --cluster "sbatch  -p {cluster.partition} --ntasks {cluster.ntasks} --mem-per-cpu={cluster.mem-per-cpu} -c {cluster.cpus-per-task} -e {cluster.error} -o {cluster.output}"
#snakemake --snakefile data_cleaning.smk --dryrun --dag --cluster-config cluster_config.yml  --configfile config.yml
snakemake  --use-singularity --use-conda --snakefile $snakePipe --configfile config.yml --jobs 200 --printshellcmds --use-envmodules 
