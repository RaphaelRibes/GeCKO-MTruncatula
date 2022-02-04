#!/bin/bash
#
#SBATCH -J DataCleaningSmk
#SBATCH -o DataCleaningSmk."%j".out
#SBATCH -e DataCleaningSmk."%j".err

#SBATCH -p agap_normal

module purge
module load snakemake/5.13.0
module load anaconda/python3.8

mkdir -p Logs_DataCleaningWorkflow

PRGR=/home/girodollej/CAPTURE_PIPELINES_SNAKEMAKE/DATA_CLEANING/WORKFLOW
CONFIG=/home/girodollej/scratch/DATA_TEST/VIR_Cap001/CONFIG


## DRYRUN ##
#snakemake --snakefile ${PRGR}/DataCleaning.smk -n --dag --forceall --cluster-config ${CONFIG}/cluster_config.yml --configfile ${CONFIG}/config.yml --printshellcmds

## GENERATE PIPELINE DIAGRAM -> à lancer hors sbatch, dot est indisponible sur les noeuds de calcul ##
#snakemake --snakefile ${PRGR}/DataCleaning.smk -n --dag --forceall --cluster-config ${CONFIG}/cluster_config.yml --configfile ${CONFIG}/mini_tmp/config_mini.yml | dot -Tsvg > data_cleaning_pipeline_mini.svg



                                        ### CLASSIC CLUSTER CONFIG ###

## RUN WITH ENVMODULES ##
#snakemake --snakefile ${PRGR}/DataCleaning.smk --cluster-config ${CONFIG}/cluster_config.yml --configfile ${CONFIG}/config.yml --jobs 200 --printshellcmds --use-envmodules --cluster "sbatch -p {cluster.partition} --job-name {cluster.job-name} -e {cluster.error} -o {cluster.output}"

## RUN WITH CONDA ##
snakemake --snakefile ${PRGR}/DataCleaning.smk --cluster-config ${CONFIG}/cluster_config.yml --configfile ${CONFIG}/config.yml --jobs 200 --printshellcmds --use-conda --cluster "sbatch -p {cluster.partition} --job-name {cluster.job-name} -e {cluster.error} -o {cluster.output}"

## GENERATE PIPELINE REPORT ##
#snakemake --snakefile ${PRGR}/DataCleaning.smk --cluster-config ${CONFIG}/cluster_config.yml --configfile ${CONFIG}/config.yml --jobs 200 --printshellcmds --use-conda --report Snakemake_Report_DataCleaning.html --cluster "sbatch -p {cluster.partition} --job-name {cluster.job-name} -e {cluster.error} -o {cluster.output}"


                                  ### CLUSTER CONFIG MAKING USE OF PROFILES ###

## RUN WITH ENVMODULES ##

## RUN WITH CONDA ##

## GENERATE PIPELINE REPORT ##
