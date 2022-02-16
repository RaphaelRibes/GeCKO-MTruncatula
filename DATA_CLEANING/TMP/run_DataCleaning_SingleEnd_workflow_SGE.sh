#!/bin/bash


mkdir -p Logs_DataCleaningWorkflow

PRGR=/home/jogirodolle/save/CAPTURE_PIPELINES_SNAKEMAKE/DATA_CLEANING/WORKFLOW
CONFIG=/home/jogirodolle/work/DATA_TEST/VIR/CONFIG
PROFILE=${PRGR}"/PROFILES"

## DRYRUN ##
#snakemake --snakefile ${PRGR}/DataCleaning_SingleEnd.smk -n --dag --forceall --cluster-config ${CONFIG}/cluster_config_SGE.json --configfile ${CONFIG}/config_file_SingleEnd.yml --printshellcmds

## GENERATE PIPELINE DIAGRAM ##
#snakemake --snakefile ${PRGR}/DataCleaning_SingleEnd.smk -n --dag --forceall --cluster-config ${CONFIG}/cluster_config_SGE.json --configfile ${CONFIG}/mini_tmp/config_mini_SGE.yml | dot -Tsvg > data_cleaning_pipeline_mini.svg




                                        ### CLASSIC CLUSTER CONFIG ###

## RUN WITH CONDA ##
#snakemake --snakefile ${PRGR}/DataCleaning_SingleEnd.smk --cluster-config ${CONFIG}/cluster_config_SGE.json --configfile ${CONFIG}/config_file_SingleEnd.yml --jobs 200 --printshellcmds --use-conda --cluster "qsub -V -b yes -cwd -q {cluster.q} -N {cluster.N} -o {cluster.o} -e {cluster.e}"




                                  ### CLUSTER CONFIG MAKING USE OF PROFILES ###

## RUN WITH CONDA ##
  # Create environment only (once)
#snakemake --snakefile ${PRGR}/DataCleaning_SingleEnd.smk --cluster-config ${CONFIG}/cluster_config_SGE.json --configfile ${CONFIG}/config_file_SingleEnd.yml --jobs 200 --printshellcmds --use-conda --conda-create-envs-only --profile ${PROFILE}/SGE

  # Run
snakemake --snakefile ${PRGR}/DataCleaning_SingleEnd.smk --cluster-config ${CONFIG}/cluster_config_SGE.json --configfile ${CONFIG}/config_file_SingleEnd.yml --jobs 200 --printshellcmds --use-conda --profile ${PROFILE}/SGE


## RUN WITH CONDA + SINGULARITY (penser à bien mettre 'singularity: "docker://condaforge/mambaforge"' dans le .smk) ##
#snakemake --snakefile ${PRGR}/DataCleaning_SingleEnd.smk --cluster-config ${CONFIG}/cluster_config_SGE.json --configfile ${CONFIG}/config_file_SingleEnd.yml --jobs 200 --printshellcmds --use-singularity --singularity-args "-B /work_home/jogirodolle:/home/jogirodolle/work,/save_home/jogirodolle:/home/jogirodolle/save" --use-conda --profile ${PROFILE}/SGE

## GENERATE PIPELINE REPORT ##
#snakemake --snakefile ${PRGR}/DataCleaning_SingleEnd.smk --cluster-config ${CONFIG}/cluster_config_SGE.json --configfile ${CONFIG}/config_file_SingleEnd.yml --jobs 200 --printshellcmds --use-conda --report Snakemake_Report_DataCleaning_SingleEnd.html --profile ${PROFILE}/SGE

## CONTAINERIZE -> marche pas pour l'instant, il faudrait passer à la dernière version de snakemake je pense ##
#snakemake --snakefile ${PRGR}/DataCleaning_SingleEnd.smk --cluster-config ${CONFIG}/cluster_config_SGE.json --configfile ${CONFIG}/config_file_SingleEnd.yml --jobs 200 --printshellcmds --use-conda --profile ${PROFILE}/SGE --containerize > Dockerfile




# Example on Migale :
#qsub -q long.q -V -b yes -cwd -N DataCleaning "conda activate snakemake-6.9.1 ; ./run_DataCleaning_SingleEnd_workflow_SGE.sh"
#qsub -q short.q -V -b yes -cwd -N DataCleaning "conda activate snakemake-6.9.1 ; ./run_DataCleaning_SingleEnd_workflow_SGE.sh"
