#!/bin/bash


mkdir -p Logs_DataCleaningWorkflow

PRGR=/home/jogirodolle/save/CAPTURE_PIPELINES_SNAKEMAKE/DATA_CLEANING/WORKFLOW
CONFIG=/home/jogirodolle/work/DATA_TEST/VIR/CONFIG



## DRYRUN ##
#snakemake --snakefile ${PRGR}/DataCleaning.smk -n --dag --forceall --cluster-config ${CONFIG}/cluster_config_SGE.json --configfile ${CONFIG}/config_SGE.yml --printshellcmds

## GENERATE PIPELINE DIAGRAM ##
#snakemake --snakefile ${PRGR}/DataCleaning.smk -n --dag --forceall --cluster-config ${CONFIG}/cluster_config_SGE.json --configfile ${CONFIG}/mini_tmp/config_mini_SGE.yml | dot -Tsvg > data_cleaning_pipeline_mini.svg




                                        ### CLASSIC CLUSTER CONFIG ###

## RUN WITH CONDA ##
#snakemake --snakefile ${PRGR}/DataCleaning.smk --cluster-config ${CONFIG}/cluster_config_SGE.json --configfile ${CONFIG}/config_SGE.yml --jobs 200 --printshellcmds --use-conda --cluster "qsub -V -b yes -cwd -q {cluster.q} -N {cluster.N} -o {cluster.o} -e {cluster.e}"




                                  ### CLUSTER CONFIG MAKING USE OF PROFILES ###

## RUN WITH CONDA ##
snakemake --snakefile ${PRGR}/DataCleaning.smk --cluster-config ${CONFIG}/cluster_config_SGE.json --configfile ${CONFIG}/config_SGE.yml --jobs 200 --printshellcmds --use-conda --profile /home/jogirodolle/save/PROFILES_SNAKEMAKE/SGE

## RUN WITH CONDA + SINGULARITY ##
#snakemake --snakefile ${PRGR}/DataCleaning.smk --cluster-config ${CONFIG}/cluster_config_SGE.json --configfile ${CONFIG}/config_SGE.yml --jobs 200 --printshellcmds --use-singularity --use-conda --profile /home/jogirodolle/save/PROFILES_SNAKEMAKE/SGE

## GENERATE PIPELINE REPORT ##
#snakemake --snakefile ${PRGR}/DataCleaning.smk --cluster-config ${CONFIG}/cluster_config_SGE.json --configfile ${CONFIG}/config_SGE.yml --jobs 200 --printshellcmds --use-conda --report Snakemake_Report_DataCleaning.html --profile /home/jogirodolle/save/PROFILES_SNAKEMAKE/SGE






# Example on Migale :
#qsub -q long.q -V -b yes -cwd -N DataCleaning "conda activate snakemake-6.9.1 ; ./run_DataCleaning_workflow_Migale.sh"
#qsub -q short.q -V -b yes -cwd -N DataCleaning "conda activate snakemake-6.9.1 ; ./run_DataCleaning_workflow_Migale.sh"
