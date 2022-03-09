#!/bin/bash
#
#SBATCH -J VcfFilteringSmk
#SBATCH -o VcfFilteringSmk."%j".out
#SBATCH -e VcfFilteringSmk."%j".err

#SBATCH -p agap_normal

module purge
module load snakemake/5.13.0
module load anaconda/python3.8
module load singularity/3.6.3

mkdir -p Logs_VcfFilteringWorkflow

PRGR=/storage/replicated/cirad_users/ardissonm/WORKFLOW/SNP_CALLING
CONFIG=./CONFIG
PROFILE=${PRGR}"/PROFILES"


## DRYRUN ##
#snakemake --snakefile ${PRGR}/DataCleaning_PairedEnd.smk -n --dag --forceall --cluster-config ${CONFIG}/cluster_config_Slurm.json --configfile ${CONFIG}/config_file_PairedEnd.yml --printshellcmds

## GENERATE PIPELINE DIAGRAM -> à lancer hors sbatch, dot est indisponible sur les noeuds de calcul ##
#snakemake --snakefile ${PRGR}/DataCleaning_PairedEnd.smk -n --dag --forceall --cluster-config ${CONFIG}/cluster_config_Slurm.json --configfile ${CONFIG}/mini_tmp/config_mini_Slurm.yml | dot -Tsvg > data_cleaning_pipeline_mini.svg



                                        ### CLASSIC CLUSTER CONFIG ###

## RUN WITH ENVMODULES ##
#snakemake --snakefile ${PRGR}/DataCleaning_PairedEnd.smk --cluster-config ${CONFIG}/cluster_config_Slurm.json --configfile ${CONFIG}/config_file_PairedEnd.yml --jobs 200 --printshellcmds --use-envmodules --cluster "sbatch -p {cluster.partition} --job-name {cluster.job-name} -e {cluster.error} -o {cluster.output}"

## RUN WITH CONDA ##
#snakemake --snakefile ${PRGR}/DataCleaning_PairedEnd.smk --cluster-config ${CONFIG}/cluster_config_Slurm.json --configfile ${CONFIG}/config_file_PairedEnd.yml --jobs 200 --printshellcmds --use-conda --cluster "sbatch -p {cluster.partition} --job-name {cluster.job-name} -e {cluster.error} -o {cluster.output}"

## GENERATE PIPELINE REPORT -> à lancer hors sbatch sinon ça ne marche pas ##
#snakemake --snakefile ${PRGR}/DataCleaning_PairedEnd.smk --cluster-config ${CONFIG}/cluster_config_Slurm.json --configfile ${CONFIG}/config_file_PairedEnd.yml --jobs 200 --printshellcmds --use-conda --report Snakemake_Report_DataCleaning.html --cluster "sbatch -p {cluster.partition} --job-name {cluster.job-name} -e {cluster.error} -o {cluster.output}"


                                  ### CLUSTER CONFIG MAKING USE OF PROFILES ###

## RUN WITH CONDA ##
snakemake  --snakefile ${PRGR}/VcfFiltering.smk --cluster-config ${CONFIG}/cluster_config_Slurm_VcfFiltering.json --configfile ${CONFIG}/config_file_VcfFiltering.yml --jobs 200 --printshellcmds --use-conda --profile ${PROFILE}/SLURM

## RUN WITH CONDA + SINGULARITY -> pour l'instant ça ne marche pas ##
#snakemake --snakefile ${PRGR}/DataCleaning_PairedEnd.smk --cluster-config ${CONFIG}/cluster_config_Slurm.json --configfile ${CONFIG}/config_file_PairedEnd.yml --jobs 200 --printshellcmds --use-conda --use-singularity --singularity-args "-B /work_home/jogirodolle:/home/jogirodolle/work,/save_home/jogirodolle:/home/jogirodolle/save" --profile ${PROFILE}/SLURM

## GENERATE PIPELINE REPORT -> à lancer hors sbatch sinon ça ne marche pas ##
#snakemake --snakefile ${PRGR}/DataCleaning_PairedEnd.smk --cluster-config ${CONFIG}/cluster_config_Slurm.json --configfile ${CONFIG}/config_file_PairedEnd.yml --jobs 200 --printshellcmds --use-conda --report Snakemake_Report_DataCleaning_PairedEnd.html --profile ${PROFILE}/SLURM

## CONTAINERIZE -> option inconnue pour l'instant ##
#snakemake --snakefile ${PRGR}/DataCleaning_PairedEnd.smk --cluster-config ${CONFIG}/cluster_config_Slurm.json --configfile ${CONFIG}/config_file_PairedEnd.yml --jobs 200 --printshellcmds --use-conda --profile ${PROFILE}/SLURM --containerize > Dockerfile
