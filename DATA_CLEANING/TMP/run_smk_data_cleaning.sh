#!/bin/bash
#
#SBATCH -J DataCleaningSmk
#SBATCH -o DataCleaningSmk."%j".out
#SBATCH -e DataCleaningSmk."%j".err

#SBATCH -p agap_normal

module purge
module load snakemake/5.13.0
module load anaconda/python3.8

mkdir -p pipeline_rules_logs





## DRYRUN ##
#snakemake --snakefile data_cleaning_pipeline.smk -n -p --dag --forceall --cluster-config cluster_config.yml --configfile config.yml #--report report-test.html

## GENERATE PIPELINE DIAGRAM ##
#snakemake --snakefile data_cleaning_pipeline.smk -n --dag --forceall --cluster-config cluster_config.yml --configfile config.yml | dot -Tsvg > data_cleaning_pipeline.svg
#snakemake --snakefile data_cleaning_pipeline.smk -n --dag --forceall --cluster-config cluster_config.yml --configfile config_mini.yml | dot -Tsvg > data_cleaning_pipeline_mini.svg

## RUN ##
#snakemake --snakefile data_cleaning_pipeline.smk --cluster-config cluster_config.yml --configfile config.yml --jobs 200 --printshellcmds --use-envmodules --cluster "sbatch -p {cluster.partition} --job-name {cluster.job-name} -e {cluster.error} -o {cluster.output}"
snakemake --snakefile data_cleaning_pipeline.smk --cluster-config cluster_config.yml --configfile config.yml --jobs 200 --printshellcmds --use-conda --cluster "sbatch -p {cluster.partition} --job-name {cluster.job-name} -e {cluster.error} -o {cluster.output}"
