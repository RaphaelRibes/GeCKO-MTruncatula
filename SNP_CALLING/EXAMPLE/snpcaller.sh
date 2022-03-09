

module purge
module load snakemake/5.13.0
module load anaconda/python3.8

snakemake --snakefile /home/ardissonm/replicated/WORKFLOW/SNP_CALLING/GATK_ENV.smk --use-conda --jobs 1 
#--configfile CONFIG/config_SnpCaller.yml