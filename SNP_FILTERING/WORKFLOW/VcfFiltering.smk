#!/usr/bin/env python

import pandas as pd
import os,sys
from itertools import compress

#singularity: "docker://condaforge/mambaforge"


####################   DEFINE CONFIG VARIABLES BASED ON CONFIG FILE   ####################

### Variables from config file
vcf_raw = config["VCF_file"]
MAX_percent_NA_per_INDIV = config["MAX_percent_NA_per_INDIV"]
outputs_dirname = config["OUTPUTS_DIRNAME"]

### Define paths
path_to_snakefile = workflow.snakefile
snakefile_dir = path_to_snakefile.rsplit('/', 1)[0]
scripts_dir = snakefile_dir+"/SCRIPTS"
working_directory = os.getcwd()

### Define outputs subfolders
outputs_directory = working_directory+"/"+outputs_dirname+"/VFC_FILTERING"
VCF_reports_dir = outputs_directory+"/REPORTS"

### PIPELINE ###

rule FinalTargets:
    input:
        #"__indivGenoFiltered_withPopStat.recode.vcf"
        "final_indivGeno_withPopStat.vcf"


 # ----------------------------------------------------------------------------------------------- #


rule Filter_genotypes:
    input:
        vcf_raw
    output:
        temp("__genoFiltered.recode.vcf")
    conda:
        "ENVS/conda_tools.yml"
    params:
        config["genotype_filter"]
    shell:
        "vcftools --gzvcf {input} {params} --recode --recode-INFO-all --out __genoFiltered"

rule Filter_individuals:
    input:
        "__genoFiltered.recode.vcf"
    output:
        "__indivGenoFiltered.recode.vcf",
        "__indiv_to_remove",
        "__out.imiss"
    conda:
        "ENVS/conda_tools.yml"
    params:
        config["MAX_percent_NA_per_INDIV"]
    shell:
        "vcftools --vcf  {input} --missing-indv;"
        "awk -v max_pc_NA={params} '{{if($5>max_pc_NA){{print $0}}}}' out.imiss | cut -f1 > __indiv_to_remove;"
        "mv out.imiss __out.imiss;"
        "vcftools --vcf {input} --remove  __indiv_to_remove --recode --recode-INFO-all --out __indivGenoFiltered" 

rule Add_PopGenStat:
    input:
        "__indivGenoFiltered.recode.vcf",
        scripts_dir+"/vcf_extra_info_header.txt"
    output:
        "__indivGenoFiltered_withPopStat.recode.vcf"
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "{scripts_dir}/add_popGenStat2VCF.sh {input}  __indivGenoFiltered_withPopStat.recode.vcf"

rule FinalFiltering:
    input:
        "__indivGenoFiltered_withPopStat.recode.vcf"
    output:
        "final_indivGeno_withPopStat.vcf"
    conda:
        "ENVS/conda_tools.yml"
    params:config["SNP_filter"]
    shell:
        "bgzip {input};"
        "tabix {input}.gz;"
        "bcftools filter -sFilterSmk -i '{params}' {input}.gz "
        "| bcftools view -f 'PASS' > {output};"

rule BuildReport:
    input:
        "final_indivGeno_withPopStat.vcf"
    output:
        "SNP_filtering_report.html"
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "bcftools stats {input}"