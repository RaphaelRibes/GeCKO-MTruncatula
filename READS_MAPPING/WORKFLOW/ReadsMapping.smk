#!/usr/bin/env python

#import pandas as pd
import os,sys
from itertools import compress

#singularity: "docker://condaforge/mambaforge"


####################   DEFINE CONFIG VARIABLES BASED ON CONFIG FILE   ####################

### Define paths
path_to_snakefile = workflow.snakefile
snakefile_dir = path_to_snakefile.rsplit('/', 1)[0]
scripts_dir = snakefile_dir+"/SCRIPTS"
working_directory = os.getcwd()

### Variables from config file
outputs_dirname = config["OUTPUTS_DIRNAME"]
paired_end = config["PAIRED_END"]
trim_dir = config["TRIM_DIR"]
if (len(trim_dir) == 0):
    trim_dir = working_directory+"/"+outputs_dirname+"/DATA_CLEANING/DEMULT_TRIM"

remove_dup = config["REMOVE_DUP"]
if remove_dup:
    rm_dup = True
else:
    rm_dup = False

ref_name = config["REFERENCE_NAME"]

### Samples list
if paired_end:
    fastqs_R1_list = [fastq for fastq in os.listdir(trim_dir) if '.R1.fastq.gz' in fastq]

    fastqs_R2_list = []
    for fastq_R1 in fastqs_R1_list:
        fastq_R2 = fastq_R1.replace(".R1.", ".R2.")
        fastqs_R2_list.append(fastq_R2)

    samples = []
    for fastq_R1 in fastqs_R1_list:
        sample = fastq_R1.replace(".R1.fastq.gz", "")
        samples.append(sample)

else:
    fastqs_list = [fastq for fastq in os.listdir(trim_dir) if '.fastq.gz' in fastq]

    samples = []
    for fastq in fastqs_list:
        sample = fastq.replace(".fastq.gz", "")
        samples.append(sample)


### Define outputs subfolders
outputs_directory = working_directory+"/"+outputs_dirname+"/READS_MAPPING"
mapping_dir = outputs_directory+"/MAPPING_TO_"+ref_name
bams_dir = mapping_dir+"/BAMS"
bams_reports_dir = bams_dir+"/REPORTS"



### FUNCTIONS

def buildExpectedFiles(filesNames, isExpected):
    expectedFiles = list(compress(filesNames, isExpected))
    return(expectedFiles)


### PIPELINE ###

rule FinalTargets:
    input:
        #expand("{bams_reports_dir}/flagstat_{sample}", sample=samples, bams_reports_dir=bams_reports_dir),
        #bams_reports_dir+"/reads_histogram.pdf",
        bams_reports_dir+"/reads_histograms.pdf"

 # ----------------------------------------------------------------------------------------------- #

# ajouter une règle optionnelle qui indexe la ref si besoin (il faut qu'un fichier .bwt.2bit.64 existe)


rule Mapping_PairedEndFastqs:
    input:
        fastq_paired_R1 = trim_dir+"/{base}.R1.fastq.gz",
        fastq_paired_R2 = trim_dir+"/{base}.R2.fastq.gz",
        ref = config["REFERENCE"]
    output:
        buildExpectedFiles([bams_dir+"/{base}.bam"],[paired_end])
    conda:
        "ENVS/conda_tools.yml"
    params:
        mapper = config["MAPPER"],
        extra_mapper_params = config["EXTRA_MAPPER_PARAMS"],
        technology = config["SEQUENCING_TECHNOLOGY"],
        markdup_params = config["MARKDUP_PARAMS"],
        index_params = config["INDEX_PARAMS"]
    shell:
        "{scripts_dir}/mapping.sh --paired_end --fastq_R1 \"{input.fastq_paired_R1}\" --fastq_R2 \"{input.fastq_paired_R2}\" "
        "--ref {input.ref} --mapper {params.mapper} --mapper_params \"{params.extra_mapper_params}\" --technology \"{params.technology}\" "
        "--output_dir {bams_dir} --reports_dir {bams_reports_dir} --sample {wildcards.base} --rm_dup {rm_dup} --markdup_params \"{params.markdup_params}\" "
        "--index_params \"{params.index_params}\""


rule Mapping_SingleEndFastqs:
    input:
        fastq_single = trim_dir+"/{base}.fastq.gz",
        ref = config["REFERENCE"]
    output:
        buildExpectedFiles([bams_dir+"/{base}.bam"],[not paired_end])
    conda:
        "ENVS/conda_tools.yml"
    params:
        mapper = config["MAPPER"],
        extra_mapper_params = config["EXTRA_MAPPER_PARAMS"],
        technology = config["SEQUENCING_TECHNOLOGY"],
        markdup_params = config["MARKDUP_PARAMS"],
        index_params = config["INDEX_PARAMS"]
    shell:
        "{scripts_dir}/mapping.sh --single_end --fastq \"{input.fastq_single}\" "
        "--ref {input.ref} --mapper {params.mapper} --mapper_params \"{params.extra_mapper_params}\" --technology \"{params.technology}\" "
        "--output_dir {bams_dir} --reports_dir {bams_reports_dir} --sample {wildcards.base} --rm_dup {rm_dup} --markdup_params \"{params.markdup_params}\" "
        "--index_params \"{params.index_params}\""


rule CountReads_Bams:
    input:
        bams_dir+"/{base}.bam"
    output:
        bams_reports_dir+"/flagstat_{base}"
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "samtools flagstat {input} > {output}"


rule Summarize_PairedEndReadsCount:
    input:
        expand("{bams_reports_dir}/flagstat_{sample}", sample=samples, bams_reports_dir=bams_reports_dir)
    output:
        buildExpectedFiles([bams_reports_dir+"/summary_flagstat"],[paired_end])
    shell:
        "{scripts_dir}/summarize_flagstat.sh --pairedEnd --flagstat_folder {bams_reports_dir}"


rule Summarize_SingleEndReadsCount:
    input:
        expand("{bams_reports_dir}/flagstat_{sample}", sample=samples, bams_reports_dir=bams_reports_dir)
    output:
        buildExpectedFiles([bams_reports_dir+"/summary_flagstat"],[not paired_end])
    shell:
        "{scripts_dir}/summarize_flagstat.sh --singleEnd --flagstat_folder {bams_reports_dir}"



rule Histogram_ReadsCount:
    input:
        bams_reports_dir+"/summary_flagstat"
    output:
        bams_reports_dir+"/reads_histograms.pdf"
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "python3 {scripts_dir}/plot_reads_histogram.py {bams_reports_dir}"



#front.migale.inrae.fr
#muse-login.hpc-lr.univ-montp2.fr
