#!/usr/bin/env python

import pandas as pd
import os,sys
from itertools import compress

#singularity: "docker://condaforge/mambaforge"


####################   DEFINE CONFIG VARIABLES BASED ON CONFIG FILE   ####################

### Variables from config file
samples = pd.read_csv(config["ADAPT_FILE"], sep='\t', header=None).iloc[:, 0]
fastq_R1_raw = config["FASTQ_R1"]
fastq_R2_raw = config["FASTQ_R2"]
outputs_dirname = config["OUTPUTS_DIRNAME"]
user_demult_dir = config["DEMULT_DIR"]

if (len(user_demult_dir) == 0):
    performDemultiplexing = True
else:
    performDemultiplexing = False

### Raw fastq files path and base names
fastq_R1_raw_base = fastq_R2_raw_base = raw_data_dir = fastq_raw_base = ""
if performDemultiplexing:
    fastq_R1_raw_base = fastq_R1_raw.rsplit('/', 1)[1].replace('.fastq','').replace('.fq','').replace('.gz','')
    fastq_R2_raw_base = fastq_R2_raw.rsplit('/', 1)[1].replace('.fastq','').replace('.fq','').replace('.gz','')
    raw_data_dir = fastq_R1_raw.rsplit('/', 1)[0]
    fastq_raw_base = fastq_R1_raw_base.replace('_R1','')

### Define paths
path_to_snakefile = workflow.snakefile
snakefile_dir = path_to_snakefile.rsplit('/', 1)[0]
scripts_dir = snakefile_dir+"/SCRIPTS"
working_directory = os.getcwd()

### Define outputs subfolders
outputs_directory = working_directory+"/"+outputs_dirname+"/DATA_CLEANING"
rawdata_reports_dir = outputs_directory+"/RAWDATA/REPORTS"

if performDemultiplexing:
    demult_dir = outputs_directory+"/DEMULT"
    demult_cutadapt_reports_dir = demult_dir+"/REPORTS/CUTADAPT_INFOS"
else:
    demult_dir = user_demult_dir
    demult_cutadapt_reports_dir = ""

demult_reports_dir = outputs_directory+"/DEMULT/REPORTS"
demult_trim_dir = outputs_directory+"/DEMULT_TRIM"
demult_trim_reports_dir = outputs_directory+"/DEMULT_TRIM/REPORTS"
demult_trim_cutadapt_reports_dir = demult_trim_reports_dir+"/CUTADAPT_INFOS"
demult_trim_fastqc_reports_dir = demult_trim_reports_dir+"/FASTQC"


### FUNCTIONS

def buildExpectedFiles(filesNames, isExpected):
    expectedFiles = list(compress(filesNames, isExpected))
    return(expectedFiles)


### PIPELINE ###

rule FinalTargets:
    input:
        buildExpectedFiles(
        [rawdata_reports_dir+"/Reads_Count_RawData.txt",
        demult_reports_dir+"/Reads_Count_Demult.txt",
        demult_trim_reports_dir+"/Reads_Count_DemultTrim.txt",
        demult_trim_reports_dir+"/multiQC_Trimming_Report.html",
        outputs_directory+"/multiQC_DataCleaning_Overall_Report.html"],

        [performDemultiplexing, True, True, True, True]
        )


 # ----------------------------------------------------------------------------------------------- #


rule Fastqc_RawFastqs:
    input:
        raw_data_dir+"/{base}.fastq.gz"
    output:
        rawdata_reports_dir+"/{base}_fastqc.zip"
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "fastqc -o {rawdata_reports_dir} {input}"


rule CountReads_RawFastqs:
    input:
        fastq_R1_raw,
        fastq_R2_raw
    output:
        rawdata_reports_dir+"/Reads_Count_RawData.txt"
    shell:
        "{scripts_dir}/fastq_read_count.sh --fastq_dir {raw_data_dir} --output {output}"


rule Demultiplex_RawFastqs:
    input:
        fastq_R1_raw = fastq_R1_raw,
        fastq_R2_raw = fastq_R2_raw,
        barcode_file = config["BARCODE_FILE"]
    output:
        expand("{demult_dir}/{sample}.R1.fastq.gz", sample=samples, demult_dir=demult_dir),
        expand("{demult_dir}/{sample}.R2.fastq.gz", sample=samples, demult_dir=demult_dir),
        demult_cutadapt_reports_dir+"/demultiplexing_cutadapt.info"
    params:
        substitutions = config["DEMULT_SUBSTITUTIONS"],
        threads = config["DEMULT_THREADS"]
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "{scripts_dir}/demultiplex_with_cutadapt_PE.sh --demultdir {demult_dir} --R1 {input.fastq_R1_raw} "
        "--R2 {input.fastq_R2_raw} --barcode_file {input.barcode_file} --nodes {params.threads} "
        "--substitutions {params.substitutions};"
        "mv {demult_dir}/demultiplexing_cutadapt.info {demult_cutadapt_reports_dir}"

rule CountReads_DemultFastqs:
    input:
        expand("{demult_dir}/{sample}.R1.fastq.gz", sample=samples, demult_dir=demult_dir),
        expand("{demult_dir}/{sample}.R2.fastq.gz", sample=samples, demult_dir=demult_dir)
    output:
        demult_reports_dir+"/Reads_Count_Demult.txt"
    shell:
        "{scripts_dir}/fastq_read_count.sh --fastq_dir {demult_dir} --output {output}"


rule Trimming_DemultFastqs:
    input:
        fastqs_R1_demult = demult_dir+"/{base}.R1.fastq.gz",
        fastqs_R2_demult = demult_dir+"/{base}.R2.fastq.gz",
        adapt_file = config["ADAPT_FILE"]
    output:
        demult_trim_dir+"/{base}_trimmed.R1.fastq.gz",
        demult_trim_dir+"/{base}_trimmed.R2.fastq.gz",
        demult_trim_cutadapt_reports_dir+"/trimming_cutadapt_{base}.info",
        temp(demult_trim_cutadapt_reports_dir+"/tmp_trimming_cutadapt_{base}.info")
    params:
        threads = config["TRIMMING_THREADS"],
        quality_cutoff = config["TRIMMING_QUAL"],
        minimum_length = config["TRIMMING_MIN_LENGTH"]
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "{scripts_dir}/trimming_with_cutadapt_PE.sh --sample {wildcards.base} --trimdir {demult_trim_dir} "
        "--R1 {input.fastqs_R1_demult} --R2 {input.fastqs_R2_demult} --adapt_file {input.adapt_file} "
        "--nodes {params.threads} --quality_cutoff {params.quality_cutoff} --minimum_length {params.minimum_length};"
        "mv {demult_trim_dir}/trimming_cutadapt_{wildcards.base}.info {demult_trim_cutadapt_reports_dir}/trimming_cutadapt_{wildcards.base}.info;"
        "sed 's/R1//g' {demult_trim_cutadapt_reports_dir}/trimming_cutadapt_{wildcards.base}.info | sed 's/R2//g' > {demult_trim_cutadapt_reports_dir}/tmp_trimming_cutadapt_{wildcards.base}.info"


rule CountReads_TrimmedFastqs:
    input:
        expand("{demult_trim_dir}/{sample}_trimmed.R1.fastq.gz", sample=samples, demult_trim_dir=demult_trim_dir),
        expand("{demult_trim_dir}/{sample}_trimmed.R2.fastq.gz", sample=samples, demult_trim_dir=demult_trim_dir)
    output:
        demult_trim_reports_dir+"/Reads_Count_DemultTrim.txt"
    shell:
        "{scripts_dir}/fastq_read_count.sh --fastq_dir {demult_trim_dir} --output {output}"


rule Fastqc_TrimmedFastqs:
    input:
        demult_trim_dir+"/{base}_trimmed.{R}.fastq.gz"
    output:
        demult_trim_fastqc_reports_dir+"/{base}_trimmed.{R}_fastqc.zip"
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "fastqc -o {demult_trim_fastqc_reports_dir} {input}"


rule MultiQC_TrimmedFastqs:
    input:
        expand("{demult_trim_cutadapt_reports_dir}/tmp_trimming_cutadapt_{sample}.info", sample=samples, demult_trim_cutadapt_reports_dir=demult_trim_cutadapt_reports_dir),
        expand("{demult_trim_fastqc_reports_dir}/{sample}_trimmed.R1_fastqc.zip", sample=samples, demult_trim_fastqc_reports_dir=demult_trim_fastqc_reports_dir),
        expand("{demult_trim_fastqc_reports_dir}/{sample}_trimmed.R2_fastqc.zip", sample=samples, demult_trim_fastqc_reports_dir=demult_trim_fastqc_reports_dir)

    output:
        demult_trim_reports_dir+"/multiQC_Trimming_Report.html"
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "multiqc {input} -o {demult_trim_reports_dir} -n multiQC_Trimming_Report;"
        "sed -i -e '/header_mqc-generalstats-cutadapt-percent_trimmed/s/>\\([a-zA-Z][^>]*\\)_trimmed.R/>\\1.R/g' {demult_trim_reports_dir}/multiQC_Trimming_Report.html"


rule Concatenate_TrimmedFastqs:
    input:
        fastqs_trim = expand("{demult_trim_dir}/{sample}_trimmed.{{R}}.fastq.gz", sample=samples, demult_trim_dir=demult_trim_dir)
    output:
        temp(demult_trim_dir+"/"+fastq_raw_base+"_trimmed.{R}.fastq.gz")
    shell:
        "cat {input.fastqs_trim} > {output}"


rule Fastqc_ConcatTrimmedFastqs:
    input:
        demult_trim_dir+"/"+fastq_raw_base+"_trimmed.{R}.fastq.gz"
    output:
        demult_trim_fastqc_reports_dir+"/All_Samples_Concat_trimmed.{R}_fastqc.zip"
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "fastqc -o {demult_trim_fastqc_reports_dir} {input} ;"
        "mv {demult_trim_fastqc_reports_dir}/{fastq_raw_base}_trimmed.{wildcards.R}_fastqc.html {demult_trim_fastqc_reports_dir}/All_Samples_Concat_trimmed.{wildcards.R}_fastqc.html ;"
        "mv {demult_trim_fastqc_reports_dir}/{fastq_raw_base}_trimmed.{wildcards.R}_fastqc.zip {demult_trim_fastqc_reports_dir}/All_Samples_Concat_trimmed.{wildcards.R}_fastqc.zip"


rule MultiQC_Global:
    input:
        buildExpectedFiles(
        [rawdata_reports_dir+"/"+fastq_R1_raw_base+"_fastqc.zip",
        rawdata_reports_dir+"/"+fastq_R2_raw_base+"_fastqc.zip",
        demult_trim_fastqc_reports_dir+"/All_Samples_Concat_trimmed.R1_fastqc.zip",
        demult_trim_fastqc_reports_dir+"/All_Samples_Concat_trimmed.R2_fastqc.zip"],

        [performDemultiplexing, performDemultiplexing, True, True]
        )

    output:
        outputs_directory+"/multiQC_DataCleaning_Overall_Report.html"
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "multiqc {input} -o {outputs_directory} -n multiQC_DataCleaning_Overall_Report"



#front.migale.inrae.fr
#muse-login.hpc-lr.univ-montp2.fr
