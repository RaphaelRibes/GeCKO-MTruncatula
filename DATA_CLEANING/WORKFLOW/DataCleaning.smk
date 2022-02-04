#!/usr/bin/env python

import pandas as pd
import os,sys


### Get variables from config file
samples = pd.read_csv(config["TAG_FILE"], sep='\t', header=None).iloc[:, 0]
scripts_dir = config["WORKFLOW_DIR"]+"/DATA_CLEANING/WORKFLOW/SCRIPTS"
working_directory = config["WORKDIR"] #essayer de le récupérer (ce serait là où le snakemake est lancé)
fastq_R1_raw = config["FASTQ_R1"]
fastq_R2_raw = config["FASTQ_R2"]
outputs_dirname = config["OUTPUTS_DIRNAME"]

### Create variables for outputs subfolders
outputs_directory = working_directory+"/"+outputs_dirname+"/DATA_CLEANING"
rawdata_reports_dir = outputs_directory+"/RAWDATA/REPORTS"
demult_dir = outputs_directory+"/DEMULT"
demult_reports_dir = outputs_directory+"/DEMULT/REPORTS"
demult_trim_dir = outputs_directory+"/DEMULT_TRIM"
demult_trim_reports_dir = outputs_directory+"/DEMULT_TRIM/REPORTS"
demult_trim_fastqc_reports_dir = demult_trim_reports_dir+"/fastqc"


### Extract folder path and base names from raw files
fastq_R1_raw_base = fastq_R1_raw.rsplit('/', 1)[1].replace('.fastq','').replace('.fq','').replace('.gz','')
fastq_R2_raw_base = fastq_R2_raw.rsplit('/', 1)[1].replace('.fastq','').replace('.fq','').replace('.gz','')
raw_data_dir = fastq_R1_raw.rsplit('/', 1)[0]

fastq_raw_base = fastq_R1_raw_base.replace('_R1','')



### PIPELINE ###

rule FinalTargets:
    input:
        # RAW DATA -> fastqc outputs -> pour mémoire, HS
        #rawdata_reports_dir+"/"+fastq_R1_raw_base+"_fastqc.zip",
        #rawdata_reports_dir+"/"+fastq_R2_raw_base+"_fastqc.zip",
        #
        # RAW DATA -> reads count output
        rawdata_reports_dir+"/Reads_Count_RawData.txt",

        # DEMULTIPLEXED FASTQ outputs -> pour mémoire, HS
        #expand("{working_directory}/DEMULT/{library_name}.R1.fastq.gz", library_name=library_names, working_directory=working_directory),
        #expand("{working_directory}/DEMULT/{library_name}.R2.fastq.gz", library_names=library_names, working_directory=working_directory),

        # DEMULTIPLEXED FASTQ -> reads count output
        demult_reports_dir+"/Reads_Count_Demult.txt",

        # TRIMMED FASTQ outputs -> pour mémoire, HS
        #expand("{working_directory}/DEMULT_TRIM/{library_name}_trimmed.R1.fastq.gz", library_name=library_names, working_directory=working_directory),
        #expand("{working_directory}/DEMULT_TRIM/{library_name}_trimmed.R2.fastq.gz", library_name=library_names, working_directory=working_directory),

        # TRIMMED FASTQ -> reads count output
        demult_trim_reports_dir+"/Reads_Count_DemultTrim.txt",

        # TRIMMED FASTQ -> multiQC output
        demult_trim_reports_dir+"/multiQC_Trimming_Report.html",

        # CONCATENATED TRIMMED FASTQ output -> pour mémoire, HS
        #working_directory+"/DEMULT_TRIM/smk_all_concat_trimmed.R1.fastq.gz",
        #working_directory+"/DEMULT_TRIM/smk_all_concat_trimmed.R2.fastq.gz"

        # CONCATENATED TRIMMED FASTQ -> fastqc outputs -> pour mémoire, HS
        #demult_trim_fastqc_reports_dir+"/all_samples_concat_trimmed.R1_fastqc.zip",
        #demult_trim_fastqc_reports_dir+"/all_samples_concat_trimmed.R2_fastqc.zip",

        # RAW DATA and CONCATENATED TRIMMED FASTQ -> multiQC output
        outputs_directory+"/multiQC_DataCleaning_Overall_Report.html"

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
        tag_file = config["TAG_FILE"]
    output:
        expand("{demult_dir}/{sample}.R1.fastq.gz", sample=samples, demult_dir=demult_dir),
        expand("{demult_dir}/{sample}.R2.fastq.gz", sample=samples, demult_dir=demult_dir)
    params:
        substitutions = config["DEMULT_SUBSTITUTIONS"],
        threads = config["DEMULT_THREADS"]
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "{scripts_dir}/demultiplex_with_cutadapt.sh --demultdir {demult_dir} --R1 {input.fastq_R1_raw} "
        "--R2 {input.fastq_R2_raw} --tag_file {input.tag_file} --nodes {params.threads} "
        "--substitutions {params.substitutions}"



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
        demult_trim_dir+"/trimming_cutadapt_{base}.info",
        temp(demult_trim_dir+"/tmp_trimming_cutadapt_{base}.info")
    params:
        threads = config["TRIMMING_THREADS"],
        qual = config["TRIMMING_QUAL"],
        min_length = config["TRIMMING_MIN_LENGTH"]
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "{scripts_dir}/trimming_with_cutadapt.sh --ind {wildcards.base} --trimdir {demult_trim_dir} "
        "--R1 {input.fastqs_R1_demult} --R2 {input.fastqs_R2_demult} --adapt_file {input.adapt_file} "
        "--nodes {params.threads} --qual {params.qual} --min_length {params.min_length};"
        "sed 's/R1//g' {demult_trim_dir}/trimming_cutadapt_{wildcards.base}.info | sed 's/R2//g' > {demult_trim_dir}/tmp_trimming_cutadapt_{wildcards.base}.info"


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
        expand("{demult_trim_dir}/tmp_trimming_cutadapt_{sample}.info", sample=samples, demult_trim_dir=demult_trim_dir),
        expand("{demult_trim_fastqc_reports_dir}/{sample}_trimmed.R1_fastqc.zip", sample=samples, demult_trim_fastqc_reports_dir=demult_trim_fastqc_reports_dir),
        expand("{demult_trim_fastqc_reports_dir}/{sample}_trimmed.R2_fastqc.zip", sample=samples, demult_trim_fastqc_reports_dir=demult_trim_fastqc_reports_dir)

    output:
        demult_trim_reports_dir+"/multiQC_Trimming_Report.html"
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "multiqc {input} -o {demult_trim_reports_dir} -n multiQC_Trimming_Report"



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
        rawdata_reports_dir+"/"+fastq_R1_raw_base+"_fastqc.zip",
        rawdata_reports_dir+"/"+fastq_R2_raw_base+"_fastqc.zip",
        demult_trim_fastqc_reports_dir+"/All_Samples_Concat_trimmed.R1_fastqc.zip",
        demult_trim_fastqc_reports_dir+"/All_Samples_Concat_trimmed.R2_fastqc.zip"
    output:
        outputs_directory+"/multiQC_DataCleaning_Overall_Report.html"
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "multiqc {input} -o {outputs_directory} -n multiQC_DataCleaning_Overall_Report"



#muse-login.hpc-lr.univ-montp2.fr
