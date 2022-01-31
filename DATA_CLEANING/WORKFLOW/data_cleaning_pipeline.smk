#!/usr/bin/env python

import pandas as pd
import os,sys


### Get variables from config file
library_names = pd.read_csv(config["TAG_FILE"], sep='\t', header=None).iloc[:, 0]
prgr = config["PRGR"]
working_directory = config["WORKDIR"]
fastq_R1_raw = config["FASTQ_R1"]
fastq_R2_raw = config["FASTQ_R2"]

### Extract folder path and base names from raw files
fastq_R1_raw_base = fastq_R1_raw.rsplit('/', 1)[1].replace('.fastq','').replace('.fq','').replace('.gz','')
fastq_R2_raw_base = fastq_R2_raw.rsplit('/', 1)[1].replace('.fastq','').replace('.fq','').replace('.gz','')
raw_data_folder = fastq_R1_raw.rsplit('/', 1)[0]



### PIPELINE ###

rule FinalTargets:
    input:
        # RAW DATA -> fastqc outputs
        working_directory+"/RESULTS/RAWDATA/"+fastq_R1_raw_base+"_fastqc.html",
        working_directory+"/RESULTS/RAWDATA/"+fastq_R2_raw_base+"_fastqc.html",

        # RAW DATA -> reads count output

        # DEMULTIPLEXED FASTQ outputs
        #expand("{working_directory}/DEMULT/{library_name}.R1.fastq.gz", library_name=library_names, working_directory=working_directory),
        #expand("{working_directory}/DEMULT/{library_name}.R2.fastq.gz", library_names=library_names, working_directory=working_directory),

        # DEMULTIPLEXED FASTQ -> reads count output
        working_directory+"/RESULTS/DEMULT/fastq_read_count_demult.txt",

        # TRIMMED FASTQ outputs
        #expand("{working_directory}/DEMULT_TRIM/{library_name}_trimmed.R1.fastq.gz", library_name=library_names, working_directory=working_directory),
        #expand("{working_directory}/DEMULT_TRIM/{library_name}_trimmed.R2.fastq.gz", library_name=library_names, working_directory=working_directory),

        # TRIMMED FASTQ -> reads count output
        working_directory+"/RESULTS/DEMULT_TRIM/fastq_read_count_demult_trim.txt",

        # TRIMMED FASTQ -> multiQC output
        working_directory+"/RESULTS/DEMULT_TRIM/trimming_multiqc_report.html",

        # CONCATENATED TRIMMED FASTQ output
        #working_directory+"/DEMULT_TRIM/smk_all_concat_trimmed.R1.fastq.gz",
        #working_directory+"/DEMULT_TRIM/smk_all_concat_trimmed.R2.fastq.gz"

        # CONCATENATED TRIMMED FASTQ -> fastqc outputs
        working_directory+"/RESULTS/DEMULT_TRIM/all_trimmed.R1_fastqc.html",
        working_directory+"/RESULTS/DEMULT_TRIM/all_trimmed.R2_fastqc.html"


 # ----------------------------------------------------------------------------------------------- #


rule Fastqc_RawFastqs:
    input:
        raw_data_folder+"/{base}.fastq.gz"
    output:
        working_directory+"/RESULTS/RAWDATA/{base}_fastqc.html"
    envmodules:
        config["modules"]["fastqc"]
    shell:
        "fastqc -o {working_directory}/RESULTS/RAWDATA/ {input}"



rule Demultiplex_RawFastqs:
    input:
        fastq_R1_raw = fastq_R1_raw,
        fastq_R2_raw = fastq_R2_raw,
        tag_file = config["TAG_FILE"]
    output:
        expand("{working_directory}/DEMULT/{library_name}.R1.fastq.gz", library_name=library_names, working_directory=working_directory),
        expand("{working_directory}/DEMULT/{library_name}.R2.fastq.gz", library_name=library_names, working_directory=working_directory)
    params:
        auth_subst_fraction = config["AUTH_SUBST_FRACTION"],
        threads = config["DEMULT_THREADS"]
    envmodules:
        config["modules"]["cutadapt"]
    shell:
        "{prgr}/demultiplex_with_cutadapt.sh --demultdir {working_directory}/DEMULT --R1 {input.fastq_R1_raw} "
        "--R2 {input.fastq_R2_raw} --tag_file {input.tag_file} --nodes {params.threads} "
        "--auth_subst {params.auth_subst_fraction}"



rule CountReads_DemultFastqs:
    input:
        fastqs_R1_demult = expand("{working_directory}/DEMULT/{library_name}.R1.fastq.gz", library_name=library_names, working_directory=working_directory),
        fastqs_R2_demult = expand("{working_directory}/DEMULT/{library_name}.R2.fastq.gz", library_name=library_names, working_directory=working_directory)
    output:
        working_directory+"/RESULTS/DEMULT/fastq_read_count_demult.txt"
    shell:
        "{prgr}/fastq_read_count.sh --fastq_dir {working_directory}/DEMULT --output {output}"



rule Trimming_DemultFastqs:
    input:
        fastqs_R1_demult = working_directory+"/DEMULT/{base}.R1.fastq.gz",
        fastqs_R2_demult = working_directory+"/DEMULT/{base}.R2.fastq.gz",
        adapt_file = config["ADAPT_FILE"]
    output:
        working_directory+"/DEMULT_TRIM/{base}_trimmed.R1.fastq.gz",
        working_directory+"/DEMULT_TRIM/{base}_trimmed.R2.fastq.gz",
        working_directory+"/DEMULT_TRIM/trimming_cutadapt_report_{base}.txt"
    params:
        threads = config["TRIMMING_THREADS"],
        qual = config["TRIMMING_QUAL"],
        min_length = config["TRIMMING_MIN_LENGTH"]
    envmodules:
        config["modules"]["cutadapt"]
    shell:
        "{prgr}/trimming_with_cutadapt.sh --ind {wildcards.base} --trimdir {working_directory}/DEMULT_TRIM "
        "--R1 {input.fastqs_R1_demult} --R2 {input.fastqs_R2_demult} --adapt_file {input.adapt_file} "
        "--nodes {params.threads} --qual {params.qual} --min_length {params.min_length}"


rule CountReads_TrimmedFastqs:
    input:
        fastqs_R1_trim = expand("{working_directory}/DEMULT_TRIM/{library_name}_trimmed.R1.fastq.gz", library_name=library_names, working_directory=working_directory),
        fastqs_R2_trim = expand("{working_directory}/DEMULT_TRIM/{library_name}_trimmed.R1.fastq.gz", library_name=library_names, working_directory=working_directory)
    output:
        working_directory+"/RESULTS/DEMULT_TRIM/fastq_read_count_demult_trim.txt"
    shell:
        "{prgr}/fastq_read_count.sh --fastq_dir {working_directory}/DEMULT_TRIM --output {output}"


rule MultiQC_TrimmedFastqs:
    input:
        expand("{working_directory}/DEMULT_TRIM/trimming_cutadapt_report_{library_name}.txt", library_name=library_names, working_directory=working_directory)
    output:
        working_directory+"/RESULTS/DEMULT_TRIM/trimming_multiqc_report.html"
    envmodules:
        config["modules"]["multiqc"]
    shell:
        "multiqc {input} -o {working_directory}/RESULTS/DEMULT_TRIM -n trimming_multiqc_report"


rule Concatenate_TrimmedFastqs:
    input:
        fastqs_trim = expand("{working_directory}/DEMULT_TRIM/{library_name}_trimmed.{{R}}.fastq.gz", library_name=library_names, working_directory=working_directory)
    output:
        temp(working_directory+"/DEMULT_TRIM/tmp_all_concat_trimmed.{R}.fastq.gz")
    shell:
        "cat {input.fastqs_trim} > {working_directory}/DEMULT_TRIM/tmp_all_concat_trimmed.{wildcards.R}.fastq.gz"


rule Fastqc_ConcatTrimmedFastqs:
    input:
        working_directory+"/DEMULT_TRIM/tmp_all_concat_trimmed.{R}.fastq.gz"
    output:
        working_directory+"/RESULTS/DEMULT_TRIM/all_trimmed.{R}_fastqc.html"
    envmodules:
        config["modules"]["fastqc"]
    shell:
        "fastqc -o {working_directory}/RESULTS/DEMULT_TRIM/ {input} ;"
        "mv {working_directory}/RESULTS/DEMULT_TRIM/tmp_all_concat_trimmed.{wildcards.R}_fastqc.html {working_directory}/RESULTS/DEMULT_TRIM/all_trimmed.{wildcards.R}_fastqc.html ;"
        "mv {working_directory}/RESULTS/DEMULT_TRIM/tmp_all_concat_trimmed.{wildcards.R}_fastqc.zip {working_directory}/RESULTS/DEMULT_TRIM/all_trimmed.{wildcards.R}_fastqc.zip"


#muse-login.hpc-lr.univ-montp2.fr
