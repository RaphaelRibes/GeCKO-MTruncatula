#!/usr/bin/env python

import os,sys
from itertools import compress

#singularity: "docker://condaforge/mambaforge"


####################   DEFINE CONFIG VARIABLES BASED ON CONFIG FILE   ####################

### Variables from config file
samples = []
with open(config["ADAPT_FILE"], "r") as file:
    for row in file:
        samples.append(row.split("\t")[0])

fastq_raw = config["FASTQ"]
#outputs_dirname = config["OUTPUTS_DIRNAME"]
user_demult_dir = config["DEMULT_DIR"]

if (len(user_demult_dir) == 0):
    performDemultiplexing = True
else:
    performDemultiplexing = False

### Raw fastq files path and base names
fastq_raw_base = raw_data_dir = ""
if performDemultiplexing:
    fastq_raw_base = fastq_raw.rsplit('/', 1)[::-1][0].replace('.fastq','').replace('.fq','').replace('.gz','')
    raw_data_dir = fastq_raw.rsplit('/', 1)[0]

### Define paths
path_to_snakefile = workflow.snakefile
snakefile_dir = path_to_snakefile.rsplit('/', 1)[0]
scripts_dir = snakefile_dir+"/SCRIPTS"
working_directory = os.getcwd()

### Define outputs subfolders
outputs_directory = working_directory+"/WORKFLOWS_OUTPUTS/DATA_CLEANING"
rawdata_reports_dir = outputs_directory+"/RAWDATA/REPORTS"

if performDemultiplexing:
    demult_dir = outputs_directory+"/DEMULT"
    demult_dir_output = demult_dir
else:
    demult_dir = user_demult_dir
    demult_dir_output = outputs_directory+"/DEMULT"

demult_reports_dir = demult_dir_output+"/REPORTS"
demult_fastqc_reports_dir = demult_reports_dir+"/FASTQC"

if performDemultiplexing:
    demult_cutadapt_reports_dir = demult_reports_dir+"/CUTADAPT_INFOS"
else:
    demult_cutadapt_reports_dir = ""


demult_trim_dir = outputs_directory+"/DEMULT_TRIM"
demult_trim_reports_dir = demult_trim_dir+"/REPORTS"
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
        demult_reports_dir+"/multiQC_Demult_Report.html",
        demult_trim_reports_dir+"/Reads_Count_DemultTrim.txt",
        demult_trim_reports_dir+"/multiQC_Trimming_Report.html",
        outputs_directory+"/multiQC_DataCleaning_Report.html",
        outputs_directory+"/workflow_info.txt"],

        [performDemultiplexing, True, True, True, True, True, True]
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
        fastq_raw
    output:
        rawdata_reports_dir+"/Reads_Count_RawData.txt"
    shell:
        "{scripts_dir}/fastq_read_count.sh --fastq_dir {raw_data_dir} --output {output}"


rule Demultiplex_RawFastqs:
    input:
        fastq_raw = fastq_raw,
        barcode_file = config["BARCODE_FILE"]
    output:
        expand("{demult_dir}/{sample}.fastq.gz", sample=samples, demult_dir=demult_dir),
        demult_cutadapt_reports_dir+"/demultiplexing_cutadapt.info"
    params:
        substitutions = config["DEMULT_SUBSTITUTIONS"],
        cores = config["DEMULT_CORES"]
    conda:
        "ENVS/conda_tools.yml"
    threads: config["DEMULT_CPUS_PER_TASK"]
    shell:
        "{scripts_dir}/demultiplex_with_cutadapt_SE.sh --demultdir {demult_dir} --R {input.fastq_raw} "
        "--barcode_file {input.barcode_file} --cores {params.cores} "
        "--substitutions {params.substitutions};"
        "mv {demult_dir}/demultiplexing_cutadapt.info {demult_cutadapt_reports_dir}"


rule CountReads_DemultFastqs:
    input:
        expand("{demult_dir}/{sample}.fastq.gz", sample=samples, demult_dir=demult_dir)
    output:
        demult_reports_dir+"/Reads_Count_Demult.txt"
    shell:
        "{scripts_dir}/fastq_read_count.sh --fastq_dir {demult_dir} --output {output}"

rule Fastqc_DemultFastqs:
    input:
        demult_dir+"/{base}.fastq.gz"
    output:
        demult_fastqc_reports_dir+"/{base}_fastqc.zip"
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "fastqc -o {demult_fastqc_reports_dir} {input}"


rule MultiQC_DemultFastqs:
    input:
        fastqc = expand("{demult_fastqc_reports_dir}/{sample}_fastqc.zip", sample=samples, demult_fastqc_reports_dir=demult_fastqc_reports_dir),
        nb_reads = demult_reports_dir+"/Reads_Count_Demult.txt"
    output:
        demult_reports_dir+"/multiQC_Demult_Report.html",
        temp(demult_reports_dir+"/config_multiQC.yaml")
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "mean_nb_reads=$(awk 'BEGIN{{T=0}}{{T=T+$2}}END{{print T/NR}}' {input.nb_reads} | sed 's/\..*//');"
        "{scripts_dir}/make_multiQC_config_file.sh --config_file_base {scripts_dir}/config_multiQC_classic.yaml --nb_reads ${{mean_nb_reads}} --output_dir {demult_reports_dir};"
        "multiqc {input.fastqc} -o {demult_reports_dir} -n multiQC_Demult_Report -c {demult_reports_dir}/config_multiQC.yaml"



rule Concatenate_DemultFastqs:
    input:
        fastqs_demult = expand("{demult_dir}/{sample}.fastq.gz", sample=samples, demult_dir=demult_dir),
        demult_reads_count = demult_reports_dir+"/Reads_Count_Demult.txt" # necessary to exclude concatenated fastq from read count
    output:
        temp(demult_dir_output+"/"+fastq_raw_base+"_demultiplexed.fastq.gz")
    shell:
        "cat {input.fastqs_demult} > {output}"


rule Fastqc_ConcatDemultFastqs:
    input:
        demult_dir_output+"/"+fastq_raw_base+"_demultiplexed.fastq.gz"
    output:
        demult_fastqc_reports_dir+"/All_Samples_Concat_demultiplexed_fastqc.zip"
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "fastqc -o {demult_fastqc_reports_dir} {input} ;"
        "mv {demult_fastqc_reports_dir}/{fastq_raw_base}_demultiplexed_fastqc.html {demult_fastqc_reports_dir}/All_Samples_Concat_demultiplexed_fastqc.html ;"
        "mv {demult_fastqc_reports_dir}/{fastq_raw_base}_demultiplexed_fastqc.zip {demult_fastqc_reports_dir}/All_Samples_Concat_demultiplexed_fastqc.zip"


rule Trimming_DemultFastqs:
    input:
        fastqs_demult = demult_dir+"/{base}.fastq.gz",
        adapt_file = config["ADAPT_FILE"]
    output:
        demult_trim_dir+"/{base}_trimmed.fastq.gz",
        demult_trim_cutadapt_reports_dir+"/trimming_cutadapt_{base}.info"
    params:
        quality_cutoff = config["TRIMMING_QUAL"],
        minimum_length = config["TRIMMING_MIN_LENGTH"],
        cores = config["TRIMMING_CORES"]
    conda:
        "ENVS/conda_tools.yml"
    threads: config["TRIMMING_CPUS_PER_TASK"]
    shell:
        "{scripts_dir}/trimming_with_cutadapt_SE.sh --sample {wildcards.base} --trimdir {demult_trim_dir} "
        "--R {input.fastqs_demult} --adapt_file {input.adapt_file} "
        "--cores {params.cores} --quality_cutoff {params.quality_cutoff} --minimum_length {params.minimum_length};"
        "mv {demult_trim_dir}/trimming_cutadapt_{wildcards.base}.info {demult_trim_cutadapt_reports_dir}/trimming_cutadapt_{wildcards.base}.info"


rule CountReads_TrimmedFastqs:
    input:
        expand("{demult_trim_dir}/{sample}_trimmed.fastq.gz", sample=samples, demult_trim_dir=demult_trim_dir)
    output:
        demult_trim_reports_dir+"/Reads_Count_DemultTrim.txt"
    shell:
        "{scripts_dir}/fastq_read_count.sh --fastq_dir {demult_trim_dir} --output {output}"


rule Fastqc_TrimmedFastqs:
    input:
        demult_trim_dir+"/{base}_trimmed.fastq.gz"
    output:
        demult_trim_fastqc_reports_dir+"/{base}_trimmed_fastqc.zip"
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "fastqc -o {demult_trim_fastqc_reports_dir} {input}"


rule MultiQC_TrimmedFastqs:
    input:
        cutadapt = expand("{demult_trim_cutadapt_reports_dir}/trimming_cutadapt_{sample}.info", sample=samples, demult_trim_cutadapt_reports_dir=demult_trim_cutadapt_reports_dir),
        fastqc = expand("{demult_trim_fastqc_reports_dir}/{sample}_trimmed_fastqc.zip", sample=samples, demult_trim_fastqc_reports_dir=demult_trim_fastqc_reports_dir),
        nb_reads =  demult_trim_reports_dir+"/Reads_Count_DemultTrim.txt"
    output:
        demult_trim_reports_dir+"/multiQC_Trimming_Report.html",
        temp(demult_trim_reports_dir+"/config_multiQC.yaml")
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "mean_nb_reads=$(awk 'BEGIN{{T=0}}{{T=T+$2}}END{{print T/NR}}' {input.nb_reads} | sed 's/\..*//');"
        "{scripts_dir}/make_multiQC_config_file.sh --config_file_base {scripts_dir}/config_multiQC_classic.yaml --nb_reads ${{mean_nb_reads}} --output_dir {demult_trim_reports_dir};"
        "multiqc {input.cutadapt} {input.fastqc} -o {demult_trim_reports_dir} -n multiQC_Trimming_Report -c {demult_trim_reports_dir}/config_multiQC.yaml"


rule Concatenate_TrimmedFastqs:
    input:
        fastqs_trim = expand("{demult_trim_dir}/{sample}_trimmed.fastq.gz", sample=samples, demult_trim_dir=demult_trim_dir),
        demult_trim_reads_count = demult_trim_reports_dir+"/Reads_Count_DemultTrim.txt" # necessary to exclude concatenated fastq from read count
    output:
        temp(demult_trim_dir+"/"+fastq_raw_base+"_trimmed.fastq.gz")
    shell:
        "cat {input.fastqs_trim} > {output}"

rule Fastqc_ConcatTrimmedFastqs:
    input:
        demult_trim_dir+"/"+fastq_raw_base+"_trimmed.fastq.gz"
    output:
        demult_trim_fastqc_reports_dir+"/All_Samples_Concat_trimmed_fastqc.zip"
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "fastqc -o {demult_trim_fastqc_reports_dir} {input} ;"
        "mv {demult_trim_fastqc_reports_dir}/{fastq_raw_base}_trimmed_fastqc.html {demult_trim_fastqc_reports_dir}/All_Samples_Concat_trimmed_fastqc.html ;"
        "mv {demult_trim_fastqc_reports_dir}/{fastq_raw_base}_trimmed_fastqc.zip {demult_trim_fastqc_reports_dir}/All_Samples_Concat_trimmed_fastqc.zip"


rule MultiQC_Global:
    input:
        buildExpectedFiles(
        [rawdata_reports_dir+"/"+fastq_raw_base+"_fastqc.zip",
        demult_fastqc_reports_dir+"/All_Samples_Concat_demultiplexed_fastqc.zip",
        demult_trim_fastqc_reports_dir+"/All_Samples_Concat_trimmed_fastqc.zip"],

        [performDemultiplexing, True, True]
        )

    output:
        outputs_directory+"/multiQC_DataCleaning_Report.html",
        temp(outputs_directory+"/config_multiQC.yaml")
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "sum_nb_reads=$(awk 'BEGIN{{T=0}}{{T=T+$2}}END{{print T}}' {demult_trim_reports_dir}/Reads_Count_DemultTrim.txt | sed 's/\..*//');"
        "{scripts_dir}/make_multiQC_config_file.sh --config_file_base {scripts_dir}/config_multiQC_keepTrim.yaml --nb_reads ${{sum_nb_reads}} --output_dir {outputs_directory};"
        "multiqc {input} -o {outputs_directory} -n multiQC_DataCleaning_Report -c {outputs_directory}/config_multiQC.yaml"



rule Metadata:
    output:
        outputs_directory+"/workflow_info.txt"
    shell:
        "echo -e \"Date and time:\" > {outputs_directory}/workflow_info.txt;"
        "Date=$(date);"
        "echo -e \"${{Date}}\\n\" >> {outputs_directory}/workflow_info.txt;"
        "echo -e \"Workflow:\" >> {outputs_directory}/workflow_info.txt;"
        "echo -e \"https://github.com/BioInfo-GE2POP-BLE/CAPTURE_SNAKEMAKE_WORKFLOWS/tree/main/DATA_CLEANING\\n\" >> {outputs_directory}/workflow_info.txt;"
        "echo -e \"Commit ID:\" >> {outputs_directory}/workflow_info.txt;"
        "cd {snakefile_dir};"
        "git rev-parse HEAD >> {outputs_directory}/workflow_info.txt"


#front.migale.inrae.fr
#muse-login.hpc-lr.univ-montp2.fr
