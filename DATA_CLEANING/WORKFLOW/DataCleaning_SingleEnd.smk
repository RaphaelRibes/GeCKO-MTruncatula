#!/usr/bin/env python

import os,sys,glob
from itertools import compress
from datetime import datetime

default_threads = 1 #this will be erased by the user's specifications for each rule in the profile yaml

####################   DEFINE CONFIG VARIABLES BASED ON CONFIG FILE   ####################

### Variables from config file
samples = []
with open(config["ADAPTER_FILE"], "r") as file:
    for row in file:
        samples.append(row.split("\t")[0])

fastq_raw = config["FASTQ"]
user_demult_dir = config["DEMULT_DIR"]

if (len(user_demult_dir) == 0):
    performDemultiplexing = True
else:
    performDemultiplexing = False

if config["UMI"]:
    extractUMI = True
else:
    extractUMI = False

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
rawdata_fastqc_reports_dir = rawdata_reports_dir+"/FASTQC"


if performDemultiplexing:
    demult_dir = outputs_directory+"/DEMULT"
    demult_dir_output = demult_dir
    demult_dir_input = demult_dir
else:
    demult_dir = user_demult_dir
    if extractUMI:
        demult_dir_output = outputs_directory+"/DEMULT_UMI"
        demult_dir_input = demult_dir_output
    else:
        demult_dir_output = outputs_directory+"/DEMULT"
        demult_dir_input = demult_dir

demult_reports_dir = demult_dir_output+"/REPORTS"
demult_fastqc_reports_dir = demult_reports_dir+"/FASTQC"
demult_umitools_reports_dir = demult_reports_dir+"/UMITOOLS_INFOS"
demult_cutadapt_reports_dir = demult_reports_dir+"/CUTADAPT_INFOS"



demult_trim_dir = outputs_directory+"/DEMULT_TRIM"
demult_trim_reports_dir = demult_trim_dir+"/REPORTS"
demult_trim_cutadapt_reports_dir = demult_trim_reports_dir+"/CUTADAPT_INFOS"
demult_trim_fastqc_reports_dir = demult_trim_reports_dir+"/FASTQC"


### Generate the workflow_info name
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
workflow_info_file = f"{outputs_directory}/workflow_info_{timestamp}.txt"


### FUNCTIONS

def buildExpectedFiles(filesNames, isExpected):
    expectedFiles = list(compress(filesNames, isExpected))
    return(expectedFiles)

def find_latest_info_file(directory):
    files = glob.glob(f"{directory}/workflow_info_*.txt")
    if not files:
        return ""
    latest_file = max(files, key=os.path.getctime)
    return latest_file


### PIPELINE ###

ruleorder: Fastqc_ConcatTrimmedFastqs > Fastqc_TrimmedFastqs
ruleorder: Fastqc_ConcatDemultFastqs > Fastqc_DemultFastqs


 # ----------------------------------------------------------------------------------------------- #
 # ----------------------------------------------------------------------------------------------- #

rule FinalTargets:
    input:
        outputs_directory+"/summary.sentinel"

 # ----------------------------------------------------------------------------------------------- #


rule Fastqc_RawFastqs:
    input:
        raw_data_dir+"/{base}.fastq.gz"
    output:
        rawdata_fastqc_reports_dir+"/{base}_fastqc.zip"
    conda:
        "ENVS/conda_tools.yml"
    threads: default_threads
    shell:
        "fastqc -o {rawdata_fastqc_reports_dir} {input}"


rule CountReads_RawFastqs:
    input:
        fastq_raw
    output:
        rawdata_reports_dir+"/Reads_Count_RawData.txt"
    threads: default_threads
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
    threads: default_threads
    shell:
        "{scripts_dir}/demultiplex_with_cutadapt_SE.sh --demultdir {demult_dir} --R {input.fastq_raw} "
        "--barcode_file {input.barcode_file} --cores {params.cores} "
        "--substitutions {params.substitutions};"
        "mv {demult_dir}/demultiplexing_cutadapt.info {demult_cutadapt_reports_dir}"


rule ExtractUMI_DemultFastqs:
    input:
        fastq = demult_dir+"/{base}.fastq.gz"
    output:
        fastq = demult_dir_output+"/{base}.fastq.gz"
    params:
        umitools_extract_options = config["UMITOOLS_EXTRACT_OPTIONS"]
    conda:
        "ENVS/conda_tools.yml"
    threads: default_threads
    shell:
        "umi_tools extract --stdin {input.fastq} --stdout {output.fastq} {params.umitools_extract_options}"


rule CountReads_DemultFastqs:
    input:
        expand("{demult_dir_input}/{sample}.fastq.gz", sample=samples, demult_dir_input=demult_dir_input)
    output:
        demult_reports_dir+"/Reads_Count_Demult.txt"
    threads: default_threads
    shell:
        "{scripts_dir}/fastq_read_count.sh --fastq_dir {demult_dir} --output {output}"



rule Fastqc_DemultFastqs:
    input:
        demult_dir_input+"/{base}.fastq.gz"
    output:
        demult_fastqc_reports_dir+"/{base}_fastqc.zip"
    conda:
        "ENVS/conda_tools.yml"
    threads: default_threads
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
    threads: default_threads
    shell:
        "mean_nb_reads=$(awk 'BEGIN{{T=0}}{{T=T+$2}}END{{print T/NR}}' {input.nb_reads} | sed 's/\..*//');"
        "{scripts_dir}/make_multiQC_config_file.sh --config_file_base {scripts_dir}/config_multiQC_classic.yaml --nb_reads ${{mean_nb_reads}} --output_dir {demult_reports_dir};"
        "multiqc {input.fastqc} -o {demult_reports_dir} -n multiQC_Demult_Report -c {demult_reports_dir}/config_multiQC.yaml -i Demultiplexing_Report -f"



rule Concatenate_DemultFastqs:
    input:
        fastqs_demult = expand("{demult_dir_input}/{sample}.fastq.gz", sample=samples, demult_dir_input=demult_dir_input),
        demult_reads_count = demult_reports_dir+"/Reads_Count_Demult.txt" # necessary to exclude concatenated fastq from read count
    output:
        temp(demult_dir_output+"/All_Samples_Concat_demultiplexed.fastq.gz")
    threads: default_threads
    shell:
        "cat {input.fastqs_demult} > {output}"


rule Fastqc_ConcatDemultFastqs:
    input:
        demult_dir_output+"/All_Samples_Concat_demultiplexed.fastq.gz"
    output:
        demult_fastqc_reports_dir+"/All_Samples_Concat_demultiplexed_fastqc.zip"
    conda:
        "ENVS/conda_tools.yml"
    threads: default_threads
    shell:
        "fastqc -o {demult_fastqc_reports_dir} {input} ;"


rule Trimming_DemultFastqs:
    input:
        fastqs_demult = demult_dir_input+"/{base}.fastq.gz",
        adapter_file = config["ADAPTER_FILE"]
    output:
        demult_trim_dir+"/{base}.fastq.gz",
        demult_trim_cutadapt_reports_dir+"/trimming_cutadapt_{base}.info"
    params:
        quality_cutoff = config["TRIMMING_QUAL"],
        minimum_length = config["TRIMMING_MIN_LENGTH"],
        cores = config["TRIMMING_CORES"]
    conda:
        "ENVS/conda_tools.yml"
    threads: default_threads
    shell:
        "{scripts_dir}/trimming_with_cutadapt_SE.sh --sample {wildcards.base} --trimdir {demult_trim_dir} "
        "--R {input.fastqs_demult} --adapter_file {input.adapter_file} "
        "--cores {params.cores} --quality_cutoff {params.quality_cutoff} --minimum_length {params.minimum_length};"
        "mv {demult_trim_dir}/trimming_cutadapt_{wildcards.base}.info {demult_trim_cutadapt_reports_dir}/trimming_cutadapt_{wildcards.base}.info"


rule CountReads_TrimmedFastqs:
    input:
        expand("{demult_trim_dir}/{sample}.fastq.gz", sample=samples, demult_trim_dir=demult_trim_dir)
    output:
        demult_trim_reports_dir+"/Reads_Count_DemultTrim.txt"
    threads: default_threads
    shell:
        "{scripts_dir}/fastq_read_count.sh --fastq_dir {demult_trim_dir} --output {output}"


rule Fastqc_TrimmedFastqs:
    input:
        demult_trim_dir+"/{base}.fastq.gz"
    output:
        demult_trim_fastqc_reports_dir+"/{base}_fastqc.zip"
    conda:
        "ENVS/conda_tools.yml"
    threads: default_threads
    shell:
        "fastqc -o {demult_trim_fastqc_reports_dir} {input}"


rule MultiQC_TrimmedFastqs:
    input:
        cutadapt = expand("{demult_trim_cutadapt_reports_dir}/trimming_cutadapt_{sample}.info", sample=samples, demult_trim_cutadapt_reports_dir=demult_trim_cutadapt_reports_dir),
        fastqc = expand("{demult_trim_fastqc_reports_dir}/{sample}_fastqc.zip", sample=samples, demult_trim_fastqc_reports_dir=demult_trim_fastqc_reports_dir),
        nb_reads =  demult_trim_reports_dir+"/Reads_Count_DemultTrim.txt"
    output:
        demult_trim_reports_dir+"/multiQC_Trimming_Report.html",
        temp(demult_trim_reports_dir+"/config_multiQC.yaml")
    conda:
        "ENVS/conda_tools.yml"
    threads: default_threads
    shell:
        "mean_nb_reads=$(awk 'BEGIN{{T=0}}{{T=T+$2}}END{{print T/NR}}' {input.nb_reads} | sed 's/\..*//');"
        "{scripts_dir}/make_multiQC_config_file.sh --config_file_base {scripts_dir}/config_multiQC_classic.yaml --nb_reads ${{mean_nb_reads}} --output_dir {demult_trim_reports_dir};"
        "multiqc {input.cutadapt} {input.fastqc} -o {demult_trim_reports_dir} -n multiQC_Trimming_Report -c {demult_trim_reports_dir}/config_multiQC.yaml -i Trimming_Report -f"


rule Concatenate_TrimmedFastqs:
    input:
        fastqs_trim = expand("{demult_trim_dir}/{sample}.fastq.gz", sample=samples, demult_trim_dir=demult_trim_dir),
        demult_trim_reads_count = demult_trim_reports_dir+"/Reads_Count_DemultTrim.txt" # necessary to exclude concatenated fastq from read count
    output:
        temp(demult_trim_dir+"/All_Samples_Concat_trimmed.fastq.gz")
    threads: default_threads
    shell:
        "cat {input.fastqs_trim} > {output}"

rule Fastqc_ConcatTrimmedFastqs:
    input:
        demult_trim_dir+"/All_Samples_Concat_trimmed.fastq.gz"
    output:
        demult_trim_fastqc_reports_dir+"/All_Samples_Concat_trimmed_fastqc.zip"
    conda:
        "ENVS/conda_tools.yml"
    threads: default_threads
    shell:
        "fastqc -o {demult_trim_fastqc_reports_dir} {input} ;"


rule MultiQC_Global:
    input:
        buildExpectedFiles(
        [rawdata_fastqc_reports_dir+"/"+fastq_raw_base+"_fastqc.zip",
        demult_fastqc_reports_dir+"/All_Samples_Concat_demultiplexed_fastqc.zip",
        demult_trim_fastqc_reports_dir+"/All_Samples_Concat_trimmed_fastqc.zip"],

        [performDemultiplexing, True, True]
        )

    output:
        outputs_directory+"/multiQC_DataCleaning_Report.html",
        temp(outputs_directory+"/config_multiQC.yaml")
    conda:
        "ENVS/conda_tools.yml"
    threads: default_threads
    shell:
        "sum_nb_reads=$(awk 'BEGIN{{T=0}}{{T=T+$2}}END{{print T}}' {demult_trim_reports_dir}/Reads_Count_DemultTrim.txt | sed 's/\..*//');"
        "{scripts_dir}/make_multiQC_config_file.sh --config_file_base {scripts_dir}/config_multiQC_keepTrim.yaml --nb_reads ${{sum_nb_reads}} --output_dir {outputs_directory};"
        "multiqc {input} -o {outputs_directory} -n multiQC_DataCleaning_Report -c {outputs_directory}/config_multiQC.yaml -i DataCleaning_Report -f -b 'WARNING: In this report the %Dups is overestimated and not really meaningful because the data results from merging several individual samples.'"


rule Write_Summary:
    input:
        buildExpectedFiles(
        [rawdata_reports_dir+"/Reads_Count_RawData.txt",
        demult_reports_dir+"/Reads_Count_Demult.txt",
        demult_reports_dir+"/multiQC_Demult_Report.html",
        demult_trim_reports_dir+"/Reads_Count_DemultTrim.txt",
        demult_trim_reports_dir+"/multiQC_Trimming_Report.html",
        outputs_directory+"/multiQC_DataCleaning_Report.html"],

        [performDemultiplexing, True, True, True, True, True]
        )
    output:
        temp(outputs_directory+"/summary.sentinel")
    params:
        latest_info_file = lambda wildcards: find_latest_info_file(outputs_directory),
        new_info_file = workflow_info_file
    threads: default_threads
    shell:
        """
        if [ ! -z "{params.latest_info_file}" ]; then mv {params.latest_info_file} {params.new_info_file} ; fi

        echo -e \"\\t\\t-----------------------------------------------------------------------------\\n\" >> {params.new_info_file}
        echo -e \">>>DATE AND TIME:\" >> {params.new_info_file}
        Date=$(date)
        echo -e \"${{Date}}\\n\" >> {params.new_info_file}
        echo -e \">>>WORKFLOW:\" >> {params.new_info_file}
        echo -e \"https://github.com/GE2POP/GeCKO/tree/main/DATA_CLEANING\\n\" >> {params.new_info_file}
        cd {snakefile_dir}
        if git rev-parse --git-dir > /dev/null 2>&1; then echo -e \">>>COMMIT ID:\" >> {params.new_info_file}; git rev-parse HEAD >> {params.new_info_file} ; fi
        cd -
        echo -e \"\\n>>>CONFIG FILE:\" >> {params.new_info_file}
        sed 's/#.*//' {config[configfile_name]} | grep -vP "^\s*$" >> {params.new_info_file}
        if [ {config[clusterprofile_name]} != "NULL" ] ; then
            echo -e \"\\n>>>CLUSTER PROFILE FILE:\" >> {params.new_info_file}
            sed 's/#.*//' {config[clusterprofile_name]} | grep -vP "^\s*$" >> {params.new_info_file}
        fi
        echo -e \"\\n>>>SUMMARY:\" >> {params.new_info_file}
        snakemake --snakefile {snakefile_dir}/DataCleaning_SingleEnd.smk --configfile {config[configfile_name]} --summary >> {params.new_info_file}
        echo -e \"\\n\" >> {params.new_info_file}

        touch {output}
        """
