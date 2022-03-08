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

ref = config["REFERENCE"]
ref_name = config["REFERENCE_NAME"]

bed = config["BED"]
if (len(bed) == 0):
    count_reads_zones = False
else:
    count_reads_zones = True

create_sub_bams = config["CREATE_SUB_BAMS"]
mapper = config["MAPPER"]


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
bams_stats_reports_dir = bams_reports_dir+"/STATS"
subref_dir = mapping_dir+"/SUB_REFERENCE"
subbams_dir = mapping_dir+"/SUB_BAMS"
subbams_reports_dir = subbams_dir+"/REPORTS"
subbams_stats_reports_dir = subbams_reports_dir+"/STATS"

### Expected reference index
if (mapper == "bwa-mem2_mem"):
    ref_index=ref+".bwt.2bit.64"

if (mapper == "bwa_mem"):
    ref_index=ref+".bwt"

if (mapper == "bowtie2"):
    ref_base=ref.replace(".fasta","").replace(".fa","")
    ref_index=ref_base+".1.bt2"

if mapper == "minimap2":
    ref_base=ref.replace(".fasta","").replace(".fa","")
    ref_index=ref_base+".mmi"



### FUNCTIONS
def buildExpectedFiles(filesNames, isExpected):
    expectedFiles = list(compress(filesNames, isExpected))
    return(expectedFiles)



### PIPELINE ###

rule FinalTargets:
    input:
        buildExpectedFiles(
        [ bams_reports_dir+"/bams_multiqc_report.html",
        bams_reports_dir+"/mean_depth_per_zone_per_bam_heatmap.pdf",
        expand("{subbams_dir}/sub_{sample}.bam.bai", sample=samples, subbams_dir=subbams_dir),
        subbams_reports_dir+"/subBams_multiqc_report.html" ],
        [ True, count_reads_zones, create_sub_bams, create_sub_bams ])


# ------------------------------------------------------------------------------------------------------------- #



rule Index_Reference:
    input:
        ref
    output:
        ref_index
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "{scripts_dir}/index_ref.sh {input} {mapper}"


rule Mapping_PairedEndFastqs:
    input:
        fastq_paired_R1 = trim_dir+"/{base}.R1.fastq.gz",
        fastq_paired_R2 = trim_dir+"/{base}.R2.fastq.gz",
        ref = ref,
        ref_index = ref_index
    output:
        buildExpectedFiles([bams_dir+"/{base}.bam"],[paired_end])
    conda:
        "ENVS/conda_tools.yml"
    params:
        mapper = mapper,
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
        ref = ref,
        ref_index = ref_index
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


#rule CountReads_Bams:
#    input:
#        bams_dir+"/{base}.bam"
#    output:
#        bams_reports_dir+"/flagstat_{base}"
#    conda:
#        "ENVS/conda_tools.yml"
#    shell:
#        "samtools flagstat {input} > {output}"


#rule Summarize_PairedEndReadsCount:
#    input:
#        expand("{bams_reports_dir}/flagstat_{sample}", sample=samples, bams_reports_dir=bams_reports_dir)
#    output:
#        buildExpectedFiles([bams_reports_dir+"/summary_flagstat.tsv"],[paired_end])
#    shell:
#        "{scripts_dir}/summarize_flagstat.sh --pairedEnd --flagstat_folder {bams_reports_dir}"


#rule Summarize_SingleEndReadsCount:
#    input:
#        expand("{bams_reports_dir}/flagstat_{sample}", sample=samples, bams_reports_dir=bams_reports_dir)
#    output:
#        buildExpectedFiles([bams_reports_dir+"/summary_flagstat.tsv"],[not paired_end])
#    shell:
#        "{scripts_dir}/summarize_flagstat.sh --singleEnd --flagstat_folder {bams_reports_dir}"


#rule Histogram_ReadsCount:
#    input:
#        bams_reports_dir+"/summary_flagstat.tsv"
#    output:
#        bams_reports_dir+"/reads_histograms.pdf"
#    conda:
#        "ENVS/conda_tools.yml"
#    shell:
#        "python3 {scripts_dir}/plot_reads_histogram.py {bams_reports_dir}"

rule Stats_Bams:
    input:
        bams_dir+"/{base}.bam"
    output:
        bams_stats_reports_dir+"/stats_{base}"
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "samtools stats {input} > {output}"


rule MultiQC_Bams:
    input:
        expand("{bams_stats_reports_dir}/stats_{sample}", sample=samples, bams_stats_reports_dir=bams_stats_reports_dir)
    output:
        bams_reports_dir+"/bams_multiqc_report.html"
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "multiqc {input} -n {output}"


rule CountReadsZones_Bams:
    input:
        expand("{bams_dir}/{sample}.bam", sample=samples, bams_dir=bams_dir)
    output:
        buildExpectedFiles([bams_reports_dir+"/mean_depth_per_zone_per_bam.tsv"], [count_reads_zones])
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "echo ZONE {samples} | sed 's/ /\t/g' > {bams_reports_dir}/mean_depth_per_zone_per_bam.tsv ;"
        "samtools bedcov {bed} {input} | awk '{{l=$3-$2+1; printf $1\"_\"$2\"_\"$3\"\t\"; for (i=4;i<=NF;i++)printf $i/l\"\t\"; print \"\"}}' >> {bams_reports_dir}/mean_depth_per_zone_per_bam.tsv"


rule Heatmap_ZonesReadsCount:
    input:
        bams_reports_dir+"/mean_depth_per_zone_per_bam.tsv"
    output:
        bams_reports_dir+"/mean_depth_per_zone_per_bam_heatmap.pdf"
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "python3 {scripts_dir}/plot_reads_heatmap.py {bams_reports_dir}"


rule Create_SubReference:
    input:
        ref = ref,
        bed = bed
    output:
        subref = subref_dir+"/sub_"+ref_name+".fasta",
        tmp_bed = temp(subref_dir+"/regions.bed"),
        dict = subref_dir+"/sub_"+ref_name+".dict"
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "awk '{{print $1\":\"$2\"-\"$3}}' {input.bed} > {output.tmp_bed} ;"
        "samtools faidx {input.ref} --region-file {subref_dir}/regions.bed > {output.subref};"
        "picard CreateSequenceDictionary R={output.subref} O={output.dict}"


rule Extract_Reads:
    input:
        dict = subref_dir+"/sub_"+ref_name+".dict",
        bams = bams_dir+"/{base}.bam"
    output:
        subbams_dir+"/sub_{base}.bam"
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "picard ReorderSam INPUT={input.bams} OUTPUT={output} SEQUENCE_DICTIONARY={input.dict} S=true VERBOSITY=WARNING"


rule Index_Subbams:
    input:
        subbams_dir+"/{base}.bam"
    output:
        buildExpectedFiles([subbams_dir+"/{base}.bam.bai"],[create_sub_bams])
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "samtools index {input}"


rule Stats_Subbams:
    input:
        subbams_dir+"/{base}.bam"
    output:
        subbams_stats_reports_dir+"/stats_{base}"
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "samtools stats {input} > {output}"


rule MultiQC_Subbams:
    input:
        expand("{subbams_stats_reports_dir}/stats_sub_{sample}", sample=samples, subbams_stats_reports_dir=subbams_stats_reports_dir)
    output:
        buildExpectedFiles([subbams_reports_dir+"/subBams_multiqc_report.html"],[create_sub_bams])
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "multiqc {input} -n {output}"


## + Créer liste des bams pour appel des SNPs

#front.migale.inrae.fr
#muse-login.hpc-lr.univ-montp2.fr
