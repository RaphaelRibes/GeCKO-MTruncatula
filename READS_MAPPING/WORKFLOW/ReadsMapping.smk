#!/usr/bin/env python

import pandas as pd
import os,sys
from itertools import compress
import gzip

#singularity: "docker://condaforge/mambaforge"


####################   DEFINE CONFIG VARIABLES BASED ON CONFIG FILE   ####################

### Define paths
path_to_snakefile = workflow.snakefile
snakefile_dir = path_to_snakefile.rsplit('/', 1)[0]
scripts_dir = snakefile_dir+"/SCRIPTS"
working_directory = os.getcwd()

### Variables from config file
paired_end = config["PAIRED_END"]
trim_dirs = list(config["TRIM_DIRS"].split(" "))

if (len(config["TRIM_DIRS"]) == 0):
    trim_dirs = [working_directory+"/WORKFLOWS_OUTPUTS/DATA_CLEANING/DEMULT_TRIM"]


if config["REMOVE_DUP"]:
    rm_dup = "TRUE"
else:
    rm_dup = "FALSE"

ref = os.path.abspath(config["REFERENCE"])
mapping_subfolder = ""
if (len(config["MAPPING_SUBFOLDER"]) > 0):
    mapping_subfolder = "/"+config["MAPPING_SUBFOLDER"]

ref_name = ref.rsplit('/', 1)[::-1][0].replace('.fasta','').replace('.fas','').replace('.fa','')

bed = config["BED"]
if (len(bed) == 0):
    count_reads_zones = False
else:
    count_reads_zones = True

if config["CREATE_SUB_BAMS"]:
    create_sub_bams = True
else:
    create_sub_bams = False

mapper = config["MAPPER"]




### Samples list
if paired_end:
    end="--pairedEnd"

    samples = []
    fastqs_list_dict = {}

    for trim_dir in trim_dirs:
        fastqs_R1_list = [fastq for fastq in os.listdir(trim_dir) if '.R1.fastq.gz' in fastq]

        fastqs_R2_list = []
        for fastq_R1 in fastqs_R1_list:
            fastq_R2 = fastq_R1.replace(".R1.", ".R2.")
            fastqs_R2_list.append(fastq_R2)
            sample = fastq_R1.replace(".R1.fastq.gz", "")
            samples.append(sample)
            fastqs_list_dict[sample] = [trim_dir+"/"+str(fastq_R1), trim_dir+"/"+str(fastq_R2)]

else:
    end="--singleEnd"

    samples = []
    fastqs_list_dict = {}

    for trim_dir in trim_dirs:
        fastqs_list = [fastq for fastq in os.listdir(trim_dir) if '.fastq.gz' in fastq]

        for fastq in fastqs_list:
            sample = fastq.replace(".fastq.gz", "")
            samples.append(sample)
            fastqs_list_dict[sample] = trim_dir+"/"+str(fastq)


### Define outputs subfolders
outputs_directory = working_directory+"/WORKFLOWS_OUTPUTS/READS_MAPPING"
mapping_dir = outputs_directory+mapping_subfolder
bams_dir = mapping_dir+"/BAMS"
bams_reports_dir = bams_dir+"/REPORTS"
bams_stats_reports_dir = bams_reports_dir+"/STATS"
zones_stats_dir = mapping_dir+"/ZONES_STATS"
subref_dir = mapping_dir+"/EXTRACTED_BAMS/REFERENCE_zones"
subbams_dir = mapping_dir+"/EXTRACTED_BAMS/BAMS_zones"
subbams_reports_dir = subbams_dir+"/REPORTS"
subbams_stats_reports_dir = subbams_reports_dir+"/STATS"

### Expected reference index
if (mapper == "bwa-mem2_mem"):
    ref_index=ref+".bwt.2bit.64"

if (mapper == "bwa_mem"):
    ref_index=ref+".bwt"

if (mapper == "bowtie2"):
    ref_base=ref.replace(".fasta","").replace(".fas","").replace(".fa","")
    ref_index=ref_base+".1.bt2"

if mapper == "minimap2":
    ref_base=ref.replace(".fasta","").replace(".fas","").replace(".fa","")
    ref_index=ref_base+".mmi"



### FUNCTIONS
def buildExpectedFiles(filesNames, isExpected):
    expectedFiles = list(compress(filesNames, isExpected))
    return(expectedFiles)



### PIPELINE ###

rule FinalTargets:
    input:
        buildExpectedFiles(
        [ mapping_dir+"/workflow_info.txt",
        bams_reports_dir+"/multiQC_ReadsMapping_Bams_Report.html",
        bams_reports_dir+"/nb_reads_per_sample.tsv",
        mapping_dir+"/bams_list.txt",
        zones_stats_dir+"/mean_depth_per_zone_per_sample_heatmap.pdf",
        subref_dir+"/"+ref_name+"_zones.fasta",
        expand("{subbams_dir}/{sample}_zones.bam.bai", sample=samples, subbams_dir=subbams_dir),
        subbams_reports_dir+"/nb_reads_per_sample.tsv",
        subbams_reports_dir+"/multiQC_ReadsMapping_SubBams_Report.html",
        mapping_dir+"/subbams_list.txt" ],
        [ True, True, True, True, count_reads_zones, create_sub_bams, create_sub_bams, create_sub_bams, create_sub_bams, create_sub_bams ])

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
        fastq_paired_R1 = lambda wildcards: fastqs_list_dict[wildcards.base][0],
        fastq_paired_R2 = lambda wildcards: fastqs_list_dict[wildcards.base][1],
        ref = ref,
        ref_index = ref_index
    output:
        buildExpectedFiles([bams_dir+"/{base}.bam"],[paired_end])
    conda:
        "ENVS/conda_tools.yml"
    threads: config["MAPPING_CPUS_PER_TASK"]
    params:
        mapper = mapper,
        extra_mapper_options = config["EXTRA_MAPPER_OPTIONS"],
        technology = config["SEQUENCING_TECHNOLOGY"],
        picard_markduplicates_options = config["PICARD_MARKDUPLICATES_OPTIONS"],
        picard_markduplicates_java_options = config["PICARD_MARKDUPLICATES_JAVA_OPTIONS"],
        samtools_index_options = config["SAMTOOLS_INDEX_OPTIONS"]
    shell:
        "{scripts_dir}/mapping.sh --paired_end --fastq_R1 \"{input.fastq_paired_R1}\" --fastq_R2 \"{input.fastq_paired_R2}\" "
        "--ref {input.ref} --mapper {params.mapper} --mapper_options \"{params.extra_mapper_options}\" --technology \"{params.technology}\" "
        "--output_dir {bams_dir} --reports_dir {bams_reports_dir} --sample {wildcards.base} --rm_dup {rm_dup} --picard_markduplicates_options \"{params.picard_markduplicates_options}\" "
        "--picard_markduplicates_java_options \"{params.picard_markduplicates_java_options}\" --samtools_index_options \"{params.samtools_index_options}\""


rule Mapping_SingleEndFastqs:
    input:
        fastq_single = lambda wildcards: fastqs_list_dict[wildcards.base],
        ref = ref,
        ref_index = ref_index
    output:
        buildExpectedFiles([bams_dir+"/{base}.bam"],[not paired_end])
    conda:
        "ENVS/conda_tools.yml"
    threads: config["MAPPING_CPUS_PER_TASK"]
    params:
        mapper = config["MAPPER"],
        extra_mapper_options = config["EXTRA_MAPPER_OPTIONS"],
        technology = config["SEQUENCING_TECHNOLOGY"],
        picard_markduplicates_options = config["PICARD_MARKDUPLICATES_OPTIONS"],
        picard_markduplicates_java_options = config["PICARD_MARKDUPLICATES_JAVA_OPTIONS"],
        samtools_index_options = config["SAMTOOLS_INDEX_OPTIONS"]
    shell:
        "{scripts_dir}/mapping.sh --single_end --fastq \"{input.fastq_single}\" "
        "--ref {input.ref} --mapper {params.mapper} --mapper_options \"{params.extra_mapper_options}\" --technology \"{params.technology}\" "
        "--output_dir {bams_dir} --reports_dir {bams_reports_dir} --sample {wildcards.base} --rm_dup {rm_dup} --picard_markduplicates_options \"{params.picard_markduplicates_options}\" "
        "--picard_markduplicates_java_options \"{params.picard_markduplicates_java_options}\" --samtools_index_options \"{params.samtools_index_options}\""


rule Stats_Bams:
    input:
        bams_dir+"/{base}.bam"
    output:
        bams_stats_reports_dir+"/stats_{base}"
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "samtools stats {input} > {output}"


rule Summarize_BamsReadsCount:
    input:
        expand("{bams_stats_reports_dir}/stats_{sample}", sample=samples, bams_stats_reports_dir=bams_stats_reports_dir)
    output:
        bams_reports_dir+"/nb_reads_per_sample.tsv"
    shell:
        "{scripts_dir}/summarize_stats.sh {end} --stats_folder {bams_stats_reports_dir} --output {output};"
	    "rm -r {bams_dir}/TMP"



rule MultiQC_Bams:
    input:
        stats_files = expand("{bams_stats_reports_dir}/stats_{sample}", sample=samples, bams_stats_reports_dir=bams_stats_reports_dir),
        nb_reads = bams_reports_dir+"/nb_reads_per_sample.tsv"
    output:
        bams_reports_dir+"/multiQC_ReadsMapping_Bams_Report.html",
        temp(bams_reports_dir+"/config_multiQC.yaml")
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "mean_nb_reads=$(awk 'BEGIN{{T=0}}{{T=T+$2}}END{{print T/NR}}' {input.nb_reads} | sed 's/\..*//') ;"
        "{scripts_dir}/make_multiQC_config_file.sh --config_file_base {scripts_dir}/config_multiQC_clean_names.yaml --nb_reads ${{mean_nb_reads}} --output_dir {bams_reports_dir};"
        "multiqc {input.stats_files} -c {bams_reports_dir}/config_multiQC.yaml -o {bams_reports_dir} -n multiQC_ReadsMapping_Bams_Report"

        #"if [[ $mean_nb_reads -lt 1000000 ]] ; then kReads=\"_kReads\" ; else kReads=\"\" ; fi ;"
        #"multiqc {input.stats_files} -c {scripts_dir}/config_multiQC_clean_names${{kReads}}.yaml -o {bams_reports_dir} -n multiQC_ReadsMapping_Bams_Report"


rule Create_BamsList:
    input:
        expand("{bams_dir}/{sample}.bam", sample=samples, bams_dir=bams_dir)
    output:
        mapping_dir+"/bams_list.txt"
    shell:
        "ls -d {bams_dir}/*.bam > {output}"


rule CountReadsZones_Bams:
    input:
        expand("{bams_dir}/{sample}.bam", sample=samples, bams_dir=bams_dir)
    output:
        zones_stats_dir+"/mean_depth_per_zone_per_sample.tsv"
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "echo ZONE {samples} | sed 's/ /\t/g' > {zones_stats_dir}/mean_depth_per_zone_per_sample.tsv ;"
        "samtools bedcov {bed} {input} | awk '{{l=$3-$2+1; printf $1\"_\"$2\"_\"$3\"\t\"; for (i=4;i<NF;i++) printf $i/l\"\t\"; printf $NF/l\"\\n\"}}' >> {zones_stats_dir}/mean_depth_per_zone_per_sample.tsv"


rule Heatmap_ZonesReadsCount:
    input:
        zones_stats_dir+"/mean_depth_per_zone_per_sample.tsv"
    output:
        zones_stats_dir+"/mean_depth_per_zone_per_sample_heatmap.pdf"
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "python3 {scripts_dir}/plot_reads_heatmap.py {zones_stats_dir}"


rule Create_SubReference:
    input:
        ref = ref,
        bed = bed
    output:
        subref = subref_dir+"/"+ref_name+"_zones.fasta",
        tmp_bed = temp(subref_dir+"/regions.bed")
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "awk '{{print $1\":\"$2\"-\"$3}}' {input.bed} > {output.tmp_bed} ;"
        "samtools faidx {input.ref} --region-file {subref_dir}/regions.bed > {output.subref};"
        "sed -i 's/:/_/g ; s/-/_/g' {output.subref}"


rule Create_ChainFile:
    input:
        ref = ref,
        bed = bed
    output:
        chain = temp(subref_dir+"/zones.chain")
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "ln -s {input.ref} {subref_dir}/ref.fasta && samtools faidx {subref_dir}/ref.fasta;"
        "awk '{{if(NR==FNR){{lg[$1]=$2}} else{{ print \"chain 4900 \"$1\" \"lg[$1]\" + \"$2\" \"$3\" \"$1\"_\"$2\"_\"$3\" \"$3-$2+1\" + 1 \"$3-$2+1\"\\n\"$3-$2+1 }}}}' {subref_dir}/ref.fasta.fai {input.bed} > {subref_dir}/zones.chain ;"
        "rm {subref_dir}/ref.fasta*"


rule Extract_Reads:
    input:
        bams = bams_dir+"/{base}.bam",
        chain = subref_dir+"/zones.chain"
    output:
        bam = subbams_dir+"/{base}_zones.bam",
        bai = subbams_dir+"/{base}_zones.bam.bai",
        metrics = subbams_reports_dir+"/DUPLICATES/{base}.bam.metrics"
    conda:
        "ENVS/conda_tools.yml"
    params:
        picard_markduplicates_options = config["PICARD_MARKDUPLICATES_OPTIONS"],
        picard_markduplicates_java_options = config["PICARD_MARKDUPLICATES_JAVA_OPTIONS"]
    shell:
        "samtools view -b -F 0x800 -F 0x100 {input.bams} > {subbams_dir}/{wildcards.base}_primary.bam ;"
        "CrossMap.py bam -a {input.chain} {subbams_dir}/{wildcards.base}_primary.bam {subbams_dir}/{wildcards.base}_zones;"
        "rm {subbams_dir}/{wildcards.base}_zones.sorted.bam.bai;"
        "picard SortSam --TMP_DIR {subbams_dir}/TMP -I {subbams_dir}/{wildcards.base}_zones.sorted.bam -O {subbams_dir}/{wildcards.base}_zones.sortname.bam -SO queryname -VALIDATION_STRINGENCY SILENT;"
        "samtools fixmate {subbams_dir}/{wildcards.base}_zones.sortname.bam {subbams_dir}/{wildcards.base}_zones.fix.bam ;"
        "picard SortSam --TMP_DIR {subbams_dir}/TMP -I {subbams_dir}/{wildcards.base}_zones.fix.bam -O {subbams_dir}/{wildcards.base}_zones.sortcoord.bam -SO coordinate -VALIDATION_STRINGENCY SILENT;"
        "picard {params.picard_markduplicates_java_options} MarkDuplicates -I {subbams_dir}/{wildcards.base}_zones.sortcoord.bam -O {output.bam} -VALIDATION_STRINGENCY SILENT {params.picard_markduplicates_options} -REMOVE_DUPLICATES FALSE -M {output.metrics};"
        "samtools index {output.bam};"
        "rm {subbams_dir}/{wildcards.base}_primary.bam {subbams_dir}/{wildcards.base}_zones.sorted.bam {subbams_dir}/{wildcards.base}_zones.sortname.bam {subbams_dir}/{wildcards.base}_zones.sortcoord.bam {subbams_dir}/{wildcards.base}_zones.fix.bam"


rule Stats_Subbams:
    input:
        subbams_dir+"/{base}.bam"
    output:
        subbams_stats_reports_dir+"/stats_{base}"
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "samtools stats {input} > {output}"


rule Summarize_SubbamsReadsCount:
    input:
        expand("{subbams_stats_reports_dir}/stats_{sample}_zones", sample=samples, subbams_stats_reports_dir=subbams_stats_reports_dir)
    output:
        subbams_reports_dir+"/nb_reads_per_sample.tsv"
    shell:
        "{scripts_dir}/summarize_stats.sh {end} --stats_folder {subbams_stats_reports_dir} --output {output};"
	    "rm -r {subbams_dir}/TMP"


rule MultiQC_Subbams:
    input:
        stats_files = expand("{subbams_stats_reports_dir}/stats_{sample}_zones", sample=samples, subbams_stats_reports_dir=subbams_stats_reports_dir),
        nb_reads = subbams_reports_dir+"/nb_reads_per_sample.tsv"
    output:
        subbams_reports_dir+"/multiQC_ReadsMapping_SubBams_Report.html",
        temp(subbams_reports_dir+"/config_multiQC.yaml")
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "mean_nb_reads=$(awk 'BEGIN{{T=0}}{{T=T+$2}}END{{print T/NR}}' {input.nb_reads}| sed 's/\..*//') ;"
        "{scripts_dir}/make_multiQC_config_file.sh --config_file_base {scripts_dir}/config_multiQC_clean_names.yaml --nb_reads ${{mean_nb_reads}} --output_dir {subbams_reports_dir};"
        "multiqc {input.stats_files} -c {subbams_reports_dir}/config_multiQC.yaml -o {subbams_reports_dir} -n multiQC_ReadsMapping_SubBams_Report"

        #"if [[ $mean_nb_reads -lt 1000000 ]] ; then kReads=\"_kReads\" ; else kReads=\"\" ; fi ;"
        #"multiqc {input.stats_files} -c {scripts_dir}/config_multiQC_clean_names${{kReads}}.yaml -o {subbams_reports_dir} -n multiQC_ReadsMapping_SubBams_Report"


rule Create_SubbamsList:
    input:
        expand("{subbams_dir}/{sample}_zones.bam", sample=samples, subbams_dir=subbams_dir)
    output:
        mapping_dir+"/subbams_list.txt"
    shell:
        "ls -d {subbams_dir}/*.bam > {output}"


rule Metadata:
    output:
        mapping_dir+"/workflow_info.txt"
    shell:
        "echo -e \"Date and time:\" > {mapping_dir}/workflow_info.txt;"
        "Date=$(date);"
        "echo -e \"${{Date}}\\n\" >> {mapping_dir}/workflow_info.txt;"
        "echo -e \"Workflow:\" >> {mapping_dir}/workflow_info.txt;"
        "echo -e \"https://github.com/BioInfo-GE2POP-BLE/CAPTURE_SNAKEMAKE_WORKFLOWS/tree/main/READS_MAPPING\\n\" >> {mapping_dir}/workflow_info.txt;"
        "echo -e \"Commit ID:\" >> {mapping_dir}/workflow_info.txt;"
        "cd {snakefile_dir};"
        "git rev-parse HEAD >> {mapping_dir}/workflow_info.txt"



#front.migale.inrae.fr
#muse-login.hpc-lr.univ-montp2.fr
