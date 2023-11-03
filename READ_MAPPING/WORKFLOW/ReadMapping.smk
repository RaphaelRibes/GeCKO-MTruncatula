#!/usr/bin/env python

import os,sys,glob
from itertools import compress
from datetime import datetime

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

existing_bed = config["BED"]

if (len(config["BED_MIN_MEAN_COV"]) > 0):
    create_bed = True
else:
    create_bed = False

if (len(existing_bed) == 0 and len(config["BED_MIN_MEAN_COV"]) == 0):
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
outputs_directory = working_directory+"/WORKFLOWS_OUTPUTS/READ_MAPPING"
mapping_dir = outputs_directory+mapping_subfolder
bams_dir = mapping_dir+"/BAMS"
bams_reports_dir = bams_dir+"/REPORTS"
bams_stats_reports_dir = bams_reports_dir+"/STATS"
zones_stats_dir = mapping_dir+"/ZONES_STATS"
subref_dir = mapping_dir+"/EXTRACTED_BAMS/REFERENCE_zones"
subbams_dir = mapping_dir+"/EXTRACTED_BAMS/BAMS_zones"
subbams_reports_dir = subbams_dir+"/REPORTS"
subbams_stats_reports_dir = subbams_reports_dir+"/STATS"

### If the bed file needs to be generated
bed_to_create = subref_dir+"/auto_zones.bed"
clean_bed = subref_dir+"/user_clean.bed"

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


### Generate the workflow_info name
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
workflow_info_file = f"{mapping_dir}/workflow_info_{timestamp}.txt"


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
ruleorder: MarkDuplicates_Bams > Filter_Bams
ruleorder: Mapping_PairedEndFastqs > Filter_Bams
ruleorder: Mapping_SingleEndFastqs > Filter_Bams
ruleorder: Stats_Bams > Filter_Bams
ruleorder: Remapping_PairedEndExtractedFastqs > Filter_Subbams
ruleorder: Remapping_SingleEndExtractedFastqs > Filter_Subbams
ruleorder: Stats_Subbams > Filter_Subbams


rule FinalTargets:
    input:
        mapping_dir+"/summary.sentinel"

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
        buildExpectedFiles([ temp(bams_dir+"/{base}_markedDup.bam"),
        temp(bams_dir+"/{base}_raw.sam"),
        temp(bams_dir+"/{base}_sortcoord.bam"),
        bams_reports_dir+"/DUPLICATES/{base}.bam.metrics" ],

        [ paired_end, paired_end, paired_end, paired_end, paired_end ])
    conda:
        "ENVS/conda_tools.yml"
    threads: config["MAPPING_CPUS_PER_TASK"]
    params:
        mapper = mapper,
        extra_mapper_options = config["EXTRA_MAPPER_OPTIONS"],
        technology = config["SEQUENCING_TECHNOLOGY"],
        picard_markduplicates_options = config["PICARD_MARKDUPLICATES_OPTIONS"],
        picard_markduplicates_java_options = config["PICARD_MARKDUPLICATES_JAVA_OPTIONS"]
    shell:
        "{scripts_dir}/mapping.sh --paired_end --fastq_R1 \"{input.fastq_paired_R1}\" --fastq_R2 \"{input.fastq_paired_R2}\" "
        "--ref {input.ref} --mapper {params.mapper} --mapper_options \"{params.extra_mapper_options}\" --technology \"{params.technology}\" "
        "--output_dir {bams_dir} --sample {wildcards.base} --MD_options \"{params.picard_markduplicates_options}\" "
        "--MD_java_options \"{params.picard_markduplicates_java_options}\" --reports_dir {bams_reports_dir}"


rule Mapping_SingleEndFastqs:
    input:
        fastq_single = lambda wildcards: fastqs_list_dict[wildcards.base],
        ref = ref,
        ref_index = ref_index
    output:
        buildExpectedFiles([ temp(bams_dir+"/{base}_markedDup.bam"),
        temp(bams_dir+"/{base}_raw.sam"),
        temp(bams_dir+"/{base}_sortcoord.bam"),
        bams_reports_dir+"/DUPLICATES/{base}.bam.metrics" ],

        [ not paired_end, not paired_end, not paired_end, not paired_end ])
    conda:
        "ENVS/conda_tools.yml"
    threads: config["MAPPING_CPUS_PER_TASK"]
    params:
        mapper = config["MAPPER"],
        extra_mapper_options = config["EXTRA_MAPPER_OPTIONS"],
        technology = config["SEQUENCING_TECHNOLOGY"],
        picard_markduplicates_options = config["PICARD_MARKDUPLICATES_OPTIONS"],
        picard_markduplicates_java_options = config["PICARD_MARKDUPLICATES_JAVA_OPTIONS"]
    shell:
        "{scripts_dir}/mapping.sh --single_end --fastq \"{input.fastq_single}\" "
        "--ref {input.ref} --mapper {params.mapper} --mapper_options \"{params.extra_mapper_options}\" --technology \"{params.technology}\" "
        "--output_dir {bams_dir} --sample {wildcards.base} --MD_options \"{params.picard_markduplicates_options}\" "
        "--MD_java_options \"{params.picard_markduplicates_java_options}\" --reports_dir {bams_reports_dir}"


rule Stats_Bams:
    input:
        bam = bams_dir+"/{base}_markedDup.bam",
    output:
        bams_stats_reports_dir+"/stats_{base}"
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "samtools stats {input.bam} > {output}"


rule MarkDuplicates_Bams:
    input:
        bams_dir+"/{base}_markedDup.bam"
    output:
        MD_bam = temp(bams_dir+"/{base}_rmDup.bam"),
        metrics = bams_reports_dir+"/DUPLICATES_TMP/{base}.bam.metrics"
    conda:
        "ENVS/conda_tools.yml"
    params:
        picard_markduplicates_options = config["PICARD_MARKDUPLICATES_OPTIONS"],
        picard_markduplicates_java_options = config["PICARD_MARKDUPLICATES_JAVA_OPTIONS"]
    shell:
        "picard {params.picard_markduplicates_java_options} MarkDuplicates -I {input} -O {output.MD_bam} -VALIDATION_STRINGENCY SILENT {params.picard_markduplicates_options} -REMOVE_DUPLICATES TRUE -M {output.metrics}"


rule Filter_Bams:
    input:
        bams_dir+"/{base}_rmDup.bam" if config["REMOVE_DUP"] else bams_dir+"/{base}_markedDup.bam"
    output:
        bams_dir+"/{base}.bam"
    conda:
        "ENVS/conda_tools.yml"
    params:
        samtools_view_filters = config["SAMTOOLS_VIEW_FILTERS1"]
    shell:
        "samtools view -b {params.samtools_view_filters} -o {output} {input}"


rule Summarize_BamsReadsCount:
    input:
        buildExpectedFiles([
        expand("{bams_stats_reports_dir}/stats_{sample}", sample=samples, bams_stats_reports_dir=bams_stats_reports_dir),
        expand("{bams_dir}/{sample}_rmDup.bam", sample=samples, bams_dir=bams_dir),
        expand("{bams_dir}/{sample}.bam", sample=samples, bams_dir=bams_dir) ],

        [ True, config["REMOVE_DUP"], True ])
    output:
        bams_reports_dir+"/nb_reads_per_sample.tsv"
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "{scripts_dir}/summarize_stats.sh --stats_folder {bams_stats_reports_dir} --bams_folder {bams_dir} --rmdup {rm_dup} --output {output};"
        "rm -rf {bams_reports_dir}/DUPLICATES_TMP"


rule MultiQC_Bams:
    input:
        stats_files = expand("{bams_stats_reports_dir}/stats_{sample}", sample=samples, bams_stats_reports_dir=bams_stats_reports_dir),
        nb_reads = bams_reports_dir+"/nb_reads_per_sample.tsv"
    output:
        bams_reports_dir+"/multiQC_ReadMapping_Bams_Report.html",
        temp(bams_reports_dir+"/config_multiQC.yaml")
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "mean_nb_reads=$(awk 'BEGIN{{T=0}}{{T=T+$3}}END{{print T/NR}}' {input.nb_reads} | sed 's/\..*//') ;"
        "{scripts_dir}/make_multiQC_config_file.sh --config_file_base {scripts_dir}/config_multiQC_clean_names.yaml --nb_reads ${{mean_nb_reads}} --output_dir {bams_reports_dir};"
        "multiqc {input.stats_files} -c {bams_reports_dir}/config_multiQC.yaml -o {bams_reports_dir} -n multiQC_ReadMapping_Bams_Report -i ReadMapping_Bams_Report"


rule Index_Bams:
    input:
        bam = bams_dir+"/{base}.bam"
    output:
        csi = bams_dir+"/{base}.bam.csi"
    conda:
        "ENVS/conda_tools.yml"
    params:
        samtools_index_options = config["SAMTOOLS_INDEX_OPTIONS"]
    shell:
        "samtools index -c {params.samtools_index_options} {input.bam}"


rule Create_BamsList:
    input:
        expand("{bams_dir}/{sample}.bam", sample=samples, bams_dir=bams_dir)
    output:
        mapping_dir+"/bams_list.txt"
    shell:
        "ls -d {bams_dir}/*.bam > {output}"


rule Clean_BedFile:
    input:
        existing_bed
    output:
        clean_bed
    conda:
        "ENVS/conda_tools.yml"
    shell:
        """
        sort -k1,1 -k2,2n {input} | awk 'BEGIN{{OFS="\t"; chr=0; start=0; end=0}}{{
          if ($1==chr && $2<=end) {{
            end=$3
          }}
          else {{
            if(end !=0) {{print chr, start, end}} ;
            chr=$1 ; start=$2 ; end=$3
          }}
        }} END{{if(end !=0) {{print chr, start, end}}}}' > {output}
        """


rule Create_BedFile:
    input:
        expand("{bams_dir}/{sample}.bam", sample=samples, bams_dir=bams_dir),
        expand("{bams_dir}/{sample}.bam.csi", sample=samples, bams_dir=bams_dir)
    output:
        bed_to_create
    conda:
        "ENVS/conda_tools.yml"
    params:
        min_cov = float(config["BED_MIN_MEAN_COV"]) * len(expand("{bams_dir}/{sample}.bam", sample=samples, bams_dir=bams_dir)) if len(config["BED_MIN_MEAN_COV"]) > 0 else 0,
        min_dist = config["BED_MIN_DIST"],
        min_length = config["BED_MIN_LENGTH"]
    shell:
        "{scripts_dir}/make_bed_file.sh --input_bams_dir {bams_dir} --output_dir {subref_dir} --min_cov {params.min_cov} --min_dist {params.min_dist} --min_length {params.min_length}"


rule CountReadsZones_Bams:
    input:
        bams = expand("{bams_dir}/{sample}.bam", sample=samples, bams_dir=bams_dir),
        csis = expand("{bams_dir}/{sample}.bam.csi", sample=samples, bams_dir=bams_dir),
        bed = bed_to_create if create_bed else clean_bed
    output:
        zones_stats_dir+"/mean_depth_per_zone_per_sample.tsv"
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "echo ZONE {samples} | sed 's/ /\t/g' > {zones_stats_dir}/mean_depth_per_zone_per_sample.tsv ;"
        "samtools bedcov {input.bed} {input.bams} | awk '{{l=$3-$2+1; printf $1\"_\"$2\"_\"$3\"\t\"; for (i=4;i<NF;i++) printf $i/l\"\t\"; printf $NF/l\"\\n\"}}' >> {zones_stats_dir}/mean_depth_per_zone_per_sample.tsv"


rule Create_SubReference:
    input:
        ref = ref,
        bed = bed_to_create if create_bed else clean_bed
    output:
        subref = subref_dir+"/"+ref_name+"_zones.fasta",
        tmp_bed = temp(subref_dir+"/tmp_zones.bed")
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "awk '{{print $1\":\"$2\"-\"$3}}' {input.bed} > {output.tmp_bed} ;"
        "samtools faidx {input.ref} --region-file {output.tmp_bed} > {output.subref};"
        "sed -i 's/:/_/g ; s/-/_/g' {output.subref} ;"
        "{scripts_dir}/index_ref.sh {output.subref} {mapper}"


rule Extract_PairedEndReads:
    input:
        bams = bams_dir+"/{base}.bam",
        bed = bed_to_create if create_bed else clean_bed
    output:
        tmp_extracted_fastq_paired_R1 = temp(subbams_dir+"/{base}_extract.R1.fastq.gz"),
        tmp_extracted_fastq_paired_R2 = temp(subbams_dir+"/{base}_extract.R2.fastq.gz"),
        tmp_extracted_fastq_unpaired = temp(subbams_dir+"/{base}_extract.U.fastq.gz")
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "{scripts_dir}/extract_PEreads.sh --bam {input.bams} --sample {wildcards.base} --bed_file {input.bed} --output_dir {subbams_dir}"


rule Remapping_PairedEndExtractedFastqs:
    input:
        fastq_paired_R1 = subbams_dir+"/{base}_extract.R1.fastq.gz",
        fastq_paired_R2 = subbams_dir+"/{base}_extract.R2.fastq.gz",
        fastq_unpaired = subbams_dir+"/{base}_extract.U.fastq.gz",
        subref = subref_dir+"/"+ref_name+"_zones.fasta"
    output:
        buildExpectedFiles([ temp(subbams_dir+"/{base}_markedDup.bam"),
        temp(subbams_dir+"/{base}_paired.sam"),
        temp(subbams_dir+"/{base}_unpaired.sam"),
        temp(subbams_dir+"/{base}_raw.sam"),
        temp(subbams_dir+"/{base}_sortcoord.bam"),
        subbams_reports_dir+"/DUPLICATES/{base}.bam.metrics" ],

        [ paired_end, paired_end, paired_end, paired_end, paired_end, paired_end ])
    conda:
        "ENVS/conda_tools.yml"
    threads: config["MAPPING_CPUS_PER_TASK"]
    params:
        mapper = mapper,
        extra_mapper_options = config["EXTRA_MAPPER_OPTIONS"],
        technology = config["SEQUENCING_TECHNOLOGY"],
        picard_markduplicates_options = config["PICARD_MARKDUPLICATES_OPTIONS"],
        picard_markduplicates_java_options = config["PICARD_MARKDUPLICATES_JAVA_OPTIONS"]
    shell:
        "{scripts_dir}/mapping.sh --paired_end --fastq_R1 \"{input.fastq_paired_R1}\" --fastq_R2 \"{input.fastq_paired_R2}\" --fastq_U \"{input.fastq_unpaired}\" "
        "--ref {input.subref} --mapper {params.mapper} --mapper_options \"{params.extra_mapper_options}\" --technology \"{params.technology}\" "
        "--output_dir {subbams_dir} --sample {wildcards.base} --MD_options \"{params.picard_markduplicates_options}\" "
        "--MD_java_options \"{params.picard_markduplicates_java_options}\" --reports_dir {subbams_reports_dir}"


rule Extract_SingleEndReads:
    input:
        bams = bams_dir+"/{base}.bam",
        bed = bed_to_create if create_bed else clean_bed
    output:
        tmp_extract_bams = temp(subbams_dir+"/{base}_extract.bam"),
        tmp_extracted_fastq = temp(subbams_dir+"/{base}_extract.fastq.gz"),
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "samtools view -F 4 -b -L {input.bed} {input.bams} > {output.tmp_extract_bams} ;"
        "picard SamToFastq -I {output.tmp_extract_bams} -F {subbams_dir}/{wildcards.base}_extract.fastq -VALIDATION_STRINGENCY SILENT ;"
        "gzip {subbams_dir}/{wildcards.base}_extract.fastq"


rule Remapping_SingleEndExtractedFastqs:
    input:
        fastq_single = subbams_dir+"/{base}_extract.fastq.gz",
        subref = subref_dir+"/"+ref_name+"_zones.fasta"
    output:
        buildExpectedFiles([ temp(subbams_dir+"/{base}_markedDup.bam"),
        temp(subbams_dir+"/{base}_raw.sam"),
        temp(subbams_dir+"/{base}_sortcoord.bam"),
        subbams_reports_dir+"/DUPLICATES/{base}.bam.metrics" ],

        [ not paired_end, not paired_end, not paired_end, not paired_end ])
    conda:
        "ENVS/conda_tools.yml"
    threads: config["MAPPING_CPUS_PER_TASK"]
    params:
        mapper = mapper,
        extra_mapper_options = config["EXTRA_MAPPER_OPTIONS"],
        technology = config["SEQUENCING_TECHNOLOGY"],
        picard_markduplicates_options = config["PICARD_MARKDUPLICATES_OPTIONS"],
        picard_markduplicates_java_options = config["PICARD_MARKDUPLICATES_JAVA_OPTIONS"]
    shell:
        "{scripts_dir}/mapping.sh --single_end --fastq \"{input.fastq_single}\" "
        "--ref {input.subref} --mapper {params.mapper} --mapper_options \"{params.extra_mapper_options}\" --technology \"{params.technology}\" "
        "--output_dir {subbams_dir} --sample {wildcards.base} --MD_options \"{params.picard_markduplicates_options}\" "
        "--MD_java_options \"{params.picard_markduplicates_java_options}\" --reports_dir {subbams_reports_dir}"


rule Stats_Subbams:
    input:
        bam = subbams_dir+"/{base}_markedDup.bam"
    output:
        subbams_stats_reports_dir+"/stats_{base}"
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "samtools stats {input.bam} > {output}"


rule Filter_Subbams:
    input:
        subbams_dir+"/{base}_markedDup.bam"
    output:
        subbams_dir+"/{base}.bam"
    conda:
        "ENVS/conda_tools.yml"
    params:
        samtools_view_filters = config["SAMTOOLS_VIEW_FILTERS2"]
    shell:
        "samtools view -b {params.samtools_view_filters} -o {output} {input}"


rule Summarize_SubbamsReadsCount:
    input:
        expand("{subbams_stats_reports_dir}/stats_{sample}", sample=samples, subbams_stats_reports_dir=subbams_stats_reports_dir),
        expand("{subbams_dir}/{sample}.bam", sample=samples, subbams_dir=subbams_dir)
    output:
        subbams_reports_dir+"/nb_reads_per_sample.tsv"
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "{scripts_dir}/summarize_stats.sh --stats_folder {subbams_stats_reports_dir} --bams_folder {subbams_dir} --rmdup FALSE --output {output};"


rule MultiQC_Subbams:
    input:
        stats_files = expand("{subbams_stats_reports_dir}/stats_{sample}", sample=samples, subbams_stats_reports_dir=subbams_stats_reports_dir),
        nb_reads = subbams_reports_dir+"/nb_reads_per_sample.tsv"
    output:
        subbams_reports_dir+"/multiQC_ReadMapping_SubBams_Report.html",
        temp(subbams_reports_dir+"/config_multiQC.yaml")
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "mean_nb_reads=$(awk 'BEGIN{{T=0}}{{T=T+$3}}END{{print T/NR}}' {input.nb_reads}| sed 's/\..*//') ;"
        "{scripts_dir}/make_multiQC_config_file.sh --config_file_base {scripts_dir}/config_multiQC_clean_names.yaml --nb_reads ${{mean_nb_reads}} --output_dir {subbams_reports_dir};"
        "multiqc {input.stats_files} -c {subbams_reports_dir}/config_multiQC.yaml -o {subbams_reports_dir} -n multiQC_ReadMapping_SubBams_Report -i ReadMapping_SubBams_Report"


rule Index_Subbams:
    input:
        bam = subbams_dir+"/{base}.bam"
    output:
        csi = subbams_dir+"/{base}.bam.csi"
    conda:
        "ENVS/conda_tools.yml"
    params:
        samtools_index_options = config["SAMTOOLS_INDEX_OPTIONS"]
    shell:
        "samtools index -c {params.samtools_index_options} {input.bam}"


rule Create_SubbamsList:
    input:
        expand("{subbams_dir}/{sample}.bam", sample=samples, subbams_dir=subbams_dir)
    output:
        mapping_dir+"/subbams_list.txt"
    shell:
        "ls -d {subbams_dir}/*.bam > {output}"


rule Create_RefChrSizeFile:
    input:
        ref
    output:
        ref_fai = ref+".fai",
        chr_size = mapping_dir+"/reference_chr_size.txt"
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "samtools faidx {input};"
        "cut -f 1,2 {output.ref_fai} > {output.chr_size}"


rule Write_Summary:
    input:
        buildExpectedFiles(
        [ bams_reports_dir+"/multiQC_ReadMapping_Bams_Report.html",
        bams_reports_dir+"/nb_reads_per_sample.tsv",
        mapping_dir+"/bams_list.txt",
        zones_stats_dir+"/mean_depth_per_zone_per_sample.tsv",
        subref_dir+"/"+ref_name+"_zones.fasta",
        expand("{subbams_dir}/{sample}.bam.csi", sample=samples, subbams_dir=subbams_dir),
        subbams_reports_dir+"/nb_reads_per_sample.tsv",
        subbams_reports_dir+"/multiQC_ReadMapping_SubBams_Report.html",
        mapping_dir+"/subbams_list.txt",
        mapping_dir+"/reference_chr_size.txt" ],
        [ True, True, True, count_reads_zones, create_sub_bams, create_sub_bams, create_sub_bams, create_sub_bams, create_sub_bams, create_sub_bams ])
    output:
        temp(mapping_dir+"/summary.sentinel")
    params:
        latest_info_file = lambda wildcards: find_latest_info_file(mapping_dir),
        new_info_file = workflow_info_file
    shell:
        """
        if [ ! -z "{params.latest_info_file}" ]; then mv {params.latest_info_file} {params.new_info_file} ; fi

        echo -e \"\\t\\t-----------------------------------------------------------------------------\\n\" >> {params.new_info_file}
        echo -e \">>>DATE AND TIME:\" >> {params.new_info_file}
        Date=$(date)
        echo -e \"${{Date}}\\n\" >> {params.new_info_file}
        echo -e \">>>WORKFLOW:\" >> {params.new_info_file}
        echo -e \"https://github.com/GE2POP/GeCKO/tree/main/READ_MAPPING\\n\" >> {params.new_info_file}
        cd {snakefile_dir}
        if git rev-parse --git-dir > /dev/null 2>&1; then echo -e \">>>COMMIT ID:\" >> {params.new_info_file}; git rev-parse HEAD >> {params.new_info_file} ; fi
        cd -
        echo -e \"\\n>>>CONFIG FILE:\" >> {params.new_info_file}
        cat {config[configfile_name]} >> {params.new_info_file}
        echo -e \"\\n\" >> {params.new_info_file}
        echo -e \"\\n>>>CLUSTER CONFIG FILE:\" >> {params.new_info_file}
        cat {config[clusterconfig_name]} >> {params.new_info_file}
        echo -e \"\\n\" >> {params.new_info_file}
        echo -e \"\\n>>>SUMMARY:\" >> {params.new_info_file}
        snakemake --snakefile {snakefile_dir}/ReadMapping.smk --configfile {config[configfile_name]} --summary >> {params.new_info_file}
        echo -e \"\\n\" >> {params.new_info_file}

        touch {output}
        """
