#!/usr/bin/env python

import os,sys,glob
from itertools import compress
from datetime import datetime

default_threads = 1 #this will be erased by the user's specifications for each rule in the profile yaml

WF="VARIANT_CALLING"

### Variables from config file
reference = config["REFERENCE"]
vc_subfolder = ""
if (len(config["VARIANT_CALLING_SUBFOLDER"]) > 0):
    vc_subfolder = "/"+config["VARIANT_CALLING_SUBFOLDER"]

### Samples list
bams_list = []
with open(config["BAMS_LIST"], "r") as file:
    for row in file:
        bams_list.append(row.rstrip("\n"))

bams_list_dict = {}
samples = []
for f in bams_list:
    sample_name = f.rsplit('/', 1)[::-1][0].replace('.bam', '')
    bams_list_dict[sample_name] = str(f)
    samples.append(sample_name)

reference_chr_size = config["GENOMIC_REFERENCE_CHR_SIZE"]

if (len(reference_chr_size) > 0):
    performConvertPositions = True
else:
    performConvertPositions = False


### remove .fa .fas .fasta file extension
reference_base = reference.rsplit('.fa', 1)[0].rsplit('.fn', 1)[0]


### Define paths
path_to_snakefile = workflow.snakefile
snakefile_dir = path_to_snakefile.rsplit('/', 1)[0]
scripts_dir = snakefile_dir+"/SCRIPTS"
working_directory = os.getcwd()

## Define outputs subfolders
outputs_directory = f"{working_directory}/WORKFLOWS_OUTPUTS/{WF}"
vc_dir = outputs_directory+vc_subfolder
HaplotypeCaller_dir = vc_dir+"/HAPLOTYPE_CALLER"
GenomicsDBImport_dir = vc_dir+"/GENOMICS_DB_IMPORT"
GenotypeGVCFs_dir = vc_dir+"/GENOTYPE_GVCFS"
GenotypeGVCFs_REPORTS_dir = GenotypeGVCFs_dir+"/REPORTS"

## Expected output
if performConvertPositions:
    final_vcf = GenotypeGVCFs_dir+"/variant_calling_converted.vcf.gz"
    final_vcf_index = final_vcf+".csi"
else:
    final_vcf = GenotypeGVCFs_dir+"/variant_calling.vcf.gz"
    final_vcf_index = final_vcf+".tbi"


### Generate the workflow_info name
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
workflow_info_file = f"{vc_dir}/workflow_info_{timestamp}.txt"

### Container path
GeCKO_container = os.path.abspath(os.path.join(snakefile_dir, "../../launcher_files/container/GeCKO.sif"))


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

rule FinalTargets:
    input:
        vc_dir+"/summary.sentinel",
        final_vcf,
        final_vcf_index

# ----------------------------------------------------------------------------------------------- #


rule Index_Reference:
    input:
        reference
    output:
        reference+".fai"
    singularity:
        GeCKO_container
    threads: default_threads
    shell:
        "samtools faidx {input};"


rule Dictionary_Reference:
    input:
        reference
    output:
        reference_base+".dict"
    singularity:
        GeCKO_container
    threads: default_threads
    shell:
        "gatk CreateSequenceDictionary REFERENCE={input} OUTPUT={output}"


rule ListIntervalsReference_Dictionary:
    input:
        reference_base+".dict"
    output:
        reference_base+"_intervals_for_GATK.list"
    singularity:
        GeCKO_container
    threads: default_threads
    shell:
        "grep '@SQ' {input} | cut -f2,3 | sed 's/SN://' | sed 's/LN://' | awk '{{print $1\":1-\"$2}}' > {output}"


rule HaplotypeCaller:
    input:
        reference = reference,
        fai = reference+".fai",
        bams = lambda wildcards: bams_list_dict[wildcards.base],
        dict = reference_base+".dict"
    output:
        vcf = HaplotypeCaller_dir+"/{base}.g.vcf.gz",
        tbi = HaplotypeCaller_dir+"/{base}.g.vcf.gz.tbi"
    params:
        java_options = config["GATK_HAPLOTYPE_CALLER_JAVA_OPTIONS"],
        extra_options = config["GATK_HAPLOTYPE_CALLER_EXTRA_OPTIONS"]
    singularity:
        GeCKO_container
    threads: default_threads
    shell:
        "gatk --java-options \"{params.java_options}\" HaplotypeCaller --reference {input.reference} --input {input.bams} --output {output.vcf} {params.extra_options} -ERC GVCF"


rule List_Haplotype:
    input:
        expand("{HaplotypeCaller_dir}/{sample}.g.vcf.gz", sample=samples, HaplotypeCaller_dir=HaplotypeCaller_dir)
    output:
        HaplotypeCaller_dir+"/vcf.list.txt"
    singularity:
        GeCKO_container
    threads: default_threads
    shell:
        "for vcf in {input} ; do sample=$(basename ${{vcf}} .g.vcf.gz) ; echo ${{sample}}\"\t\"${{vcf}} ; done > {HaplotypeCaller_dir}/vcf.list.txt"


rule GenomicsDBImport:
    input:
        vcf_list = HaplotypeCaller_dir+"/vcf.list.txt",
        intervals = reference_base+"_intervals_for_GATK.list"
    output:
        DB = directory(GenomicsDBImport_dir),
        tmp_DB = temp(directory(outputs_directory+"/tmp_dir_DB"))
    params:
        java_options = config["GATK_GENOMICS_DB_IMPORT_JAVA_OPTIONS"],
        extra_options = config["GATK_GENOMICS_DB_IMPORT_EXTRA_OPTIONS"]
    singularity:
        GeCKO_container
    threads: default_threads
    shell:
        "mkdir -p {output.tmp_DB};"
        "gatk --java-options \"{params.java_options}\" GenomicsDBImport --sample-name-map {input.vcf_list} --intervals {input.intervals} {params.extra_options} --genomicsdb-workspace-path {output.DB} --tmp-dir {output.tmp_DB}"


rule GenotypeGVCFs:
    input:
        reference = reference,
        DB = GenomicsDBImport_dir
    output:
        vcf_gz = temp(GenotypeGVCFs_dir+"/variant_calling.vcf.gz") if performConvertPositions else GenotypeGVCFs_dir+"/variant_calling.vcf.gz",
        vcf_gz_tbi = temp(GenotypeGVCFs_dir+"/variant_calling.vcf.gz.tbi") if performConvertPositions else GenotypeGVCFs_dir+"/variant_calling.vcf.gz.tbi",
        tmp_GVCF = temp(directory(GenotypeGVCFs_dir+"/tmp_dir_GVCF"))
    params:
        java_options = config["GATK_GENOTYPE_GVCFS_JAVA_OPTIONS"],
        extra_options = config["GATK_GENOTYPE_GVCFS_EXTRA_OPTIONS"]
    singularity:
        GeCKO_container
    threads: default_threads
    shell:
        "mkdir -p {output.tmp_GVCF};"
        "gatk --java-options \"{params.java_options}\" GenotypeGVCFs --reference {input.reference} --variant gendb://{input.DB} {params.extra_options} --output {output.vcf_gz} --tmp-dir {output.tmp_GVCF}"


rule ConvertPositions:
    input:
        vcf_gz = GenotypeGVCFs_dir+"/variant_calling.vcf.gz",
        vcf_gz_tbi = GenotypeGVCFs_dir+"/variant_calling.vcf.gz.tbi",
        reference_chr_size = reference_chr_size
    output:
        vcf = temp(GenotypeGVCFs_dir+"/variant_calling.vcf"),
        vcf_converted_gz = GenotypeGVCFs_dir+"/variant_calling_converted.vcf.gz",
        vcf_converted_gz_csi = GenotypeGVCFs_dir+"/variant_calling_converted.vcf.gz.csi"
    singularity:
        GeCKO_container
    threads: default_threads
    shell:
        "gunzip {input.vcf_gz};"
        "python3 {scripts_dir}/convert_vcf_positions.py --vcf {output.vcf} --chr_size {input.reference_chr_size} --vcf_converted {GenotypeGVCFs_dir}/variant_calling_converted.vcf;"
        "bgzip {GenotypeGVCFs_dir}/variant_calling_converted.vcf;"
        "tabix --csi {output.vcf_converted_gz};"



rule Summarize_GVCFVariables:
    input:
        final_vcf
    output:
        stats_tsv = GenotypeGVCFs_REPORTS_dir+"/variants_stats_VC.tsv",
        DP_tsv = temp(GenotypeGVCFs_REPORTS_dir+"/genotypes_DP_VC.tsv"),
        GT_tsv = temp(GenotypeGVCFs_REPORTS_dir+"/genotypes_GT_VC.tsv"),
        pos_tsv = temp(GenotypeGVCFs_REPORTS_dir+"/variants_pos.tsv"),
        lengths_tsv = temp(GenotypeGVCFs_REPORTS_dir+"/contigs_lengths.tsv")
    singularity:
        GeCKO_container
    threads: default_threads
    shell:
        "{scripts_dir}/extract_variants_stats_from_vcf.sh {input} {output.stats_tsv} {output.DP_tsv} {output.GT_tsv} {output.pos_tsv} {output.lengths_tsv} {GenotypeGVCFs_REPORTS_dir}"


rule Plot_GVCFVariablesHistograms:
    input:
        GenotypeGVCFs_REPORTS_dir+"/variants_stats_VC.tsv"
    output:
        GenotypeGVCFs_REPORTS_dir+"/variants_stats_histograms_VC.pdf"
    singularity:
        GeCKO_container
    threads: default_threads
    shell:
        "python {scripts_dir}/plot_variants_stats_histograms.py --input {input} --output {output}"


rule Plot_GVCFDPBoxplot:
    input:
        DP_tsv = GenotypeGVCFs_REPORTS_dir+"/genotypes_DP_VC.tsv",
        GT_tsv = GenotypeGVCFs_REPORTS_dir+"/genotypes_GT_VC.tsv"
    output:
        GenotypeGVCFs_REPORTS_dir+"/genotypes_DP_boxplot_VC.pdf"
    singularity:
        GeCKO_container
    threads: default_threads
    shell:
        "python {scripts_dir}/plot_DP_boxplot.py --input-DP {input.DP_tsv} --input-GT {input.GT_tsv} --output {output}"


rule Plot_GVCFVariantsAlongGenome:
    input:
        pos_tsv = GenotypeGVCFs_REPORTS_dir+"/variants_pos.tsv",
        lengths_tsv = GenotypeGVCFs_REPORTS_dir+"/contigs_lengths.tsv"
    output:
        GenotypeGVCFs_REPORTS_dir+"/variants_along_genome_VC.pdf"
    singularity:
        GeCKO_container
    threads: default_threads
    shell:
        "python {scripts_dir}/plot_variants_along_genome.py --snp-pos {input.pos_tsv} --contigs-lengths {input.lengths_tsv} --output {output}"


rule Write_Summary:
    input:
        GenotypeGVCFs_REPORTS_dir+"/variants_stats_histograms_VC.pdf",
        GenotypeGVCFs_REPORTS_dir+"/genotypes_DP_boxplot_VC.pdf",
        GenotypeGVCFs_REPORTS_dir+"/variants_along_genome_VC.pdf"
    output:
        temp(vc_dir+"/summary.sentinel")
    params:
        latest_info_file = lambda wildcards: find_latest_info_file(vc_dir),
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
        echo -e \"https://github.com/GE2POP/GeCKO/tree/main/VARIANT_CALLING\\n\" >> {params.new_info_file}
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
        snakemake --snakefile {snakefile_dir}/VariantCalling.smk --configfile {config[configfile_name]} --summary >> {params.new_info_file}
        echo -e \"\\n\" >> {params.new_info_file}

        touch {output}
        """
