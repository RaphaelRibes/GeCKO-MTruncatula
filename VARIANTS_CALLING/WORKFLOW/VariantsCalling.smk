#!/usr/bin/env python

import pandas as pd
import os,sys
from itertools import compress

#singularity: "docker://condaforge/mambaforge"


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
reference_base = reference.rsplit('.fa', 1)[0]

### Define paths
path_to_snakefile = workflow.snakefile
snakefile_dir = path_to_snakefile.rsplit('/', 1)[0]
scripts_dir = snakefile_dir+"/SCRIPTS"
working_directory = os.getcwd()

## Define outputs subfolders
outputs_directory = working_directory+"/WORKFLOWS_OUTPUTS/VARIANTS_CALLING"
vc_dir = outputs_directory+vc_subfolder
HaplotypeCaller_dir = vc_dir+"/HAPLOTYPE_CALLER"
GenomicsDBImport_dir = vc_dir+"/GENOMICS_DB_IMPORT"
GenotypeGVCFs_dir = vc_dir+"/GENOTYPE_GVCFS"


def buildExpectedFiles(filesNames, isExpected):
    expectedFiles = list(compress(filesNames, isExpected))
    return(expectedFiles)


 # ----------------------------------------------------------------------------------------------- #

### PIPELINE ###

rule FinalTargets:
    input:
        buildExpectedFiles(
        [GenotypeGVCFs_dir+"/variants_calling.vcf.gz",
        GenotypeGVCFs_dir+"/variants_calling_converted.vcf.gz",
        vc_dir+"/workflow_info.txt"],

        [ not performConvertPositions, performConvertPositions, True ]
        )


rule Index_Reference:
    input:
        reference
    output:
        reference+".fai"
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "samtools faidx {input};"


rule Dictionary_Reference:
    input:
        reference
    output:
        reference_base+".dict"
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "gatk CreateSequenceDictionary REFERENCE={input} OUTPUT={output}"


rule ListIntervalsReference_Dictionary:
    input:
        reference_base+".dict"
    output:
        reference_base+"_intervals_for_GATK.list"
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
    conda:
        "ENVS/conda_tools.yml"
    threads: config["GATK_HAPLOTYPE_CALLER_CPUS_PER_TASK"]
    shell:
        "gatk --java-options \"{params.java_options}\" HaplotypeCaller --reference {input.reference} --input {input.bams} --output {output.vcf} {params.extra_options} -ERC GVCF"


rule List_Haplotype:
    input:
        expand("{HaplotypeCaller_dir}/{sample}.g.vcf.gz", sample=samples, HaplotypeCaller_dir=HaplotypeCaller_dir)
    output:
        HaplotypeCaller_dir+"/vcf.list.txt"
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
    conda:
        "ENVS/conda_tools.yml"
    threads: config["GATK_GENOMICS_DB_IMPORT_CPUS_PER_TASK"]
    shell:
        "mkdir -p {output.tmp_DB};"
        "gatk --java-options \"{params.java_options}\" GenomicsDBImport --sample-name-map {input.vcf_list} --intervals {input.intervals} {params.extra_options} --genomicsdb-workspace-path {output.DB} --tmp-dir {output.tmp_DB}"


rule GenotypeGVCFs:
    input:
        reference = reference,
        DB = GenomicsDBImport_dir
    output:
        vcf_gz = GenotypeGVCFs_dir+"/variants_calling.vcf.gz",
        vcf_gz_tbi = GenotypeGVCFs_dir+"/variants_calling.vcf.gz.tbi",
        tmp_GVCF = temp(directory(GenotypeGVCFs_dir+"/tmp_dir_GVCF"))
    params:
        java_options = config["GATK_GENOTYPE_GVCFS_JAVA_OPTIONS"],
        extra_options = config["GATK_GENOTYPE_GVCFS_EXTRA_OPTIONS"]
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "mkdir -p {output.tmp_GVCF};"
        "gatk --java-options \"{params.java_options}\" GenotypeGVCFs --reference {input.reference} --variant gendb://{input.DB} {params.extra_options} --output {output.vcf_gz} --tmp-dir {output.tmp_GVCF}"


rule ConvertPositions:
    input:
        vcf_gz = GenotypeGVCFs_dir+"/variants_calling.vcf.gz",
        vcf_gz_tbi = GenotypeGVCFs_dir+"/variants_calling.vcf.gz.tbi",
        reference_chr_size = reference_chr_size
    output:
        vcf = temp(GenotypeGVCFs_dir+"/variants_calling.vcf"),
        vcf_converted_gz = GenotypeGVCFs_dir+"/variants_calling_converted.vcf.gz",
        vcf_converted_gz_csi = GenotypeGVCFs_dir+"/variants_calling_converted.vcf.gz.csi"
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "gunzip {input.vcf_gz};"
        "python3 {scripts_dir}/convert_vcf_positions.py --vcf {output.vcf} --chr_size {input.reference_chr_size} --vcf_converted {GenotypeGVCFs_dir}/variants_calling_converted.vcf;"
        "bgzip {GenotypeGVCFs_dir}/variants_calling_converted.vcf;"
        "tabix --csi {output.vcf_converted_gz};"
        "rm {input.vcf_gz_tbi}"


rule Metadata:
    output:
        vc_dir+"/workflow_info.txt"
    shell:
        "echo -e \"Date and time:\" > {vc_dir}/workflow_info.txt;"
        "Date=$(date);"
        "echo -e \"${{Date}}\\n\" >> {vc_dir}/workflow_info.txt;"
        "echo -e \"Workflow:\" >> {vc_dir}/workflow_info.txt;"
        "echo -e \"https://github.com/BioInfo-GE2POP-BLE/CAPTURE_SNAKEMAKE_WORKFLOWS/tree/main/VARIANTS_CALLING\\n\" >> {vc_dir}/workflow_info.txt;"
        "cd {snakefile_dir};"
        "if git rev-parse --git-dir > /dev/null 2>&1; then echo -e \"Commit ID:\" >> {vc_dir}/workflow_info.txt; git rev-parse HEAD >> {vc_dir}/workflow_info.txt ; fi"
