#!/usr/bin/env python

import pandas as pd
import os,sys
from itertools import compress

#singularity: "docker://condaforge/mambaforge"


### Variables from config file
reference = config["REFERENCE"]

### Samples list
bams_list = pd.read_csv(config["BAMS_LIST"], header=None).iloc[:, 0]
bams_list_dict = {}
samples = []
for f in bams_list:
    sample_name = f.rsplit('/', 1)[1].replace('.bam', '')
    bams_list_dict[sample_name] = str(f)
    samples.append(sample_name)



### remove .fa .fas .fasta file extension
reference_base = reference.rsplit('.fa', 1)[0]

### Define paths
path_to_snakefile = workflow.snakefile
snakefile_dir = path_to_snakefile.rsplit('/', 1)[0]
scripts_dir = snakefile_dir+"/SCRIPTS"
working_directory = os.getcwd()

## Define outputs subfolders
outputs_directory = working_directory+"/WORKFLOWS_OUTPUTS/VARIANTS_CALLING"
HaplotypeCaller_dir = outputs_directory+"/HAPLOTYPE_CALLER"
GenomicsDBImport_dir = outputs_directory+"/GENOMICS_DB_IMPORT"
GenotypeGVCFs_dir = outputs_directory+"/GENOTYPE_GVCFS"

 # ----------------------------------------------------------------------------------------------- #

### PIPELINE ###

rule FinalTargets:
    input:
        GenotypeGVCFs_dir+"/variants_calling.vcf.gz"


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
        bams = lambda wildcards: bams_list_dict[wildcards.base]
    output:
        vcf = HaplotypeCaller_dir+"/{base}.g.vcf.gz",
        tbi = HaplotypeCaller_dir+"/{base}.g.vcf.gz.tbi"
    params:
        java_options = config["GATK_HAPLOTYPE_CALLER_JAVA_OPTIONS"],
        extra_options = config["GATK_HAPLOTYPE_CALLER_EXTRA_OPTIONS"]
    conda:
        "ENVS/conda_tools.yml"
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
    shell:
        "mkdir -p {output.tmp_DB};"
        "gatk --java-options \"{params.java_options}\" GenomicsDBImport --sample-name-map {input.vcf_list} --intervals {input.intervals} {params.extra_options} --genomicsdb-workspace-path {output.DB} --tmp-dir {output.tmp_DB}"


rule GenotypeGVCFs:
    input:
        reference = reference,
        DB = GenomicsDBImport_dir
    output:
        GVCF = GenotypeGVCFs_dir+"/variants_calling.vcf.gz",
        tmp_GVCF = temp(directory(GenotypeGVCFs_dir+"/tmp_dir_GVCF"))
    params:
        java_options = config["GATK_GENOTYPE_GVCFS_JAVA_OPTIONS"],
        extra_options = config["GATK_GENOTYPE_GVCFS_EXTRA_OPTIONS"]
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "mkdir -p {output.tmp_GVCF};"
        "gatk --java-options \"{params.java_options}\" GenotypeGVCFs --reference {input.reference} --variant gendb://{input.DB} {params.extra_options} --output {output.GVCF} --tmp-dir {output.tmp_GVCF}"
