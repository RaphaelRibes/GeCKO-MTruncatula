#!/usr/bin/env python

import os,sys
from itertools import compress

#singularity: "docker://condaforge/mambaforge"


####################   DEFINE CONFIG VARIABLES BASED ON CONFIG FILE   ####################

### Variables from config file
vcf_raw = config["VCF_FILE"]

### Define paths
path_to_snakefile = workflow.snakefile
snakefile_dir = path_to_snakefile.rsplit('/', 1)[0]
scripts_dir = snakefile_dir+"/SCRIPTS"
working_directory = os.getcwd()

### Define outputs subfolders
outputs_directory = working_directory+"/WORKFLOWS_OUTPUTS/VCF_FILTERING"
VCF_reports_dir = outputs_directory+"/REPORTS"

### PIPELINE ###

rule FinalTargets:
    input:
        VCF_reports_dir+"/multiQC_VcfFiltering_report.html",
        outputs_directory+"/workflow_info.txt"


 # ----------------------------------------------------------------------------------------------- #


rule Filter_Loci:
    input:
        vcf_raw
    output:
        tmp_Locus_Filtering = temp(outputs_directory+"/tmp_Locus_Filtered.recode.vcf"),
        Locus_Filtered = outputs_directory+"/01_Locus_Filtered.vcf"
    conda:
        "ENVS/conda_tools.yml"
    params:
        config["VCFTOOLS_LOCUS_FILTERING_OPTIONS"]
    shell:
        "vcftools --gzvcf {input} {params} --recode --recode-INFO-all --out {outputs_directory}/tmp_Locus_Filtered;"
        "grep '#' {outputs_directory}/tmp_Locus_Filtered.recode.vcf > {output.Locus_Filtered};"
        "grep -v '#' {outputs_directory}/tmp_Locus_Filtered.recode.vcf | sort -k1,1 -k2,2n >> {output.Locus_Filtered}"


rule Filter_Samples:
    input:
        outputs_directory+"/01_Locus_Filtered.vcf"
    output:
        out_imiss = temp(outputs_directory+"/SampleFilter.imiss"),
        samples_to_remove = outputs_directory+"/samples_to_remove.list",
        SampleLocus_Filtered = outputs_directory+"/02_SampleLocus_Filtered.vcf"
    conda:
        "ENVS/conda_tools.yml"
    params:
        config["MAX_RATIO_NA_PER_SAMPLE"]
    shell:
        "vcftools --vcf {input} --missing-indv --out {outputs_directory}/SampleFilter;"
        "awk -v max_pc_NA={params} '{{if($5>max_pc_NA){{print $0}}}}' {output.out_imiss} | cut -f1 > {output.samples_to_remove};"
        "vcftools --vcf {input} --remove {output.samples_to_remove} --recode --recode-INFO-all --out {outputs_directory}/02_SampleLocus_Filtered;"
        "mv {outputs_directory}/02_SampleLocus_Filtered.recode.vcf {output.SampleLocus_Filtered}"


rule Calculate_PopGenStats:
    input:
        SampleLocus_Filtered = outputs_directory+"/02_SampleLocus_Filtered.vcf"
    output:
        outputs_directory+"/02_SampleLocus_Filtered_withPopGenStats.vcf"
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "python {scripts_dir}/egglib_PopGenStats.py --input {input} --output {output}"


rule Filter_PopGenStats:
    input:
        outputs_directory+"/02_SampleLocus_Filtered_withPopGenStats.vcf"
    output:
        temp_SampleLocus_gz = temp(outputs_directory+"/02_SampleLocus_Filtered_withPopGenStats.vcf.gz"),
        temp_SampleLocus_csi = temp(outputs_directory+"/02_SampleLocus_Filtered_withPopGenStats.vcf.gz.csi"),
        PopGenStatsSampleLocus_Filtered = outputs_directory+"/03_PopGenStatsSampleLocus_Filtered.vcf"
    conda:
        "ENVS/conda_tools.yml"
    params:
        config["BCFTOOLS_FILTERING_OPTIONS"]
    shell:
        "bgzip -c {input} > {output.temp_SampleLocus_gz};"
        "tabix --csi {output.temp_SampleLocus_gz};"
        "bcftools filter -sFilterSmk -i '{params}' {output.temp_SampleLocus_gz} | bcftools view -f 'PASS' > {output.PopGenStatsSampleLocus_Filtered}"


rule Build_StatsReports:
    input:
        vcf_raw = vcf_raw,
        vcf_Locus_Filtered = outputs_directory+"/01_Locus_Filtered.vcf",
        vcf_SampleLocus_Filtered = outputs_directory+"/02_SampleLocus_Filtered.vcf",
        vcf_PopGenStatsSampleLocus_Filtered = outputs_directory+"/03_PopGenStatsSampleLocus_Filtered.vcf"
    output:
        stats_raw = VCF_reports_dir+"/00_variants_raw_vcf.stats",
        stats_Locus = VCF_reports_dir+"/01_Locus_Filtered_vcf.stats",
        stats_SampleLocus = VCF_reports_dir+"/02_SampleLocus_Filtered_vcf.stats",
        stats_PopGenStatsSampleLocus = VCF_reports_dir+"/03_PopGenStatsSampleLocus_Filtered_vcf.stats"
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "bcftools stats {input.vcf_raw} > {output.stats_raw};"
        "bcftools stats {input.vcf_Locus_Filtered} > {output.stats_Locus};"
        "bcftools stats {input.vcf_SampleLocus_Filtered} > {output.stats_SampleLocus};"
        "bcftools stats {input.vcf_PopGenStatsSampleLocus_Filtered} > {output.stats_PopGenStatsSampleLocus}"


rule Build_Report:
    input:
        VCF_reports_dir+"/00_variants_raw_vcf.stats",
        VCF_reports_dir+"/01_Locus_Filtered_vcf.stats",
        VCF_reports_dir+"/02_SampleLocus_Filtered_vcf.stats",
        VCF_reports_dir+"/03_PopGenStatsSampleLocus_Filtered_vcf.stats"
    output:
        VCF_reports_dir+"/multiQC_VcfFiltering_report.html"
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "multiqc {input} -o {VCF_reports_dir} -n multiQC_VcfFiltering_report"


rule Metadata:
    output:
        outputs_directory+"/workflow_info.txt"
    shell:
        "echo -e \"Date and time:\" > {outputs_directory}/workflow_info.txt;"
        "Date=$(date);"
        "echo -e \"${{Date}}\\n\" >> {outputs_directory}/workflow_info.txt;"
        "echo -e \"Workflow:\" >> {outputs_directory}/workflow_info.txt;"
        "echo -e \"https://github.com/BioInfo-GE2POP-BLE/CAPTURE_SNAKEMAKE_WORKFLOWS/tree/main/READS_MAPPING\\n\" >> {outputs_directory}/workflow_info.txt;"
        "cd {snakefile_dir};"
        "if git rev-parse --git-dir > /dev/null 2>&1; then echo -e \"Commit ID:\" >> {outputs_directory}/workflow_info.txt; git rev-parse HEAD >> {outputs_directory}/workflow_info.txt ; fi"
