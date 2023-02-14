#!/usr/bin/env python

import os,sys
from itertools import compress


####################   DEFINE CONFIG VARIABLES BASED ON CONFIG FILE   ####################

### Variables from config file
vcf_raw = config["VCF_FILE"]
filtering_subfolder = ""
if (len(config["FILTERING_SUBFOLDER"]) > 0):
    filtering_subfolder = "/"+config["FILTERING_SUBFOLDER"]


### Define paths
path_to_snakefile = workflow.snakefile
snakefile_dir = path_to_snakefile.rsplit('/', 1)[0]
scripts_dir = snakefile_dir+"/SCRIPTS"
working_directory = os.getcwd()

### Define outputs subfolders
outputs_directory = working_directory+"/WORKFLOWS_OUTPUTS/VCF_FILTERING"+filtering_subfolder
VCF_reports_dir = outputs_directory+"/REPORTS"

### PIPELINE ###

rule FinalTargets:
    input:
        VCF_reports_dir+"/variants_stats_histograms_VF.pdf",
        VCF_reports_dir+"/genotypes_DP_boxplot_VF.pdf",
        VCF_reports_dir+"/variants_along_genome_VF.pdf",
        VCF_reports_dir+"/multiQC_VcfFiltering_report.html",
        VCF_reports_dir+"/missing_data_per_sample.txt",
        outputs_directory+"/workflow_info.txt"



 # ----------------------------------------------------------------------------------------------- #



rule Filter_Genotypes:
    input:
        vcf_raw
    output:
        outputs_directory+"/01__Genotype_Filtered.vcf"
    conda:
        "ENVS/conda_tools.yml"
    params:
        config["BCFTOOLS_GENOTYPES_FILTERING_OPTIONS"]
    shell:
        "bcftools filter -Ou -i '{params}' -S . {input} | bcftools view -Ou --exclude-uncalled --trim-alt-alleles | bcftools view -m2 -o {output}"



rule Filter_Loci_1:
    input:
        outputs_directory+"/01__Genotype_Filtered.vcf"
    output:
        #tmp_Locus_Filtered = temp(outputs_directory+"/tmp_Locus_Filtered.vcf.gz"),
        outputs_directory+"/02__Genotype_Locus1_Filtered.vcf"
    conda:
        "ENVS/conda_tools.yml"
    params:
        locus_filters = config["BCFTOOLS_LOCUS_FILTERING1_OPTIONS"]
    shell:
        "bcftools filter -Ou -sFilter1 -i '{params}' {input} | bcftools view -f 'PASS' > {output}"



rule Filter_Samples:
    input:
        outputs_directory+"/02__Genotype_Locus1_Filtered.vcf"
    output:
        out_imiss = temp(outputs_directory+"/SampleFilter.imiss"),
        samples_to_remove = outputs_directory+"/samples_to_remove.list",
        SampleLocus_Filtered = outputs_directory+"/03__Genotype_Locus1_Sample_Filtered.vcf"
    conda:
        "ENVS/conda_tools.yml"
    params:
        config["MAX_NA_PER_SAMPLE"]
    shell:
        "paste <(bcftools query -f '[%SAMPLE\t]\n' {input} | head -1 | tr '\t' '\n')"
        " <(bcftools query -f '[%GT\t]\n' {input} | awk -v OFS=\"\t\" '{{for (i=1;i<=NF;i++) if (($i == \"./.\") || ($i == \".|.\")) sum[i]+=1 }} END {{for (i in sum) print i, sum[i] / NR }}' | sort -k1,1n | cut -f 2)"
        " > {output.out_imiss};"
        "awk -v max_pc_NA={params} '{{if($2>max_pc_NA){{print $1}}}}' {output.out_imiss} > {output.samples_to_remove};"
        "bcftools view -Ou -S ^{output.samples_to_remove} {input} | bcftools view -Ou --exclude-uncalled --trim-alt-alleles | bcftools view -m2 -o {output.SampleLocus_Filtered}"



rule Calculate_LocusExtraStats:
    input:
        outputs_directory+"/03__Genotype_Locus1_Sample_Filtered.vcf"
    output:
        outputs_directory+"/03__Genotype_Locus1_Sample_Filtered__withExtraStats.vcf"
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "sed -i 's/=nan/=./g' {input};"
        "module purge ; python {scripts_dir}/egglib_LocusExtraStats.py --input {input} --output {output}"


rule Filter_Loci_2:
    input:
        outputs_directory+"/03__Genotype_Locus1_Sample_Filtered__withExtraStats.vcf"
    output:
        outputs_directory+"/04__Genotype_Locus1_Sample_Locus2_Filtered.vcf"
    conda:
        "ENVS/conda_tools.yml"
    params:
        config["BCFTOOLS_LOCUS_FILTERING2_OPTIONS"]
    shell:
        "bcftools filter -Ou -sFilter2 -i '{params}' {input} | bcftools view -f 'PASS' > {output}"


rule Build_StatsReports:
    input:
        vcf_raw = vcf_raw,
        vcf_Genotype_Filtered = outputs_directory+"/01__Genotype_Filtered.vcf",
        vcf_GenotypeLocus1_Filtered = outputs_directory+"/02__Genotype_Locus1_Filtered.vcf",
        vcf_GenotypeLocus1Sample_Filtered = outputs_directory+"/03__Genotype_Locus1_Sample_Filtered.vcf",
        vcf_GenotypeLocus1SampleLocus2_Filtered = outputs_directory+"/04__Genotype_Locus1_Sample_Locus2_Filtered.vcf"
    output:
        stats_raw = VCF_reports_dir+"/00__Raw_Variants.stats",
        stats_Genotype = VCF_reports_dir+"/01__Genotype_Filtered.stats",
        stats_GenotypeLocus1 = VCF_reports_dir+"/02__Genotype_Locus1_Filtered.stats",
        stats_GenotypeLocus1Sample = VCF_reports_dir+"/03__Genotype_Locus1_Sample_Filtered.stats",
        stats_GenotypeLocus1SampleLocus2 = VCF_reports_dir+"/04__Genotype_Locus1_Sample_Locus2_Filtered.stats"
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "bcftools stats {input.vcf_raw} > {output.stats_raw};"
        "bcftools stats {input.vcf_Genotype_Filtered} > {output.stats_Genotype};"
        "bcftools stats {input.vcf_GenotypeLocus1_Filtered} > {output.stats_GenotypeLocus1};"
        "bcftools stats {input.vcf_GenotypeLocus1Sample_Filtered} > {output.stats_GenotypeLocus1Sample};"
        "bcftools stats {input.vcf_GenotypeLocus1SampleLocus2_Filtered} > {output.stats_GenotypeLocus1SampleLocus2}"


rule Build_Report:
    input:
        VCF_reports_dir+"/00__Raw_Variants.stats",
        VCF_reports_dir+"/01__Genotype_Filtered.stats",
        VCF_reports_dir+"/02__Genotype_Locus1_Filtered.stats",
        VCF_reports_dir+"/03__Genotype_Locus1_Sample_Filtered.stats",
        VCF_reports_dir+"/04__Genotype_Locus1_Sample_Locus2_Filtered.stats"
    output:
        VCF_reports_dir+"/multiQC_VcfFiltering_report.html"
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "multiqc {input} -o {VCF_reports_dir} -n multiQC_VcfFiltering_report -i VcfFiltering_report"


rule Summarize_FinalVCFVariables:
    input:
        outputs_directory+"/04__Genotype_Locus1_Sample_Locus2_Filtered.vcf"
    output:
        stats_tsv = VCF_reports_dir+"/variants_stats_VF.tsv",
        DP_tsv = temp(VCF_reports_dir+"/genotypes_DP_VF.tsv"),
        GT_tsv = temp(VCF_reports_dir+"/genotypes_GT_VF.tsv"),
        pos_tsv = temp(VCF_reports_dir+"/variants_pos.tsv"),
        lengths_tsv = temp(VCF_reports_dir+"/contigs_lengths.tsv")
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "{scripts_dir}/extract_variants_stats_from_vcf.sh {input} {output.stats_tsv} {output.DP_tsv} {output.GT_tsv} {output.pos_tsv} {output.lengths_tsv} {VCF_reports_dir}"


rule Plot_FinalVCFVariablesHistograms:
    input:
        VCF_reports_dir+"/variants_stats_VF.tsv"
    output:
        VCF_reports_dir+"/variants_stats_histograms_VF.pdf"
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "python {scripts_dir}/plot_variants_stats_histograms.py --input {input} --output {output}"


rule Plot_FinalVCFDPBoxplot:
    input:
        DP_tsv = VCF_reports_dir+"/genotypes_DP_VF.tsv",
        GT_tsv = VCF_reports_dir+"/genotypes_GT_VF.tsv"
    output:
        VCF_reports_dir+"/genotypes_DP_boxplot_VF.pdf"
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "python {scripts_dir}/plot_DP_boxplot.py --input-DP {input.DP_tsv} --input-GT {input.GT_tsv} --output {output}"


rule Plot_FinalVCFVariantsAlongGenome:
    input:
        pos_tsv = VCF_reports_dir+"/variants_pos.tsv",
        lengths_tsv = VCF_reports_dir+"/contigs_lengths.tsv"
    output:
        VCF_reports_dir+"/variants_along_genome_VF.pdf"
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "python {scripts_dir}/plot_variants_along_genome.py --snp-pos {input.pos_tsv} --contigs-lengths {input.lengths_tsv} --output {output}"


rule Compute_MissingDataPerSample:
    input:
        outputs_directory+"/04__Genotype_Locus1_Sample_Locus2_Filtered.vcf"
    output:
        VCF_reports_dir+"/missing_data_per_sample.txt"
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "echo -e \"Sample\tNA_fraction\" > {output};"
        "paste <(bcftools query -f '[%SAMPLE\t]\n' {input} | head -1 | tr '\t' '\n')"
        " <(bcftools query -f '[%GT\t]\n' {input} | awk -v OFS=\"\t\" '{{for (i=1;i<=NF;i++) if (($i == \"./.\") || ($i == \".|.\")) sum[i]+=1 }} END {{for (i in sum) print i, sum[i] / NR }}' | sort -k1,1n | cut -f 2)"
        " | awk 'NF' >> {output};"

rule Metadata:
    output:
        outputs_directory+"/workflow_info.txt"
    shell:
        "echo -e \"Date and time:\" > {outputs_directory}/workflow_info.txt;"
        "Date=$(date);"
        "echo -e \"${{Date}}\\n\" >> {outputs_directory}/workflow_info.txt;"
        "echo -e \"Workflow:\" >> {outputs_directory}/workflow_info.txt;"
        "echo -e \"https://github.com/GE2POP/GeCKO/tree/main/READS_MAPPING\\n\" >> {outputs_directory}/workflow_info.txt;"
        "cd {snakefile_dir};"
        "if git rev-parse --git-dir > /dev/null 2>&1; then echo -e \"Commit ID:\" >> {outputs_directory}/workflow_info.txt; git rev-parse HEAD >> {outputs_directory}/workflow_info.txt ; fi"
