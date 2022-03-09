#!/usr/bin/env python

import pandas as pd
import os,sys
from itertools import compress

#singularity: "docker://condaforge/mambaforge"


####################   DEFINE CONFIG VARIABLES BASED ON CONFIG FILE   ####################

### Variables from config file

bam_dir = config["BAM_DIR"] # ou par defaut si pipelin complet?
outputs_dirname = config["OUTPUTS_DIRNAME"]
reference = config["REFERENCE"]

gatk_HaplotypeCaller_java_mem = config["GATK_HAPLOTYPE_CALLER_JAVA_MEM"]
gatk_HaplotypeCaller_extra_options = config["GATK_HAPLOTYPE_CALLER_EXTRA_OPTIONS"]

gatk_GenomicsDBImport_java_mem = config["GATK_GENOMICS_DB_IMPORT_JAVA_MEM"]
gatk_GenomicsDBImport_extra_options = config["GATK_GENOMICS_DB_IMPORT_EXTRA_OPTIONS"]

gatk_Genotype_GVCFs_java_mem = config["GATK_GENOTYPE_GVCFS_JAVA_MEM"]
gatk_Genotype_GVCFs_extra_options = config["GATK_GENOTYPE_GVCFS_EXTRA_OPTIONS"]

### Define paths
path_to_snakefile = workflow.snakefile
snakefile_dir = path_to_snakefile.rsplit('/', 1)[0]
scripts_dir = snakefile_dir+"/SCRIPTS"
working_directory = os.getcwd()

### base names
#récupérer base names des .bam


## Define outputs subfolders
outputs_directory = working_directory+"/"+outputs_dirname+"/SNP_CALLING"
snp_caller_hc_dir = outputs_directory+"/HAPLOTYPE_CALLER"
snp_caller_gdbi_dir = outputs_directory+"/GENOMICS_DB_IMPORT"

### PIPELINE ###

rule FinalTargets:
    input:



 # ----------------------------------------------------------------------------------------------- #

# création de l'index de la référence si besoin
#sbatch --partition=agap_normal --job-name=faidx --wrap="module load samtools/1.10-bin ; /nfs/work/agap_id-bin/img/samtools/1.10-bin/bin/samtools faidx /home/ardissonm/scratch/TMP/GATK/REFERENCES/Svevo_baits_04_20_2019_0pident_40nident_0length.fasta"

# création du dictionnaire pour la DB
#sbatch --partition=agap_normal --mem=10G --job-name=dict --wrap=" module load jre/jre.8_x64 ; module load gatk/4.2.0.0 ; gatk CreateSequenceDictionary REFERENCE=/home/ardissonm/scratch/TMP/GATK/REFERENCES/Svevo_baits_04_20_2019_0pident_40nident_0length.fasta OUTPUT=/home/ardissonm/scratch/TMP/GATK/REFERENCES/Svevo_baits_04_20_2019_0pident_40nident_0length.dict"

# création d'une liste de contigs pour la DB
#grep "@SQ" Svevo_baits_04_20_2019_0pident_40nident_0length.dict |cut -f2,3 | sed 's/SN:\([^\t]*\)\tLN:\([0-9]*\)/\1:1-\2/'>contigs_intervals_Svevo_baits_04_20_2019_40nident_GenomicsDB_GATK.list



### HaplotypeCaller
#récupérer la base des noms de génotype sans ".bam"
#for bam in *.bam
#sbatch --partition=agap_normal --mem=10G --job-name=HC10_P2 --wrap="module load jre/jre.8_x64; module load gatk/4.2.0.0 ; gatk --java-options "-Xmx4g" HaplotypeCaller --reference /home/ardissonm/scratch/TMP/GATK/REFERENCES/Svevo_baits_04_20_2019_0pident_40nident_0length.fasta --input "$bam" --output /home/ardissonm/scratch/TMP/GATK/HAPLOTYPE_GVCF_SIMUL10_PLOIDY2/"$bam"_ploidy2.g.vcf.gz --sample-ploidy 2 -ERC GVCF"
#done

rule haplotypeCaller:
    input:
        bam_dir+"/{base}.bam"
    output:
        vcf = snp_caller_hc_dir+"/{base}.g.vcf.gz",
        tbi = snp_caller_hc_dir+"/{base}.g.vcf.gz.tbi"
    params:
        java_options = config["GATK_HAPLOTYPE_CALLER_JAVA_OPTIONS"],
        extra_options = config["GATK_HAPLOTYPE_CALLER_EXTRA_OPTIONS"]
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "gatk --java-options \"{params.java_options}\" HaplotypeCaller --reference {reference} --input {input} --output {output.vcf} {params.extra_options} -ERC GVCF"


rule create_HaplotypeList:
    input:
        expand("{snp_caller_hc_dir}/{sample}.g.vcf.gz", sample=samples, snp_caller_hc_dir=snp_caller_hc_dir)
    output:
        snp_caller_hc_dir+"/vcf.list.txt"
    shell:
        #"ls {input} | awk '{print $1\"\t\"$1}' | sed 's/\.g\.vcf.gz//' > vcf.list.txt"
        "for vcf in {input} ; do sample=$(basename ${{vcf}} .g.vcf.gz) ; echo ${{sample}}\"\t\"${{vcf}} ; done > {snp_caller_hc_dir}/vcf.list.txt"


rule GenomicsDBImport:
    input:
        vcf_list = snp_caller_hc_dir+"/vcf.list.txt"
    output:
        snp_caller_gdbi_dir
    params:
        java_options = config["GATK_GENOMICS_DB_IMPORT_JAVA_MEM"],
        extra_options = config["GATK_GENOMICS_DB_IMPORT_EXTRA_OPTIONS"]
    shell:
        "gatk --java-options \"{params.java_options}\" GenomicsDBImport --sample-name-map {input.vcf_list} --genomicsdb-workspace-path {snp_caller_gdbi_dir} --intervals /home/ardissonm/scratch/TMP/GATK/REFERENCES/contigs_intervals_Svevo_baits_04_20_2019_40nident_GenomicsDB_GATK.list  --tmp-dir {snp_caller_gdbi_dir}/tmp_file_DB"






# créer une liste des vcf
#ls Gt*.vcf.gz | awk '{print $1"\t"$1}' | sed 's/_p10_ploidy2\.g\.vcf.gz//' > vcf.list.txt

### GenomicsDBImport

#sbatch --partition=agap_normal --cpus-per-task=20 --mem-per-cpu=3G --job-name=p10_p4b --wrap="module load jre/jre.8_x64; module load gatk/4.2.0.0 ; gatk --java-options "-Xmx30g" GenomicsDBImport --sample-name-map /home/ardissonm/scratch/TMP/GATK/HAPLOTYPE_GVCF_SIMUL10_PLOIDY4/vcf.list.txt --genomicsdb-workspace-path GENOMICS_DB_SIMUL10_PLOIDY4 --intervals /home/ardissonm/scratch/TMP/GATK/REFERENCES/contigs_intervals_Svevo_baits_04_20_2019_40nident_GenomicsDB_GATK.list --merge-contigs-into-num-partitions 20 --batch-size 50 --reader-threads 20 --tmp-dir /home/ardissonm/scratch/TMP/GATK/tmp_file_DB"


### GenotypeGVCFs
#sbatch --partition=agap_normal --cpus-per-task=20 --job-name=s10p2h1 --wrap="module load jre/jre.8_x64; module load gatk/4.2.0.0 ; gatk --java-options "-Xmx30g" GenotypeGVCFs --reference /home/ardissonm/scratch/TMP/GATK/REFERENCES/Svevo_baits_04_20_2019_0pident_40nident_0length.fasta --variant gendb://GENOMICS_DB_SIMUL10_PLOIDY2 --heterozygosity 0.001 --output /home/ardissonm/scratch/TMP/GATK/HAPLOTYPE_GVCF_SIMUL10_PLOIDY2/SNPs_SIMUL10_PLOIDY2_HET0001.vcf.gz --tmp-dir /home/ardissonm/scratch/TMP/GATK/tmp_file_DB"
