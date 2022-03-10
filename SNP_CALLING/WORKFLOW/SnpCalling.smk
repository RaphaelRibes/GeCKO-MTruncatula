#!/usr/bin/env python

import pandas as pd
import os,sys
from itertools import compress

#singularity: "docker://condaforge/mambaforge"


### Variables from config file
outputs_dirname = config["OUTPUTS_DIRNAME"]
bam_dir = config["BAM_DIR"] # ou mapping_dir ? > comment lui dire de prendre par defaut si pipeline complet?

reference = config["REFERENCE"]

### Samples list
bams_list = [bam for bam in os.listdir(bam_dir) if bam.endswith(".bam")]
samples = []
for bam in bams_list:
    sample = bam.replace(".bam", "")
    samples.append(sample)

### remove .fa .fas .fasta file extension
reference_base = reference.rsplit('.fa', 1)[0]

### Define paths
path_to_snakefile = workflow.snakefile
snakefile_dir = path_to_snakefile.rsplit('/', 1)[0]
scripts_dir = snakefile_dir+"/SCRIPTS"
working_directory = os.getcwd()

## Define outputs subfolders
outputs_directory = working_directory+"/"+outputs_dirname+"/SNP_CALLING"
HaplotypeCaller_dir = outputs_directory+"/HAPLOTYPE_CALLER"
GenomicsDBImport_dir = outputs_directory+"/GENOMICS_DB_IMPORT"
GenotypeGVCFs_dir = outputs_directory+"/GENOTPYE_GVCFS"

 # ----------------------------------------------------------------------------------------------- #

### PIPELINE ###

rule FinalTargets:
    input:
        #reference_base+".fai",
        #reference_base+".dict",
        #reference_base+"_intervals_for_GATK.list",
        #HaplotypeCaller_dir+"/{base}.g.vcf.gz",
        #HaplotypeCaller_dir+"/{base}.g.vcf.gz.tbi",
        #HaplotypeCaller_dir+"/vcf.list.txt",
        GenotypeGVCFs_dir+"/variants.vcf.gz"



# création de l'index de la référence si besoin
#sbatch --partition=agap_normal --job-name=faidx --wrap="module load samtools/1.10-bin ; /nfs/work/agap_id-bin/img/samtools/1.10-bin/bin/samtools faidx /home/ardissonm/scratch/TMP/GATK/REFERENCES/Svevo_baits_04_20_2019_0pident_40nident_0length.fasta"

rule Index_Reference:
    input:
        reference
    output:
        reference+".fai"
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "samtools faidx {input}"

# création du dictionnaire pour la DB
#sbatch --partition=agap_normal --mem=10G --job-name=dict --wrap=" module load jre/jre.8_x64 ; module load gatk/4.2.0.0 ; gatk CreateSequenceDictionary REFERENCE=/home/ardissonm/scratch/TMP/GATK/REFERENCES/Svevo_baits_04_20_2019_0pident_40nident_0length.fasta OUTPUT=/home/ardissonm/scratch/TMP/GATK/REFERENCES/Svevo_baits_04_20_2019_0pident_40nident_0length.dict"

rule Dictionary_Reference:
    input:
        reference # .fai?
    output:
        reference_base+".dict"
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "gatk CreateSequenceDictionary REFERENCE={input} OUTPUT={output}"


# création d'une liste de contigs pour la DB
#grep "@SQ" Svevo_baits_04_20_2019_0pident_40nident_0length.dict |cut -f2,3 | sed 's/SN:\([^\t]*\)\tLN:\([0-9]*\)/\1:1-\2/'>contigs_intervals_Svevo_baits_04_20_2019_40nident_GenomicsDB_GATK.list

rule ListIntervalsReference_Dictionary:
    input:
        reference_base+".dict"
    output:
        reference_base+"_intervals_for_GATK.list"
    shell:
        "grep '@SQ' {input} | cut -f2,3 | sed 's/SN://' | sed 's/LN://' | awk '{{print $1\":1-\"$2}}' > {output}"


### HaplotypeCaller
#récupérer la base des noms de génotype sans ".bam"
#for bam in *.bam
#sbatch --partition=agap_normal --mem=10G --job-name=HC10_P2 --wrap="module load jre/jre.8_x64; module load gatk/4.2.0.0 ; gatk --java-options "-Xmx4g" HaplotypeCaller --reference /home/ardissonm/scratch/TMP/GATK/REFERENCES/Svevo_baits_04_20_2019_0pident_40nident_0length.fasta --input "$bam" --output /home/ardissonm/scratch/TMP/GATK/HAPLOTYPE_GVCF_SIMUL10_PLOIDY2/"$bam"_ploidy2.g.vcf.gz --sample-ploidy 2 -ERC GVCF"
#done

rule HaplotypeCaller:
    input:
        reference = reference,
        reference+".fai",
        bams = bam_dir+"/{base}.bam"
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


# créer une liste des vcf
#ls Gt*.vcf.gz | awk '{print $1"\t"$1}' | sed 's/_p10_ploidy2\.g\.vcf.gz//' > vcf.list.txt

rule List_Haplotype:
    input:
        expand("{HaplotypeCaller_dir}/{sample}.g.vcf.gz", sample=samples, HaplotypeCaller_dir=HaplotypeCaller_dir)
    output:
        HaplotypeCaller_dir+"/vcf.list.txt"
    shell:
        #"ls {input} | awk '{print $1\"\t\"$1}' | sed 's/\.g\.vcf.gz//' > vcf.list.txt"
        "for vcf in {input} ; do sample=$(basename ${{vcf}} .g.vcf.gz) ; echo ${{sample}}\"\t\"${{vcf}} ; done > {HaplotypeCaller_dir}/vcf.list.txt"


### GenomicsDBImport

#sbatch --partition=agap_normal --cpus-per-task=20 --mem-per-cpu=3G --job-name=p10_p4b --wrap="module load jre/jre.8_x64; module load gatk/4.2.0.0 ; gatk --java-options "-Xmx30g" GenomicsDBImport --sample-name-map /home/ardissonm/scratch/TMP/GATK/HAPLOTYPE_GVCF_SIMUL10_PLOIDY4/vcf.list.txt --genomicsdb-workspace-path GENOMICS_DB_SIMUL10_PLOIDY4 --intervals /home/ardissonm/scratch/TMP/GATK/REFERENCES/contigs_intervals_Svevo_baits_04_20_2019_40nident_GenomicsDB_GATK.list --merge-contigs-into-num-partitions 20 --batch-size 50 --reader-threads 20 --tmp-dir /home/ardissonm/scratch/TMP/GATK/tmp_file_DB"

rule GenomicsDBImport:
    input:
        vcf_list = HaplotypeCaller_dir+"/vcf.list.txt",
        intervals = reference_base+"_intervals_for_GATK.list"
    output:
        GenomicsDBImport_dir,
        outputs_directory+"/tmp_file"
    params:
        java_options = config["GATK_GENOMICS_DB_IMPORT_JAVA_OPTIONS"],
        extra_options = config["GATK_GENOMICS_DB_IMPORT_EXTRA_OPTIONS"]
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "mkdir -p {outputs_directory}/tmp_file;"
        "gatk --java-options \"{params.java_options}\" GenomicsDBImport --sample-name-map {input.vcf_list} --intervals {input.intervals} {params.extra_options} --genomicsdb-workspace-path {GenomicsDBImport_dir} --tmp-dir {outputs_directory}/tmp_file"


### GenotypeGVCFs
#sbatch --partition=agap_normal --cpus-per-task=20 --job-name=s10p2h1 --wrap="module load jre/jre.8_x64; module load gatk/4.2.0.0 ; gatk --java-options "-Xmx30g" GenotypeGVCFs --reference /home/ardissonm/scratch/TMP/GATK/REFERENCES/Svevo_baits_04_20_2019_0pident_40nident_0length.fasta --variant gendb://GENOMICS_DB_SIMUL10_PLOIDY2 --heterozygosity 0.001 --output /home/ardissonm/scratch/TMP/GATK/HAPLOTYPE_GVCF_SIMUL10_PLOIDY2/SNPs_SIMUL10_PLOIDY2_HET0001.vcf.gz --tmp-dir /home/ardissonm/scratch/TMP/GATK/tmp_file_DB"

rule GenotypeGVCFs:
    input:
        reference,
        GenomicsDBImport_dir
    output:
        GenotypeGVCFs_dir+"/variants.vcf.gz"
    params:
        java_options = config["GATK_GENOTYPE_GVCFS_JAVA_OPTIONS"],
        extra_options = config["GATK_GENOTYPE_GVCFS_EXTRA_OPTIONS"]
    conda:
        "ENVS/conda_tools.yml"
    shell:
        "gatk --java-options \"{params.java_options}\" GenotypeGVCFs --reference {reference} --variant gendb://{GenomicsDBImport_dir} {params.extra_options} --output {output} --tmp-dir {outputs_directory}/tmp_file"
