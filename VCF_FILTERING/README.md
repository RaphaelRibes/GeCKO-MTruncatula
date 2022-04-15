# VCF FILTERING

this VCF_FILTERING workflow allows to filter the raw vcf obtained after the variants calling. Different types of filters are applied: locus-based filters, a sample-based filter and filters based on population genetics statistics. 


### The VCF_FILTERING workflow's steps
1) Variants filtering by **locus** with [vcftools](http://vcftools.sourceforge.net/man_latest.html). The filters proposed by VcfTools can be applied either at the genotype level (locus x samples) and in this case the data not respecting the filters are replaced by missing data (e.g. minGQ, minDP). or at the locus level and in this case the locus is completely deleted (e.g. minQ, max-missing)
2) Variants filtering by **samples**. After performing the first filter step (locus), it is possible to remove samples that have too many locus with missing data (proportion).
3) Variants filtering by **population genetics statistics** (e.g. FIS, He) 
4) MultiQC report is created, showing variants informations/statistics after each filtering steps.


## QUICK START

To easily launch the workflow, use our runSnakemakeWorkflow.sh launcher:  
```./runSnakemakeWorkflow.sh --workflow VcfFiltering --workflow-path PATH/TO/CAPTURE_SNAKEMAKE_WORKFLOWS```  

Needed files:  
- the full CAPTURE_SNAKEMAKE_WORKFLOWS/ folder  
- the runSnakemakeWorkflow.sh launcher  
- your variants calling .vcf.gz files
- the cluster_config_VcfFiltering.yml (in case you work on a cluster) and config_VcfFiltering.yml files in a CONFIG folder  

&nbsp;

For example, if you need to launch the workflow on our ... dataset on a Slurm job-scheduler, run the following command from the EXAMPLE/... directory:  
```./runSnakemakeWorkflow.sh --workflow VariantsCalling --workflow-path /storage/replicated/cirad_users/ardissonm/CAPTURE_SNAKEMAKE_WORKFLOWS --config-file CONFIG/config_VcfFiltering.yml --cluster-config CONFIG/cluster_config_VcfFiltering.yml --jobs 20 --job-scheduler SLURM```  


&nbsp;

![](https://github.com/BioInfo-GE2POP-BLE/CAPTURE\_PIPELINES\_SNAKEMAKE/blob/main/readme\_img/VcfFiltering\_4elements.png")


## How to use the VCF_FILTERING workflow
 
1) [Prepare your input data](#1-prepare-your-input-data)  
2) [Clone our GitHub repository](#2-clone-our-github-repository)  
3) [Prepare the CONFIG files](#3-prepare-the-config-files)  
4) [Launch the analysis](#4-launch-the-analysis)  
5) [Expected outputs](#5-expected-outputs)


### 1/ Prepare your input data

the input data is the .vcf file obtained after the variants calling ( variants.vcf.gz) and the associated index files (variants_calling.vcf.gz.tbi). 


### 2/ Clone our GitHub repository

The CAPTURE_SNAKEMAKE_WORKFLOWS folder must be fully copied in a workspace/storage of your choice.  
For example, you can clone the our repository with:  
```git clone git@github.com:BioInfo-GE2POP-BLE/CAPTURE_PIPELINES_SNAKEMAKE.git```   


### 3/ Prepare the config files

The VCF_FILTERING workflow will need information about the dataset and the analysis parameters to perform its different steps.  
These information are provided through two files: *cluster_config_VcfFiltering.yml* and *config_VcfFiltering.yml*.  
If you name them exactly as written above and place them in a folder named 'CONFIG', the bash launching script will detect them automatically. Otherwise, you will have to pass them as arguments with --config and --cluster-config (see [below](#4-launch-the-analysis) for details).

#### *cluster_config_VcfFiltering.yml file:*
This file will be needed if you run the workflow on a computer cluster and want Snakemake to submit jobs. You <ins>only need to modify the partitions or queues names</ins> to match those of your cluster. The first section of the file gives the default values for the job-scheduler's parameters that Snakemake should use for all its steps (or rules). The following sections correspond to specific Snakemake steps, with new parameters values to overwrite the defaults. If you want to assign a different partition/queue for a specific step that does not yet have its own section, you can create a new section for it, preceded by a comma:  

	specificStepName:
    	q or partition: {partitionNameForSpecificStep}

You will find [the list of the steps names](#list-of-the-snakefile-rules) along with what they do and the tools they use at the end of this page.  
Our workflows support SGE and Slurm job-schedulers. <ins>You will find cluster-config files for both in the EXAMPLE/CONFIG folder</ins>.  


#### *config_VcfFiltering.yml file:*  
This file is used to pass all the information and tools parameters that will be used by the READS_MAPPING workflow. The workflow expects it to contain a specific list of variables and their assigned values, organized in YAML format. Expected variables are:  

**GENERAL VARIABLES**  

**INPUT FILES**  
- *VCF_FILE*&nbsp;&nbsp;&nbsp;The path to the variants calling .vcf.gz file.

**VCF FILTERING PARAMETERS**  

- *VCFTOOLS_LOCUS_FILTERING_OPTIONS:*&nbsp;&nbsp;&nbsp; options for [vcftools](http://vcftools.sourceforge.net/man_latest.html) filtering (e.g. "--minDP 5 --minQ 30 --max-missing 0.5" ) 
- *MAX_RATIO_NA_PER_SAMPLE:*&nbsp;&nbsp;&nbsp; maximum proportion of missing data per sample ("1": no filtering)
- *POPGENSTATS_FILTERING_OPTIONS:*&nbsp;&nbsp;&nbsp; Filtration available on 9 calculated parameters: He, heterozygote deficit (F), homozygous allele 1 (A1A1),  homozygous allele 2 (A2A2), heterozygous (A1A2), nbre of sample genotyped(nbG), NA proportion (pcNA), frequency allele 1 (p), frequency allele 1 (q). (e.g. "INFO/F>=0.8 & INFO/A1A1>0 & INFO/A2A2>0 & INFO/pcNA<0.5" )

&nbsp;

<ins>An example of config_VcfFiltering.yml file can be found in the EXAMPLE/CONFIG folder</ins>.  

&nbsp;

### 4/ Launch the analysis

**Environment**  
You can run this workflow on a computer or on a computer cluster. You will need Snakemake and Conda to be available.

**Launching**  
To launch the VCF_FILTERING workflow, you can use our launching script runSnakemakeWorkflow.sh with the option --workflow VcfFiltering:  
```./runSnakemakeWorkflow.sh --workflow VcfFiltering --workflow-path PATH/TO/CAPTURE_SNAKEMAKE_WORKFLOWS```  

For more help on how to use it, see our GitHub's general README file or run:  
```./runSnakemakeWorkflow.sh --help --workflow-path PATH/TO/CAPTURE_SNAKEMAKE_WORKFLOWS```  

**Notes on Conda**  
The workflow will download and make available the [tools it needs](#tools) through Conda, which means you do not need to have them installed in your working environment behorehand.  
When called for the first time, the VARIANTS_CALLING Snakemake workflow will download the tools' packages in a pkgs_dirs folder, and install them in a conda environment that will be stored in a .snakemake/conda folder, in the directory you called the workflow from. Every time you call the workflow from a new directory, the Conda environment will be generated again.  

The pkgs_dirs folder however is common to your whole system or cluster personnal environment. Conda's default behaviour is to create it in your home directory, in a .conda folder. If your home space is limited or if you do not have the right to write there from your cluster's nodes, you will need to tell Conda to store its packages somewhere else, thanks to a .condarc file. Place it in your home folder and specify the directory path where you want Conda to store the packages, following this example:  
```
envs_dirs:  
    - /home/username/path/to/appropriate/folder/env  
pkgs_dirs:  
    - /home/username/path/to/appropriate/folder/pkgs  
```




### 5/ Expected outputs  
the workflow create a directory "VCF_FILTERING" in the "WORKFLOWS_OUTPUTS" directory. This directory contained:
- the vcf file after Variants filtering by **locus** : 01_Locus_Filtered.recode.vcf 
- the vcf file after Variants filtering by **sample** : 02_SamplesLocus_Filtered.recode.vcf 
- vcf.gz and vcf.gz.tbi files with population genetics statistics (F, He, ...)
- the vcf file after Variants filtering by **population genetics statistics** : 03_PopGenStatsSamplesLocus_Filtered.vcf  
And REPORTS directory cotained:
- 00_variants_raw_vcf.stats: bcftools statistics on vcf file unfiltered. 
- 01_Locus_Filtered_vcf.stats: bcftools statistics after filtering by **locus**
- 02_SamplesLocus_Filtered_vcf.stats: bcftools statistics after filtering by **sample**
- 03_PopGenStatsSamplesLocus_Filtered_vcf.stats: bcftools statistics after filtering by **population genetics statistics**
- multiQC_VcfFiltering_report.html allowing to visualize informations/statistics variants after each filtering steps.

## Tools
This workflow uses the following tools: 
- [vcftools 0.1.16](https://github.com/vcftools/vcftools)
- [bcftools 1.15](https://samtools.github.io/bcftools/bcftools.html)
- [multiqc v1.11](https://github.com/ewels/MultiQC/releases)

These tools are loaded in a CONDA environment from the conda-forge and bioconda channels.

##  List of the snakefile rules
Name, description and tools used for each of the snakemake workflow rules:

| **Rule name**        | **Description**                                                                                            | **Tools**       |
|:--------------------:|:----------------------------------------------------------------------------------------------------------:|:---------------:|
| Filters_Locus        | Filtering of variants by locus. Applied either at the genotype level (locusXsamples) or at the locus level | vcftools        |
| Filter_samples       | Filtering of variants by samples. Remove samples that have too many locus with missing data.               | vcftools        |
| Calculate_PopGenStat | Calculating population genetics statistics (e.g. FIS, He)                                                  |                 |
| Filters_PopGenStat   | Filtering of variants by population genetics statistics (e.g. FIS, He)                                     | bcftools filter |
| BuildStatReport      | Building statistics reports of unfiltering variants and each filtering step                                | bcftools stats  |
| BuildReport          | Runing MultiQC on unfiltering variants and each filtering steps                                            | MultiQC         |

