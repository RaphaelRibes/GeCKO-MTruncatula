# VCF FILTERING

This VCF_FILTERING workflow allows to filter the raw vcf obtained after the variants calling step. Different types of filters can be applied: locus-based filters, a sample-based filter and filters based on population genetics statistics. 


### The VCF_FILTERING workflow's steps
1) Variants are filtered by **locus** with [vcftools](http://vcftools.sourceforge.net/man_latest.html). The filters proposed by VcfTools can be applied either at the genotype level (locus x sample) or at the locus level. For genotype level filtering, genotypes that don't pass the given thresholds (e.g. minGQ, minDP) are replaced by missing data. For locus level filtering, loci that do not pass the given thresholds (e.g. minQ, max-missing) are completely removed.
2) Variants are filtered by **samples**. After performing the previous filtering step, it is possible to remove samples with too many missing data (proportion of ungenotyped loci above a given threshold).
3) Variants are filtered by **population genetics statistics** (e.g. FIS, He). 
4) A MultiQC report is created, showing variants informations/statistics after each filtering step.


## QUICK START

To easily launch the workflow, use our runSnakemakeWorkflow.sh launcher:  
```./runSnakemakeWorkflow.sh --workflow VcfFiltering --workflow-path PATH/TO/CAPTURE_SNAKEMAKE_WORKFLOWS```  

Needed files:  
- the full CAPTURE_SNAKEMAKE_WORKFLOWS/ folder  
- the runSnakemakeWorkflow.sh launcher  
- your variants calling .vcf.gz files
- the cluster_config_VcfFiltering.yml (in case you work on a cluster) and config_VcfFiltering.yml files in a CONFIG folder  

&nbsp;


For example, if you need to launch the workflow on our VCF example dataset on a Slurm job-scheduler, run the following command from the EXAMPLE directory:  
```../../runSnakemakeWorkflow.sh --workflow VariantsCalling --workflow-path ../../../CAPTURE_SNAKEMAKE_WORKFLOWS --config-file CONFIG/config_VcfFiltering.yml --cluster-config CONFIG/cluster_config_VcfFiltering_SLURM.yml --jobs 20 --job-scheduler SLURM```  


&nbsp;

![](https://github.com/BioInfo-GE2POP-BLE/CAPTURE_PIPELINES_SNAKEMAKE/blob/main/readme_img/VcfFiltering_4elements.png)


## How to use the VCF_FILTERING workflow
 
1) [Prepare your input data](#1-prepare-your-input-data)  
2) [Clone our GitHub repository](#2-clone-our-github-repository)  
3) [Prepare the CONFIG files](#3-prepare-the-config-files)  
4) [Launch the analysis](#4-launch-the-analysis)  
5) [Expected outputs](#5-expected-outputs)


### 1/ Prepare your input data

The input data is the .vcf file obtained after the variants calling step (variants_calling.vcf.gz) and its associated index file (variants_calling.vcf.gz.tbi). 


### 2/ Clone our GitHub repository

The CAPTURE_SNAKEMAKE_WORKFLOWS folder must be fully copied in a workspace/storage of your choice.  
For example, you can clone the repository with:  
```git clone git@github.com:BioInfo-GE2POP-BLE/CAPTURE_SNAKEMAKE_WORKFLOWS.git```   


### 3/ Prepare the config files

The VCF_FILTERING workflow will need information about the dataset and the analysis parameters to perform its different steps.  
These information are provided through two files: *cluster_config_VcfFiltering.yml* and *config_VcfFiltering.yml*.  
If you name them exactly as written above and place them in a folder named 'CONFIG', the bash launching script will detect them automatically. Otherwise, you will have to pass them as arguments with --config and --cluster-config (see [below](#4-launch-the-analysis) for details).

#### *A/ The cluster_config_VcfFiltering.yml file:*
This file will be needed if you run the workflow on a computer cluster and want Snakemake to submit jobs. You <ins>only need to modify the partitions or queues names</ins> to match those of your cluster. The first section of the file gives the default values for the job-scheduler's parameters that Snakemake should use for all its steps (or rules). The following sections correspond to specific Snakemake steps, with new parameters values to overwrite the defaults. If you want to assign a different partition/queue for a specific step that does not yet have its own section, you can create a new section for it:  

	specificStepName:
    	q or partition: {partitionNameForSpecificStep}

You will find [the list of the steps names](#list-of-the-snakefile-rules) along with what they do and the tools they use at the end of this page.  
Our workflows support SGE and Slurm job-schedulers. <ins>You will find cluster-config files for both in the EXAMPLE/CONFIG folder</ins>.  


#### *B/ The config_VcfFiltering.yml file:*  
This file is used to pass all the information and tools parameters that will be used by the VCF_FILTERING workflow. The workflow expects it to contain a specific list of variables and their assigned values, organized in YAML format. Expected variables are:  

**INPUT FILES**  
- *VCF_FILE*&nbsp;&nbsp;&nbsp;The path to the variants calling file in zipped vcf format (.vcf.gz).

**VCF FILTERING PARAMETERS**  

- *VCFTOOLS_LOCUS_FILTERING_OPTIONS:*&nbsp;&nbsp;&nbsp; The list of filtering options to pass to [vcftools](http://vcftools.sourceforge.net/man_latest.html) (e.g. "--minDP 5 --minQ 30 --max-missing 0.5"). Be careful to provide them between quotes. 
- *MAX_RATIO_NA_PER_SAMPLE:*&nbsp;&nbsp;&nbsp; The maximum proportion of allowed missing data per sample ("1": no filtering).
- *BCFTOOLS_FILTERING_OPTIONS:*&nbsp;&nbsp;&nbsp; The list of filtering options to pass to the 'bcftools filter' command. 13 parameters are calculated for each locus and can be used for filtering: 
	- SNP = Variant with at least 2 mononucleotide alleles (Type=Integer)
	- INDEL = Variant with at least 2 alleles of different lengths (Type=Integer)
	- nbGS = Number of genotyped samples (Type=Integer)
	- pNA = Frequency of missing genotypes for this site (Type=Float)
	- nbG = Number of genotypes (Type=Integer)
	- nbAll = Number of alleles (Type=Integer)
	- All = list of alleles (Type=String)
	- AllFreq = Frequency of all alleles (Type=Float)
	- MinHomo = Number of the least frequent homozygous genotype (Type=Integer)
	- MAF = Minor allele frequency (Type=Float)
	- MAFb = Minor base frequency (Type=Float)
	- He = Nei expected Heterozygosity (Type=Float)
	- Fis = Inbreeding coefficient (Type=Float)

  e.g. : "INFO/SNP==1 & INFO/INDEL==0 & INFO/nbAll=2 & INFO/pNA<0.5 & MinHomo>=1 & Fis>=0.8"
  Be careful to use quotes when passing this parameter.

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
When called for the first time, the VCF_FILTERING Snakemake workflow will download the tools' packages in a pkgs_dirs folder, and install them in a conda environment that will be stored in a .snakemake/conda folder, in the directory you called the workflow from. Every time you call the workflow from a new directory, the Conda environment will be generated again. To avoid creating the environment multiple times, which can be both time and resource-consuming, you can provide a specific folder where you want Snakemake to store all of its conda environments with the --conda-env-path option of the runSnakemakeWorkflow.sh launcher.  

The pkgs_dirs folder is by default common to your whole system or cluster personnal environment. Conda's standard behaviour is to create it in your home directory, in a .conda folder. If your home space is limited or if you do not have the right to write there from your cluster's nodes, you will need to tell Conda to store its packages somewhere else, thanks to a .condarc file. Place it in your home folder and specify the directory path where you want Conda to store the packages, following this example:  
```
envs_dirs:  
    - /home/username/path/to/appropriate/folder/env  
pkgs_dirs:  
    - /home/username/path/to/appropriate/folder/pkgs  
```




### 5/ Expected outputs  
This workflow will create a "VCF_FILTERING" directory in the "WORKFLOWS_OUTPUTS" directory. This directory is structured as follows and contains:  

<img src="https://github.com/BioInfo-GE2POP-BLE/CAPTURE_PIPELINES_SNAKEMAKE/blob/main/readme_img/OutputsTree_VcfFiltering.png" width="600"/>



<ins>Description of the main files:</ins>  

- *01_Locus_Filtered.vcf*:&nbsp;&nbsp;&nbsp;the vcf file after filtering variants by **locus**  
- *02_SampleLocus_Filtered.vcf*:&nbsp;&nbsp;&nbsp;the vcf file after filtering variants by **sample** 
- *samples_to_remove.list*:&nbsp;&nbsp;&nbsp;list of samples that were deleted in step 02 (filtering by sample)
- *02_SampleLocus_Filtered_withPopGenStats.vcf*:&nbsp;&nbsp;&nbsp;the intermediate vcf file corresponding to variants filtered after step 02 (locus + sample), with population genetics statistics (Fis, He, ...) by variants.  
- *03_PopGenStatsSampleLocus_Filtered.vcf*:&nbsp;&nbsp;&nbsp;the vcf file after filtering variants by **population genetics statistics**  
- *workflow_info.txt*:&nbsp;&nbsp;&nbsp;File that contains the date and time of the workflow launch, the link to the Github repository and the commit ID

**REPORTS directory** 
- *00_variants_raw_vcf.stats*:&nbsp;&nbsp;&nbsp;bcftools statistics of the unfiltered vcf file 
- *01_Locus_Filtered_vcf.stats*:&nbsp;&nbsp;&nbsp;bcftools statistics after filtering by **locus**
- *02_SampleLocus_Filtered_vcf.stats*:&nbsp;&nbsp;&nbsp;bcftools statistics after filtering by **sample**
- *03_PopGenStatsSampleLocus_Filtered_vcf.stats*:&nbsp;&nbsp;&nbsp;bcftools statistics after filtering by **population genetics statistics**
- *multiQC_VcfFiltering_report.html*:&nbsp;&nbsp;&nbsp;graphic representation of the variants informations/statistics after each filtering step


## Tools
This workflow uses the following tools: 
- [vcftools 0.1.16](https://github.com/vcftools/vcftools)
- [bcftools 1.15](https://samtools.github.io/bcftools/bcftools.html)
- [multiqc v1.11](https://github.com/ewels/MultiQC/releases)

These tools are loaded in a CONDA environment from the conda-forge and bioconda channels.

##  List of the snakefile rules
Name, description and tools used for each of the snakemake workflow rules:

| **Rule name**          | **Description**                                                                                            | **Tools**       |
|:----------------------:|:----------------------------------------------------------------------------------------------------------:|:---------------:|
| Filter_Loci            | Filtering variants by locus. Applied either at the genotype level (locus X sample) or at the locus level   | vcftools        |
| Filter_Samples         | Filtering variants by sample. Removing samples that have too many loci with missing data.                  | vcftools        |
| Calculate_PopGenStats  | Calculating population genetics statistics (e.g. FIS, He)                                                  |                 |
| Filter_PopGenStats     | Filtering variants based on population genetics statistics (e.g. FIS, He)                                  | bcftools filter |
| Build_StatsReport      | Building statistics reports for unfiltered variants and each filtering step                                | bcftools stats  |
| Build_Report           | Running MultiQC on unfiltered variants and each filtering step                                             | MultiQC         |



![](https://github.com/BioInfo-GE2POP-BLE/CAPTURE_PIPELINES_SNAKEMAKE/blob/main/readme_img/VcfFiltering_Workflow.jpg?raw=true)


