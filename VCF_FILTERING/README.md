# VCF FILTERING

This VCF_FILTERING workflow allows to filter the raw VCF file obtained after the variants calling step. Different types of filters can be applied sequentially: genotype level filters, site level filters, and a sample level filter. Several additional site level statistics will be computed (including some relating to population genetics) and can be used in a final step for further loci filtering.


### The VCF_FILTERING workflow's steps
1) The first filter applies to **genotypes** (Locus x Sample). What we call ‘genotype’ here is the single-locus genotype assigned to a sample at a given position. Genotypes that don't pass the given thresholds (e.g. FMT/GQ>=15; FMT/DP>=5) are replaced by missing data. Loci that have become monomorphic or fully ungenotyped as a result of this first step are eliminated.
2) Next, **loci** that do not pass the given thresholds (QUAL>=30; F_MISSING<=0.5) are removed. At this step, all site level statistics that are provided in GATK VCF outputs can be used to filter out unwanted loci.
3) You can now decide to remove poorly covered **samples** if their proportion of ungenotyped loci is above a given threshold. 
4) Additional site level statistics are then computed: general information (SNP or INDEL, number of alleles, MAF…) and population genetics statistics (Fis, He) are added to the VCF ‘INFO’ field for each variant.
5) The statistics from step 4. (along with all other GATK stats) can be used in a second **loci** filtering step.
6) A MultiQC report is created, showing variants information/statistics after each filtering step.
7) Based on the variant statistics, histograms are created to estimate the quality of the variant calling after filtering, as well as a boxplot of the observed depth at each genotype (Locus x Sample) and a plot showing the detected variants positioned along the reference genome's contigs


## QUICK START

Needed files:  
- the full GeCKO/ folder, including the runGeCKO.sh launcher  
- your variant calling .vcf.gz files
- the cluster_config_VcfFiltering.yml (in case you work on a cluster) and config_VcfFiltering.yml files 

&nbsp;

To easily launch the workflow, use the runGeCKO.sh launcher. For example, to launch the workflow on the VCF example dataset on a Slurm job-scheduler, run the following command from the EXAMPLE directory:  
```../../../runGeCKO.sh --workflow VcfFiltering --workflow-path ../../../../GeCKO --config-file CONFIG/config_VcfFiltering.yml --cluster-config CONFIG/cluster_config_VcfFiltering_SLURM.yml --jobs 20 --job-scheduler SLURM```  

To launch it on your own data, if you cloned the repository in /home/user and placed your config_DataCleaning.yml and cluster_config_DataCleaning.yml files in a CONFIG folder:  
```WORKFLOW_PATH=/home/user/GeCKO```  
```${WORKFLOW_PATH}/runGeCKO.sh --workflow VcfFiltering --workflow-path ${WORKFLOW_PATH} --jobs 50 --job-scheduler SLURM```  

&nbsp;

![](https://github.com/GE2POP/GeCKO/blob/main/readme_img/VcfFiltering_4elements.png)


## How to use the VCF_FILTERING workflow
 
1) [Clone the GitHub repository](#1-clone-the-github-repository)  
2) [Prepare your input data](#2-prepare-your-input-data)  
3) [Prepare the CONFIG files](#3-prepare-the-config-files)  
4) [Launch the analysis](#4-launch-the-analysis)  
5) [Expected outputs](#5-expected-outputs)


### 1/ Clone the GitHub repository

Follow the procedure described [here](https://github.com/GE2POP/GeCKO/tree/main#readme) to have GeCKO ready to process your data.

### 2/ Prepare your input data

The input data is the .vcf file obtained after the variants calling step (variants_calling.vcf.gz) and its associated index file (variants_calling.vcf.gz.tbi). 

### 3/ Prepare the config files

The VCF_FILTERING workflow will need information about the dataset and the analysis parameters to perform its different steps.  
These information are provided through two files: *cluster_config_VcfFiltering.yml* and *config_VcfFiltering.yml*.  
If you name them exactly as written above and place them in a folder named 'CONFIG', the bash launching script will detect them automatically. Otherwise, you will have to pass them as arguments with --config and --cluster-config (see [below](#4-launch-the-analysis) for details).

#### *A/ The cluster_config_VcfFiltering.yml file:*
This file will be needed if you run the workflow on a computer cluster and want Snakemake to submit jobs. You will <ins>only need to modify two things: the partitions or queues names</ins> to match those of your cluster, and <ins>the memory to be requested for each submitted job</ins>. The first section of the file gives the default values for the job-scheduler's parameters that Snakemake should use for all its steps (or rules). The following sections correspond to specific Snakemake steps, with new parameters values to overwrite the defaults. If you want to assign a different partition/queue or memory requirement for a specific step that does not yet have its own section, you can create a new section for it:  

	specificStepName:
    	q or partition: {partition name for specificStep}
    	mem-per-cpu or h_vmem: {needed memory for each job submitted in specificStep}

You will find [the list of the steps names](#list-of-the-snakefile-rules) along with what they do and the tools they use at the end of this page.  
Our workflows support SGE and Slurm job-schedulers. <ins>You will find cluster-config files for both in the EXAMPLE/CONFIG folder</ins>.  


#### *B/ The config_VcfFiltering.yml file:*  
This file is used to pass all the information and tools parameters that will be used by the VCF_FILTERING workflow. The workflow expects it to contain a specific list of variables and their assigned values, organized in YAML format. Expected variables are:  

**INPUT FILES**  
- *VCF_FILE*&nbsp;&nbsp;&nbsp;The path to the variant calling file in zipped vcf format (.vcf.gz).

**VCF FILTERING PARAMETERS**  

- *BCFTOOLS_GENOTYPE_FILTERING_OPTIONS:*&nbsp;&nbsp;&nbsp; Genotype filtering parameters passed to bcftools filter (e.g. “FMT/DP >= 5 & FMT/GQ >= 15”). Any genotype not meeting these conditions will be assigned a missing value. Be careful to provide them between quotes.  
- *BCFTOOLS_LOCUS_FILTERING1_OPTIONS:*&nbsp;&nbsp;&nbsp; Locus filtering parameters passed to bcftools filter (e.g. "QUAL >= 30 & F_MISSING <= 0.5"). Any locus not meeting these conditions will be removed. Be careful to provide them between quotes.  
- *MAX_NA_PER_SAMPLE:*&nbsp;&nbsp;&nbsp; The maximum proportion of allowed missing data per sample. Samples above the threshold will be removed. Between 0 and 1. Set to "1" for no filtering.  
- *BCFTOOLS_LOCUS_FILTERING2_OPTIONS:*&nbsp;&nbsp;&nbsp; Locus filtering parameters passed to bcftools filter after removing poorly covered samples and computing additional site level statistics:  
	- SNP = Variant with at least 2 mononucleotide alleles (Type=Integer) [boolean]
	- INDEL = Variant with at least 2 alleles of different lengths (Type=Integer) [boolean]
	- nbGS = Number of genotyped samples (Type=Integer)
	- pNA = Frequency of missing genotypes for this site (Type=Float)
	- nbG = Number of genotypes (Type=Integer)
	- nbAll = Number of alleles (Type=Integer)
	- All = List of alleles (Type=String)
	- AllFreq = Frequency of all alleles (Type=Float)
	- MinHomo = Number of occurrences of the least frequent homozygous genotype (Type=Integer)
	- MAF = Minor allele frequency (Type=Float)
	- MAFb = Minor base frequency (Type=Float)
	- He = Nei expected Heterozygosity (Type=Float)
	- Fis = Inbreeding coefficient (Type=Float)

  e.g.: " INFO/SNP==1 & INFO/INDEL==0 & INFO/nbAll==2 & INFO/pNA<=0.5 & INFO/MinHomo>=2". Any locus not meeting these conditions will be removed.  
  Be careful to use quotes when passing this parameter.  

&nbsp;

<ins>An example of config_VcfFiltering.yml file can be found in the EXAMPLE/CONFIG folder</ins>.  

&nbsp;

**UNDERSTANDING HOW TO USE THE PREVIOUS VARIABLES TO FILTER VARIANTS WITH BCFTOOLS**  

The genotype and locus filtering steps are performed with the [bcftools’ filter function](https://samtools.github.io/bcftools/bcftools.html#filter). The options you will pass to bcftools are the conditions that must be met for the genotypes or loci to be retained (--include). For more information on the syntax of the expressions that can be passed to bcftools filter, see [this page](https://samtools.github.io/bcftools/bcftools.html#expressions).  
In a VCF file, metrics can relate either to a genotype (sample x locus) or to a locus. Some metrics (e.g. DP = reads depth) can relate to both (e.g. genotype level reads depth or site level reads depth across all samples). To differentiate between the genotype’s DP and the locus’ DP, we need to specify **FORMAT**/DP (or **FMT**/DP) for the genotype or **INFO**/DP for the locus.  

- **Genotype filtering**  
Genotypes are filtered with the following bcftools command:  
```bcftools filter -i ‘[BCFTOOLS_GENOTYPES_FILTERING_OPTIONS]’ -S .```  
For this step you will have to use ‘FMT/’ in front of the metrics you want to use for filtering.  
For example, if you wrote in your config file:  
*BCFTOOLS_GENOTYPE_FILTERING_OPTIONS: "FMT/DP >= 5 & FMT/GQ >= 15"*  
Then the command will be:  
```bcftools filter -i ‘FMT/DP >= 5 & FMT/GQ >= 15’ -S .```  
Which will set to “.” (meaning NA in a VCF file) any genotype that does not pass DP >= 5 **AND** GQ >= 15.

- **Locus filtering (first step)**  
Loci are filtered with the following command:  
```bcftools filter -i ‘[BCFTOOLS_LOCUS_FILTERING1_OPTIONS]’```  
At this step, you can filter loci using all the common GATK metrics given in the VCF INFO field (DP, InbreedingCoeff, etc.), the QUAL column, or the metrics/info that can be computed on the fly by bcftools (F_MISSING, MAF, TYPE, etc. See the bcftools documentation for an exhaustive list).  
For example, if you want to keep only biallelic SNP sites with quality above 30 and less than 50% missing values, here is what you could write in your config file:  
*BCFTOOLS_ LOCUS_FILTERING1_OPTIONS: "N_ALT=1 & TYPE=\\"snp\\" & QUAL >= 30 & F_MISSING <= 0.5"*

- **Locus filtering (second step)**  
After filtering out samples, loci can be filtered a second time with:  
```bcftools filter -i ‘[BCFTOOLS_LOCUS_FILTERING2_OPTIONS]’```  
At this step, you can use all the previously mentioned GATK and bcftools metrics, as well as some newly computed metrics listed above in the BCFTOOLS_LOCUS_FILTERING2_OPTIONS description. Some of them are redundant with the bcftools metrics but we thought it could be useful to have them written down in the final VCF file. Also, you may want to filter loci again after possibly removing some samples from your dataset to make sure that all loci still pass your desired thresholds.  
For example, if you want to make sure all variants still have less than 50% missing data and present at least two homozygous genotypes of each allele among the remaining samples, you could write in your config file:  
*BCFTOOLS_ LOCUS_FILTERING1_OPTIONS: "F_MISSING <= 0.5 & INFO/MinHomo >= 2"*


### 4/ Launch the analysis

**Environment**  
You can run this workflow on a computer or on a computer cluster. You will need Snakemake and Conda to be available.

**Launching**  
To launch the VCF_FILTERING workflow, use the launching script runGeCKO.sh with the option --workflow VcfFiltering:  
```WORKFLOW_PATH=/home/user/GeCKO```  
```${WORKFLOW_PATH}/runGeCKO.sh --workflow VcfFiltering --workflow-path ${WORKFLOW_PATH} --jobs 50 --job-scheduler SLURM``` 
 
For more help on how to use the launcher, see GeCKO's general [README](https://github.com/GE2POP/GeCKO/tree/main#quick-start), or run:  
```${WORKFLOW_PATH}/runGeCKO.sh --help --workflow-path ${WORKFLOW_PATH}```  

 
### 5/ Expected outputs  
This workflow will create a "VCF_FILTERING" directory in the "WORKFLOWS_OUTPUTS" directory. This directory is structured as follows and contains:  

<img src="https://github.com/GE2POP/GeCKO/blob/main/readme_img/OutputsTree_VcfFiltering.png" width="600"/>



<ins>Description of the main files:</ins>  

- *01__Genotype_Filtered.vcf*:&nbsp;&nbsp;&nbsp;the vcf file after filtering variants by **genotype**
- *02__Genotype_Locus1_Filtered.vcf*:&nbsp;&nbsp;&nbsp;the vcf file after the first variants filtering by **locus**  
- *03__Genotype_Locus1_Sample_Filtered.vcf*:&nbsp;&nbsp;&nbsp;the vcf file after filtering variants by **sample** 
- *samples_to_remove.list*:&nbsp;&nbsp;&nbsp;list of samples that were deleted in step 03 (filtering by sample)
- *03__Genotype_Locus1_Sample_Filtered__withExtraStats.vcf*:&nbsp;&nbsp;&nbsp;the intermediate vcf file corresponding to variants filtered after step 03 (genotype + locus1 + sample), with additional site level statistics (Fis, He, MAF,...) by variants.  
- *04__Genotype_Locus1_Sample_Locus2_Filtered.vcf*:&nbsp;&nbsp;&nbsp;the vcf file after filtering variants by **locus** 
- *workflow_info.txt*:&nbsp;&nbsp;&nbsp;File that contains the date and time of the workflow launch, the link to the Github repository and the commit ID

**REPORTS directory** 
- *00_variants_raw_vcf.stats*:&nbsp;&nbsp;&nbsp;bcftools statistics of the unfiltered vcf file 
- *01__Genotype_Filtered.stats*:&nbsp;&nbsp;&nbsp;bcftools statistics after filtering by **genotype**
- *02__Genotype_Locus1_Filtered.stats*:&nbsp;&nbsp;&nbsp;bcftools statistics after the first filtering by **locus**
- *03__Genotype_Locus1_Sample_Filtered.stats*:&nbsp;&nbsp;&nbsp;bcftools statistics after filtering by **sample**
- *04__Genotype_Locus1_Sample_Locus2_Filtered.stats*:&nbsp;&nbsp;&nbsp;bcftools statistics after the second filtering by **locus**
- *multiQC_VcfFiltering_report.html*:&nbsp;&nbsp;&nbsp;graphic representation of the variants informations/statistics after each filtering step
- *variants_stats_VF.tsv*:&nbsp;&nbsp;&nbsp; file that summarizes the statistics per locus present in the vcf file after filtering
- *variants_stats_histograms_VF.pdf*:&nbsp;&nbsp;&nbsp; file with histograms based on locus statistics after filtering
- *genotypes_DP_boxplot_VF.pdf*:&nbsp;&nbsp;&nbsp; file with a boxplot of the observed depth at each genotype, along with the percentage of missing values in the vcf file
- *variants_along_genome_VF.pdf*:&nbsp;&nbsp;&nbsp;&nbsp; file with a plot showing the detected variants positioned along the reference genome's contigs

## Tools
This workflow uses the following tools: 
- [bcftools v1.15](https://samtools.github.io/bcftools/bcftools.html)
- [multiqc v1.11](https://github.com/ewels/MultiQC/releases)
- [egglib v3.1.0](https://www.egglib.org/)
- [seaborn v0.12.2](https://seaborn.pydata.org/)
- [matplotlib v3.2.1](https://matplotlib.org/)
- [pandas v1.3.5](https://pandas.pydata.org/)
- [numpy v1.23.1](https://numpy.org/)

These tools are loaded in a CONDA environment from the conda-forge and bioconda channels.

##  List of the snakefile rules
Name, description and tools used for each of the snakemake workflow rules:

| **Rule name**                    | **Description**                                                                           | **Tools**       |
|:--------------------------------:|:-----------------------------------------------------------------------------------------:|:---------------:|
| Filter_Genotypes                 | Filtering variants by genotype (locus X sample)                                           | bcftools filter |
| Filter_Loci_1                    | Filtering variants by locus - first step                                                  | bcftools filter |
| Filter_Samples                   | Filtering variants by sample (Removing samples that have too many loci with missing data) | bcftools view   |
| Calculate_LocusExtraStats        | Calculating additional site level statistics (e.g. FIS, He, MAF)                          | Egglib          |
| Filter_Loci_2                    | Filtering variants based on additional site level statistics (e.g. FIS, He, MAF)          | bcftools filter |
| Build_StatsReport                | Building statistics reports for unfiltered variants and each filtering step               | bcftools stats  |
| Build_Report                     | Running MultiQC on unfiltered variants and each filtering step                            | MultiQC         |
| Summarize_FinalVCFVariables      | Recovery and summarize locus statistics                                                   | bcftools query  |
| Plot_FinalVCFVariablesHistograms | Creating histograms based on locus statistics                                             | seaborn, pyplot |
| Plot_FinalVCFDPBoxplot           | Creating a boxplot of the depth at each genotype                                          | pyplot          |
| Plot_FinalVCFVariantsAlongGenome | Creating a plot of the detected variants along the genome                                 | pyplot          |



![](https://github.com/GE2POP/GeCKO/blob/main/readme_img/VcfFiltering_Workflow.jpg?raw=true)


