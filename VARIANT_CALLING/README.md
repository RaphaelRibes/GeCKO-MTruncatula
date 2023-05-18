# VARIANT CALLING

This VARIANT_CALLING workflow generates a vcf file from bam files obtained after mapping your reads to a reference genome. This workflow uses GATK to call variants.


### The VARIANT_CALLING workflow's steps
1) An index of the provided reference is created if it does not exist yet
2) A dictionary of the provided reference is created if it does not exist yet
3) The list of chromosomes or contigs in the reference is created for the GenomicsDBImport step 
4) Variant calling by sample is performed with the GATK HaplotypeCaller function
5) A database from variants calling by sample is generated with the GATK GenomicsDBImport function, and a list of the reference's chromosomes or contigs is created
6) Variants calling for all samples (population) is performed with the GATK GenotypeGVCFs function, creating a single vcf file
7) [Optional] (if extracting bams in a sub reference): convert the positions of the variants in the variant file (vcf file) with the positions given in the genomic reference
8) Based on the variant statistics calculated by GATK, histograms are created to estimate the quality of the variant calling before filtration, as well as a boxplot of the observed depth at each genotype (Locus x Sample) and a plot showing the detected variants positioned along the reference genome's contigs


## QUICK START

To easily launch the workflow, use the runGeCKO.sh launcher:  
```bash runGeCKO.sh --workflow VariantCalling --workflow-path /home/user/GeCKO```  

Needed files:  
- the full GeCKO/ folder  
- the runGeCKO.sh launcher  
- your mapped .bam files
- the reference file in fasta format that was used to map your reads
- the cluster_config_VariantCalling.yml (in case you work on a cluster) and config_VariantCalling.yml files in a CONFIG folder  

&nbsp;

For example, if you need to launch the workflow on our BAMS example dataset on a Slurm job-scheduler, run the following command from the EXAMPLE directory:  
```bash ../../runGeCKO.sh --workflow VariantCalling --workflow-path ../../../GeCKO --config-file CONFIG/config_VariantCalling.yml --cluster-config CONFIG/cluster_config_VariantCalling_SLURM.yml --jobs 20 --job-scheduler SLURM```  


&nbsp;

![](https://github.com/GE2POP/GeCKO/blob/main/readme\_img/VariantCalling\_4elements.png)


## How to use the VARIANT_CALLING workflow
 
1) [Prepare your input data](#1-prepare-your-input-data)  
2) [Clone our GitHub repository](#2-clone-our-github-repository)  
3) [Prepare the CONFIG files](#3-prepare-the-config-files)  
4) [Launch the analysis](#4-launch-the-analysis)  
5) [Expected outputs](#5-expected-outputs)


### 1/ Prepare your input data

The expected input data are .bam files and their associated index files (.bam.bai), along with the reference that was used for the mapping.  


### 2/ Clone our GitHub repository

The GeCKO folder must be fully copied in a workspace/storage of your choice.  
For example, you can clone the repository with:  
```git clone git@github.com:GE2POP/GeCKO.git```   


### 3/ Prepare the config files

The VARIANT_CALLING workflow will need information about the dataset and the analysis parameters to perform its different steps.  
These information are provided through two files: *cluster_config_VariantCalling.yml* and *config_VariantCalling.yml*.  
If you name them exactly as written above and place them in a folder named 'CONFIG', the bash launching script will detect them automatically. Otherwise, you will have to pass them as arguments with --config and --cluster-config (see [below](#4-launch-the-analysis) for details).

#### *A/ The cluster_config_VariantCalling.yml file:*
This file will be needed if you run the workflow on a computer cluster and want Snakemake to submit jobs. You will <ins>only need to modify two things: the partitions or queues names</ins> to match those of your cluster, and <ins>the memory to be requested for each submitted job</ins>. The first section of the file gives the default values for the job-scheduler's parameters that Snakemake should use for all its steps (or rules). The following sections correspond to specific Snakemake steps, with new parameters values to overwrite the defaults. If you want to assign a different partition/queue or memory requirement for a specific step that does not yet have its own section, you can create a new section for it:  

	specificStepName:
    	q or partition: {partition name for specificStep}
    	mem-per-cpu or h_vmem: {needed memory for each job submitted in specificStep}

You will find [the list of the steps names](#list-of-the-snakefile-rules) along with what they do and the tools they use at the end of this page.  
Our workflows support SGE and Slurm job-schedulers. <ins>You will find cluster-config files for both in the EXAMPLE/CONFIG folder</ins>.  


#### *B/ The config_VariantCalling.yml file:*  
This file is used to pass all the information and tools parameters that will be used by the VARIANT_CALLING workflow. The workflow expects it to contain a specific list of variables and their assigned values, organized in YAML format. Expected variables are:  

**GENERAL VARIABLES**  
- *VARIANT_CALLING_SUBFOLDER:*&nbsp;&nbsp;&nbsp;If you want to separate results from different variants calling parameters (different reference, mapping options...), provide a name for an extra folder to create in the VARIANT_CALLING output folder. Otherwise leave blank ("").  

**INPUT FILES**  
- *BAMS_LIST:*&nbsp;&nbsp;&nbsp;The path to the file containing the list of paths to the mapped bam files and index files in .bam.bai format. If the bams to be used for the VARIANT_CALLING have been mapped to the full genomic reference, use the file: bams_list.txt. If the bams to be used for the VARIANT_CALLING have been extracted on the basis of a sub reference (subbams), use the file: subbams_list.txt. These lists are generated by the READ_MAPPING workflow and are stored in: WORKFLOWS_OUTPUTS/READ_MAPPING
- *REFERENCE:*&nbsp;&nbsp;&nbsp;The path to the reference file in fasta format (must end with .fa, .fas or .fasta) used for the mapping. 
- *GENOMIC_REFERENCE_CHR_SIZE:*&nbsp;&nbsp;&nbsp;If your input bams result from an extraction of the reads mapping to specific genomic zones (i.e. CREATE_SUB_BAMS was set to TRUE during the mapping step) and you want the variants positions in this workflow's output vcf file to be given in the whole genomic reference, then please provide here the path to the reference_chr_size.txt file containing your genomic reference chromosomes sizes. This file is automatically created by the READ_MAPPING workflow when CREATE_SUB_BAMS is set to TRUE, and stored in WORKFLOWS_OUTPUTS/READ_MAPPING. Otherwise leave blank ("").  

**VARIANT CALLING PARAMETERS**  
For each of the three GATK steps, two options fields are available: options related to the use of java (JAVA_OPTIONS) and step-specific options (EXTRA_OPTIONS) , if not leave blank: ""

- *GATK_HAPLOTYPE_CALLER_JAVA_OPTIONS:*&nbsp;&nbsp;&nbsp;Java options for the GATK HaplotypeCaller function (eg: "-Xmx4g"). Be careful to provide them between quotes.
- *GATK_HAPLOTYPE_CALLER_EXTRA_OPTIONS:*&nbsp;&nbsp;&nbsp;Any list of options you would like to pass to the 'GATK Haplotypecaller' command. Be careful to provide them between quotes.
- *GATK_HAPLOTYPE_CALLER_CPUS_PER_TASK:*&nbsp;&nbsp;&nbsp;The number of CPUs to allocate for each Haplotypecaller task. Set to 1 if you are not working on a computing cluster. Be careful to never use quotes around this number.
- *GATK_GENOMICS_DB_IMPORT_JAVA_OPTIONS:*&nbsp;&nbsp;&nbsp;Java options for the GATK GenomicsDBImport function (eg: "-Xmx30g"). Be careful to provide them between quotes.
- *GATK_GENOMICS_DB_IMPORT_EXTRA_OPTIONS:*&nbsp;&nbsp;&nbsp;Any list of options you would like to pass to the 'GATK GenomicsDBImport' command (eg: "--merge-contigs-into-num-partitions 20 --batch-size 50 --reader-threads 20"). Be careful to provide them between quotes.
- *GATK_GENOMICS_DB_IMPORT_CPUS_PER_TASK:*&nbsp;&nbsp;&nbsp;The number of CPUs to allocate for the GenomicsDBImport step. Set to 1 if you are not working on a computing cluster. Be careful to never use quotes around this number.
- *GATK_GENOTYPE_GVCFS_JAVA_OPTIONS:*&nbsp;&nbsp;&nbsp;Java options for the GATK GenotypeGVCFs function (eg: "-Xmx30g"). Be careful to provide them between quotes.
- *GATK_GENOTYPE_GVCFS_EXTRA_OPTIONS:*&nbsp;&nbsp;&nbsp;Any list of options you would like to pass to the 'GATK GenotypeGVCFs' command (eg: "--include-non-variant-sites --heterozygosity 0.001). Be careful to provide them between quotes.

&nbsp;

<ins>An example of config_VariantCalling.yml file can be found in the EXAMPLE/CONFIG folder</ins>.  

&nbsp;


### 4/ Launch the analysis

**Environment**  
You can run this workflow on a computer or on a computer cluster. You will need Snakemake and Conda to be available.

**Launching**  
To launch the VARIANT_CALLING workflow, you can use our launching script runGeCKO.sh with the option --workflow VariantCalling:  
```./runGeCKO.sh --workflow VariantCalling --workflow-path PATH/TO/GeCKO```  

For more help on how to use it, see our GitHub's general README file or run:  
```./runGeCKO.sh --help --workflow-path PATH/TO/GeCKO```  

### 5/ Expected outputs  

This workflow will create a "VARIANT_CALLING" directory in the "WORKFLOWS_OUTPUTS" directory. This directory is structured as follows and contains:  

<img src="https://github.com/GE2POP/GeCKO/blob/main/readme_img/OutputsTree_VariantCalling.png" width="600"/>



<ins>Description of the main files:</ins> 

- *workflow_info.txt*:&nbsp;&nbsp;&nbsp;File that contains the date and time of the workflow launch, the link to the Github repository and the commit ID

**HAPLOTYPE_CALLER directory**  
- Two files by sample, the vcf.gz file (sample.g.vcf.gz) and the associated index file (sample.g.vcf.gz.tbi). A list of the vcf files contained in this folder will also be here (vcf.list.txt).  

**GENOMICS_DB_IMPORT directory** 
- Several directories containing the GATK data base and associated files (.json, .vcf and . tdb)  

**GENOTYPE_GVCFS directory**
- If the VARIANT_CALLING was performed on the full genomic reference: this folder contains final variant_calling.vcf.gz file and its associated index (variant_calling.vcf.gz.tbi)  
- If the VARIANT_CALLING was performed on the basis of a sub reference (subbams) and the positions of the variants (vcf file) were converted to the genomic reference, this folder contains final variant_calling_converted.vcf.gz file and its associated index (variant_calling_converted.vcf.gz.csi)
- **REPORTS directory** contains:  
    - *variants_stats_VC.tsv*:&nbsp;&nbsp;&nbsp;&nbsp; file that summarizes the statistics per locus present in the vcf file before filtering
    - *variants_stats_histograms_VC.pdf*:&nbsp;&nbsp;&nbsp;&nbsp; file with histograms based on locus statistics before filtering
    - *genotypes_DP_boxplot_VC.pdf*:&nbsp;&nbsp;&nbsp;&nbsp; file with a boxplot of the observed depth at each genotype, along with the percentage of missing values in the vcf file
    - *variants_along_genome_VC.pdf*:&nbsp;&nbsp;&nbsp;&nbsp; file with a plot showing the detected variants positioned along the reference genome's contigs

## Tools
This workflow uses the following tools: 
- [gatk v4.2.5.0](https://github.com/broadinstitute/gatk/)
- [samtools v1.14](https://github.com/samtools/samtools/)
- [bcftools v1.15](https://samtools.github.io/bcftools/bcftools.html)
- [seaborn v0.12.2](https://seaborn.pydata.org/)
- [matplotlib v3.2.1](https://matplotlib.org/)
- [pandas v1.3.5](https://pandas.pydata.org/)
- [numpy v1.23.1](https://numpy.org/)

These tools are loaded in a CONDA environment from the conda-forge and bioconda channels.

##  List of the snakefile rules
Name, description and tools used for each of the snakemake workflow rules:

| **Rule name**                     | **Description**                                                                 | **Tools**                     |
|:---------------------------------:|:-------------------------------------------------------------------------------:|:-----------------------------:|
| Index_Reference                   | Creating the reference index if needed                                          | samtools faidx                |
| Dictionary_Reference              | Creating the reference dictionnary for gatk if needed                           | gatk CreateSequenceDictionary |
| ListIntervalsReference_Dictionary | Listing chromosomes/contigs in the dictionnary for gatk GenomicsDBImport        |                               |
| HaplotypeCaller                   | Calling variants by sample                                                      | gatk HaplotypeCaller          |
| List_Haplotype                    | Listing sample files (g.vcf.gz) from HaplotypeCaller for gatk GenomicsDBImport  |                               |
| GenomicsDBImport                  | Creating data base from variants calling by sample and the intervals list       | gatk GenomicsDBImport         |
| GenotypeGVCFs                     | Calling variants for all samples (population) from GenomicsDBImport to vcf file | gatk GenotypeGVCFs            |
| ConvertPositions                  | Convert the positions of the variants (in vcf file) on the genomic reference    |                               |
| Summarize_GVCFVariables           | Recovery and summarize GATK locus statistics                                    | bcftools query                |
| Plot_GVCFVariablesHistograms      | Creating histograms based on GATK locus statistics                              | seaborn, pyplot               |
| Plot_GVCFDPBoxplot                | Creating a boxplot of the depth at each genotype                                | pyplot                        |
| Plot_GVCFVariantsAlongGenome      | Creating a plot of the detected variants along the genome                       | pyplot                        |


![Image non trouvée : https://github.com/GE2POP/GeCKO/blob/main/readme_img/VariantCalling_Workflow.jpg?raw=true](https://github.com/GE2POP/GeCKO/blob/main/readme_img/VariantCalling_Workflow.jpg?raw=true)
