# VARIANTS CALLING

This VARIANTS_CALLING workflow generates a vcf file from bam files obtained after mapping the cleaned sequences. This workflow use GATK for call the variants.


### The VARIANTS_CALLING workflow's steps
1) An index of the provided reference is created if it does not exist yet.
2) A dictionary of the provided reference is created if it does not exist yet.
3) The list of chromosomes or contigs in the reference is created for the GenomicsDBImport step 
4) Variants calling by sample > with GATK - HaplotypeCaller
5) Create data base from Variants calling by sample and the list of chromosomes or contigs in the reference > with GATK - GenomicsDBImport
6) Variants calling for all samples ( population) from GENOMICSDBImport to vcf file >  with GATK - GenotypeGVCFs

![](https://github.com/BioInfo-GE2POP-BLE/CAPTURE_PIPELINES_SNAKEMAKE/blob/main/readme_img/VariantsCalling_Workflow.jpg?raw=true)


## QUICK START

To easily launch the workflow, use our runSnakemakeWorkflow.sh launcher:  
```./runSnakemakeWorkflow.sh --workflow VariantsCalling --workflow-path PATH/TO/CAPTURE_SNAKEMAKE_WORKFLOWS```  

Needed files:  
- the full CAPTURE_SNAKEMAKE_WORKFLOWS/ folder  
- the runSnakemakeWorkflow.sh launcher  
- your mapped .bam files
- a reference file in fasta format to call variants
- the cluster_config_VariantsCalling.yml (in case you work on a cluster) and config_VariantsCalling.yml files in a CONFIG folder  

&nbsp;

For example, if you need to launch the workflow on our BAMS example dataset on a Slurm job-scheduler, run the following command from the EXAMPLE directory:  
```./runSnakemakeWorkflow.sh --workflow VariantsCalling --workflow-path /storage/replicated/cirad_users/ardissonm/CAPTURE_SNAKEMAKE_WORKFLOWS --config-file CONFIG/config_VariantsCalling.yml --cluster-config CONFIG/cluster_config_VariantsCalling.yml --jobs 20 --job-scheduler SLURM```  


&nbsp;

![](https://github.com/BioInfo-GE2POP-BLE/CAPTURE\_PIPELINES\_SNAKEMAKE/blob/main/readme\_img/VariantsCalling\_4elements.png)


## How to use the VARIANTS_CALLING workflow
 
1) [Prepare your input data](#1-prepare-your-input-data)  
2) [Clone our GitHub repository](#2-clone-our-github-repository)  
3) [Prepare the CONFIG files](#3-prepare-the-config-files)  
4) [Launch the analysis](#4-launch-the-analysis)  
5) [Expected outputs](#5-expected-outputs)


### 1/ Prepare your input data

the input data are the mapped .bam files and the associated index files (.bam.bai). 
To do the variants calling, it is necessary to fill in the reference used for the mapping.


### 2/ Clone our GitHub repository

The CAPTURE_SNAKEMAKE_WORKFLOWS folder must be fully copied in a workspace/storage of your choice.  
For example, you can clone the repository with:  
```git clone git@github.com:BioInfo-GE2POP-BLE/CAPTURE_PIPELINES_SNAKEMAKE.git```   


### 3/ Prepare the config files

The VARIANTS_CALLING workflow will need information about the dataset and the analysis parameters to perform its different steps.  
These information are provided through two files: *cluster_config_VariantsCalling.yml* and *config_VariantsCalling.yml*.  
If you name them exactly as written above and place them in a folder named 'CONFIG', the bash launching script will detect them automatically. Otherwise, you will have to pass them as arguments with --config and --cluster-config (see [below](#4-launch-the-analysis) for details).

#### *cluster_config_VariantsCalling.yml file:*
This file will be needed if you run the workflow on a computer cluster and want Snakemake to submit jobs. You <ins>only need to modify the partitions or queues names</ins> to match those of your cluster. The first section of the file gives the default values for the job-scheduler's parameters that Snakemake should use for all its steps (or rules). The following sections correspond to specific Snakemake steps, with new parameters values to overwrite the defaults. If you want to assign a different partition/queue for a specific step that does not yet have its own section, you can create a new section for it, preceded by a comma:  

	specificStepName:
    	q or partition: {partitionNameForSpecificStep}

You will find [the list of the steps names](#list-of-the-snakefile-rules) along with what they do and the tools they use at the end of this page.  
Our workflows support SGE and Slurm job-schedulers. <ins>You will find cluster-config files for both in the EXAMPLE/CONFIG folder</ins>.  


#### *config_VariantsCalling.yml file:*  
This file is used to pass all the information and tools parameters that will be used by the READS_MAPPING workflow. The workflow expects it to contain a specific list of variables and their assigned values, organized in YAML format. Expected variables are:  

**GENERAL VARIABLES**  

**INPUT FILES**  
- *BAM_DIR:*&nbsp;&nbsp;&nbsp;The path to the directory containing the mapped bam files and the index file in .bam.bai format.
- *REFERENCE:*&nbsp;&nbsp;&nbsp;The path to the reference file in fasta format (must end with .fa, .fas or .fasta) used for the mapping.   

**VARIANTS CALLING PARAMETERS**  
For each of the three GATK steps, two options fields are available: options related to the use of java (JAVA_OPTIONS) and step-specific options (EXTRA_OPTIONS) , if not leave blank: ""

- *GATK_HAPLOTYPE_CALLER_JAVA_OPTIONS:*&nbsp;&nbsp;&nbsp; options for using java at the step HaplotypeCaller by GATK (eg: "-Xmx4g").
- *GATK_HAPLOTYPE_CALLER_EXTRA_OPTIONS:*&nbsp;&nbsp;&nbsp; options for the step GATK - Haplotypecaller
- *GATK_GENOMICS_DB_IMPORT_JAVA_OPTIONS:*&nbsp;&nbsp;&nbsp; options for using java at the step GenomicsDBImport by GATK (eg: "-Xmx30g").
- *GATK_GENOMICS_DB_IMPORT_EXTRA_OPTIONS:*&nbsp;&nbsp;&nbsp; options for the step GATK - GenomicsDBImport (eg: "--merge-contigs-into-num-partitions 20 --batch-size 50 --reader-threads 20").
- *GATK_GENOTYPE_GVCFS_JAVA_OPTIONS:*&nbsp;&nbsp;&nbsp;options for using java at the step GenotypeGVCFs by GATK (eg: "-Xmx30g").
- *GATK_GENOTYPE_GVCFS_EXTRA_OPTIONS:*&nbsp;&nbsp;&nbsp; options for the step GATK - GenotypeGVCFs (eg: "--include-non-variant-sites --heterozygosity 0.001)
- 
&nbsp;

<ins>An example of config_VariantsCalling.yml file can be found in the EXAMPLE/CONFIG folder</ins>.  

&nbsp;

### 4/ Launch the analysis

**Environment**  
You can run this workflow on a computer or on a computer cluster. You will need Snakemake and Conda to be available.

**Launching**  
To launch the VARIANTS_CALLING workflow, you can use our launching script runSnakemakeWorkflow.sh with the option --workflow VariantsCalling:  
```./runSnakemakeWorkflow.sh --workflow VariantsCalling --workflow-path PATH/TO/CAPTURE_SNAKEMAKE_WORKFLOWS```  

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
the workflow create a directory "VARIANTS_CALLING" in the "WORKFLOWS_OUTPUTS" directory. 
In this directory there are 3 folders corresponding to the output of the 3 steps of variants calling with GATK:
HAPLOTYPE_CALLER : two files by sample, the vcf.gz file ( sample.g.vcf.gz) and the index file associated (sample.g.vcf.gz.tbi). there is also a list of vcf files contained in this folder (vcf.list.txt).
GENOMICS_DB_IMPORT: several directories containing the data base and associated files (.json, .vcf and . tdb)
GENOTYPE_GVCFS: the variants_calling.vcf.gz file and the index associated (variants_calling.vcf.gz.tbi)


## Tools
This workflow uses the following tools: 
- [gatk 4.2.5.0](https://github.com/broadinstitute/gatk/)
- [samtools 1.15](https://github.com/samtools/samtools/)

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

