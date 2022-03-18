# READS MAPPING

This READS_MAPPING workflow generates bams files from demultiplexed cleaned sequences.

It can be used to process:  
- Single-end sequences (SE), sequenced from only one end of each DNA fragment.  
- Paired-end sequences (PE), sequenced from both ends of each DNA fragment.  

### The READS_MAPPING workflow's steps
1) An index of the provided reference is created if it does not exist yet.
2) Reads are mapped to the reference, the resulting bams are sorted, the duplicates are removed if needed, and the final bams are indexed.
3) Bams reads are counted.
4) [Optionnal] For each genomic region provided in a bed file, reads are counted in each sample. A heatmap representing this data is generated.
5) [Optionnal] The reads that mapped to these regions are extracted and sub-bams are created. A corresponding sub-reference is also produced.
6) Two MultiQC reports are created, showing the reads numbers and quality after mapping, both before and after extracting reads from regions of interest.



## QUICK START

To easily launch the workflow, use our runSnakemakeWorkflow.sh launcher:  
```./runSnakemakeWorkflow.sh --workflow ReadsMapping --workflow-path PATH/TO/CAPTURE_SNAKEMAKE_WORKFLOWS```  

Needed files:  
- the full CAPTURE_SNAKEMAKE_WORKFLOWS/ folder  
- the runSnakemakeWorkflow.sh launcher  
- your demultiplexed and trimmed fastq.gz files
- a reference file in fasta format to map your reads unto
- the cluster_config_ReadsMapping.yml (in case you work on a cluster) and config_ReadsMapping.yml files in a CONFIG folder  
- a bed file listing genomic regions of interest
- (only for Morgane ;), temporary): the region.chain file that I sent you, placed in the folder where you run the workflow (in the future this file will be created automatically from the bed file)

&nbsp;

For example, if you need to launch the workflow on our ... dataset on a Slurm job-scheduler, run the following command from the EXAMPLE/... directory:  
```./runSnakemakeWorkflow.sh --workflow ReadsMapping --workflow-path /home/jogirodolle/save/CAPTURE_PIPELINES_SNAKEMAKE --config-file CONFIG/config_DataCleaning.yml --cluster-config CONFIG/cluster_config_Slurm_DataCleaning.json --jobs 20 --job-scheduler SLURM```  


&nbsp;


## How to use the READS_MAPPING workflow
 
1) [Prepare your input data](#1-prepare-your-input-data)  
2) [Clone our GitHub repository](#2-clone-our-github-repository)  
3) [Prepare the CONFIG files](#3-prepare-the-config-files)  
4) [Launch the analysis](#4-launch-the-analysis)  
5) [Expected outputs](#5-expected-outputs)


### 1/ Prepare your input data

The input data must be sequences from an Illumina sequencer (Miseq / Hiseq).  

Input sequences can be:  
- single-end sequences (SE): you must provide fastq files named in the format \*name\*.fastq.gz  
- paired-end sequences (PE): you must provide pairs of fastq files named in the format \*name\*.R1.fastq.gz and \*name\*.R2.fastq.gz  


### 2/ Clone our GitHub repository

The CAPTURE_SNAKEMAKE_WORKFLOWS folder must be fully copied in a workspace/storage of your choice.  
For example, you can clone the our repository with:  
```git clone git@github.com:BioInfo-GE2POP-BLE/CAPTURE_PIPELINES_SNAKEMAKE.git```   


### 3/ Prepare the config files

The READS_MAPPING workflow will need information about the dataset and the analysis parameters to perform its different steps.  
These information are provided through two files: *cluster_config_ReadsMapping.json* and *config_ReadsMapping.yml*.  
If you name them exactly as written above and place them in a folder named 'CONFIG', the bash launching script will detect them automatically. Otherwise, you will have to pass them as arguments with --config and --cluster-config (see [below](#4-launch-the-analysis) for details).

#### *cluster_config_ReadsMapping.json file:*
This file will be needed if you run the workflow on a computer cluster and want Snakemake to submit jobs. You <ins>only need to modify the partitions or queues names</ins> to match those of your cluster. The first section of the file gives the default values for the job-scheduler's parameters that Snakemake should use for all its steps (or rules). The following sections correspond to specific Snakemake steps, with new parameters values to overwrite the defaults. If you want to assign a different partition/queue for a specific step that does not yet have its own section, you can create a new section for it, preceded by a comma:  

	"specificStepName" : {
	"q" or "partition"         : "{partitionNameForSpecificStep}"
	}  

Our workflows support SGE and Slurm job-schedulers. <ins>You will find cluster-config files for both in the EXAMPLE/CONFIG folder</ins>.  

&nbsp;

#### *config_ReadsMapping.yml file:*  
This file is used to pass all the information and tools parameters that will be used by the READS_MAPPING workflow. The workflow expects it to contain a specific list of variables and their assigned values, organized in YAML format. Expected variables are:  

