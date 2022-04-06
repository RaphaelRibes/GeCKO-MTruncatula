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
- a reference file in fasta format to map your reads to
- the cluster_config_ReadsMapping.yml (in case you work on a cluster) and config_ReadsMapping.yml files in a CONFIG folder  
- a bed file listing genomic regions of interest
- (only for Morgane ;), temporary): the region.chain file that I sent you, placed in the folder where you run the workflow (in the future this file will be created automatically from the bed file)

&nbsp;

For example, if you need to launch the workflow on our ... dataset on a Slurm job-scheduler, run the following command from the EXAMPLE/... directory:  
```./runSnakemakeWorkflow.sh --workflow ReadsMapping --workflow-path /home/jogirodolle/save/CAPTURE_PIPELINES_SNAKEMAKE --config-file CONFIG/config_DataCleaning.yml --cluster-config CONFIG/cluster_config_Slurm_DataCleaning.json --jobs 20 --job-scheduler SLURM```  


&nbsp;

![](https://github.com/BioInfo-GE2POP-BLE/CAPTURE\_PIPELINES\_SNAKEMAKE/blob/main/readme\_img/ReadsMapping\_4elements.png)


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


#### *config_ReadsMapping.yml file:*  
This file is used to pass all the information and tools parameters that will be used by the READS_MAPPING workflow. The workflow expects it to contain a specific list of variables and their assigned values, organized in YAML format. Expected variables are:  

**GENERAL VARIABLES**  
- *PAIRED_END:*&nbsp;&nbsp;&nbsp;Whether your data is paired-end or single-end. [TRUE or FALSE]  
- *CREATE_SUB_BAMS:*&nbsp;&nbsp;&nbsp;Whether to extract reads from regions of interest (listed in bed file) and to create the corresponding sub-bams. Cannot be set to TRUE if the BED variable is left blank. [TRUE or FALSE]  
- *MAPPING_SUBFOLDER:*&nbsp;&nbsp;&nbsp;If you want to separate results from different mapping parameters (different reference, mapping options...), provide a name for an extra folder to create in the READS_MAPPING output folder. Otherwise leave blank ("").  

**INPUT FILES**  
- *TRIM_DIR:*&nbsp;&nbsp;&nbsp;The path to the directory containing the trimmed fastq files to be mapped. If left blank, the workflow will assume the fastq files are in WORKFLOWS_OUTPUTS/DATA_CLEANING/DEMULT_TRIM, which is the path to our DATA_CLEANING workflow output files.  
- *REFERENCE:*&nbsp;&nbsp;&nbsp;The path to the reference file in fasta format (must end with .fa, .fas or .fasta).  
- *BED:*&nbsp;&nbsp;&nbsp;The path to the bed file listing regions of interest to count reads in. Optionnal: can be left blank ("").  

**MAPPING PARAMETERS**  
- *MAPPER:*&nbsp;&nbsp;&nbsp;The name of the mapper you want to use. Currently implemented options are 'bwa-mem2_mem', 'bwa_mem', 'bowtie2' and 'minimap2'.  
- *REMOVE_DUP:*&nbsp;&nbsp;&nbsp;Whether or not to remove duplicates after mapping. [TRUE or FALSE]  
- *SEQUENCING_TECHNOLOGY:*&nbsp;&nbsp;&nbsp;The name of the sequencing technology (eg: "ILLUMINA"), which will appear in the reads names after mapping: 'PL:{SEQUENCING_TECHNOLOGY}')  
- *EXTRA_MAPPER_PARAMS:*&nbsp;&nbsp;&nbsp;Any list of options you would like to pass to the mapper command. Be careful to provide them between quotes.   
- *MARKDUP_PARAMS:*&nbsp;&nbsp;&nbsp;Any list of options you would like to pass to the 'picard MarkDuplicates' command, if REMOVE_DUP is set to TRUE. Be careful to provide them between quotes.  
- *INDEX_PARAMS:*&nbsp;&nbsp;&nbsp;Any list of options you would like to pass to the 'samtools index' command. Be careful to provide them between quotes. For example, you may need to pass the "-c" option if you need to map your reads to a very big reference file.  

&nbsp;

<ins>An example of config_ReadsMapping.yml file can be found in the EXAMPLE/CONFIG folder</ins>.  

&nbsp;

### 4/ Launch the analysis

**Environment**  
You can run this workflow on a computer or on a computer cluster. You will need Snakemake and Conda to be available.

**Launching**  
To launch the READS_MAPPING workflow, you can use our launching script runSnakemakeWorkflow.sh with the option --workflow ReadsMapping:  
```./runSnakemakeWorkflow.sh --workflow ReadsMapping --workflow-path PATH/TO/CAPTURE_SNAKEMAKE_WORKFLOWS```  

For more help on how to use it, see our GitHub's general README file or run:  
```./runSnakemakeWorkflow.sh --help --workflow-path PATH/TO/CAPTURE_SNAKEMAKE_WORKFLOWS```  

**Notes on Conda**  
The workflow will download and make available the [tools it needs](#tools) through Conda, which means you do not need to have them installed in your working environment behorehand.  
When called for the first time, the READS_MAPPING Snakemake workflow will download the tools' packages in a pkgs_dirs folder, and install them in a conda environment that will be stored in a .snakemake/conda folder, in the directory you called the workflow from. Every time you call the workflow from a new directory, the Conda environment will be generated again.  

The pkgs_dirs folder however is common to your whole system or cluster personnal environment. Conda's default behaviour is to create it in your home directory, in a .conda folder. If your home space is limited or if you do not have the right to write there from your cluster's nodes, you will need to tell Conda to store its packages somewhere else, thanks to a .condarc file. Place it in your home folder and specify the directory path where you want Conda to store the packages, following this example:  
```
envs_dirs:  
    - /home/username/path/to/appropriate/folder/env  
pkgs_dirs:  
    - /home/username/path/to/appropriate/folder/pkgs  
```




### 5/ Expected outputs  
...work in progress...


## Tools
This workflow uses the following tools: 
- [bwa-mem2 v2.2.1](https://github.com/bwa-mem2/bwa-mem2)
- [bwa v0.7.17](https://github.com/lh3/bwa)
- [bowtie2 v2.4.5](https://github.com/BenLangmead/bowtie2)
- [minimap2 v2.24](https://github.com/lh3/minimap2)
- [samtools v1.14](http://www.htslib.org/)
- [picard v2.26.10](https://broadinstitute.github.io/picard/)
- [seaborn v0.11.2](https://seaborn.pydata.org/)
- [matplotlib v3.5.1](https://matplotlib.org/)
- [multiqc v1.11](https://github.com/ewels/MultiQC/releases)
- [crossmap v0.6.3](http://crossmap.sourceforge.net/)

 
These tools are loaded in a CONDA environment from the conda-forge and bioconda channels.

