# READS MAPPING

This READS_MAPPING workflow generates bams files from demultiplexed cleaned sequences.

It can be used to process:  
- Single-end sequences (SE), sequenced from only one end of each DNA fragment.  
- Paired-end sequences (PE), sequenced from both ends of each DNA fragment.  

### The READS_MAPPING workflow's steps
1) An index of the provided reference is created if it does not exist yet.
2) Reads are mapped to the reference, the resulting bams are sorted, the duplicates are removed if needed, and the final bams are indexed.
3) Bams reads are counted.
4) [Optional] For each genomic region provided in a bed file, reads are counted in each sample.
5) [Optional] The reads that mapped to these regions are extracted and sub-bams are created. A corresponding sub-reference is also produced.
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

&nbsp;

For example, if you need to launch the workflow on our example PAIRED_END DEMULT_TRIM dataset on a Slurm job-scheduler, run the following command from the EXAMPLE/PAIRED_END directory:  
```../../../runSnakemakeWorkflow.sh --workflow ReadsMapping --workflow-path ../../../../CAPTURE_SNAKEMAKE_WORKFLOWS --config-file CONFIG/config_ReadsMapping.yml --cluster-config CONFIG/cluster_config_ReadsMapping_SLURM.yml --jobs 20 --job-scheduler SLURM```  


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
For example, you can clone the repository with:  
```git clone git@github.com:BioInfo-GE2POP-BLE/CAPTURE_SNAKEMAKE_WORKFLOWS.git```   


### 3/ Prepare the config files

The READS_MAPPING workflow will need information about the dataset and the analysis parameters to perform its different steps.  
These information are provided through two files: *cluster_config_ReadsMapping.yml* and *config_ReadsMapping.yml*.  
If you name them exactly as written above and place them in a folder named 'CONFIG', the bash launching script will detect them automatically. Otherwise, you will have to pass them as arguments with --config and --cluster-config (see [below](#4-launch-the-analysis) for details).

#### *A/ The cluster_config_ReadsMapping.yml file:*
This file will be needed if you run the workflow on a computer cluster and want Snakemake to submit jobs. You will <ins>only need to modify two things: the partitions or queues names</ins> to match those of your cluster, and <ins>the memory to be requested for each submitted job</ins>. The first section of the file gives the default values for the job-scheduler's parameters that Snakemake should use for all its steps (or rules). The following sections correspond to specific Snakemake steps, with new parameters values to overwrite the defaults. If you want to assign a different partition/queue or memory requirement for a specific step that does not yet have its own section, you can create a new section for it:  

	specificStepName:
    	q or partition: {partition name for specificStep}
    	mem-per-cpu or h_vmem: {needed memory for each job submitted in specificStep}

You will find [the list of the steps names](#list-of-the-snakefile-rules) along with what they do and the tools they use at the end of this page.  
Our workflows support SGE and Slurm job-schedulers. <ins>You will find cluster-config files for both in the EXAMPLE/.../CONFIG folders</ins>.  


#### *B/ The config_ReadsMapping.yml file:*  
This file is used to pass all the information and tools parameters that will be used by the READS_MAPPING workflow. The workflow expects it to contain a specific list of variables and their assigned values, organized in YAML format. Expected variables are:  

**GENERAL VARIABLES**  
- *PAIRED_END:*&nbsp;&nbsp;&nbsp;Whether your data is paired-end or single-end. [TRUE or FALSE]  
- *CREATE_SUB_BAMS:*&nbsp;&nbsp;&nbsp;Whether to extract reads from regions of interest (listed in bed file) and to create the corresponding sub-bams. Cannot be set to TRUE if the BED variable is left blank. [TRUE or FALSE]  
- *MAPPING_SUBFOLDER:*&nbsp;&nbsp;&nbsp;If you want to separate results from different mapping parameters (different reference, mapping options...), provide a name for an extra folder to create in the READS_MAPPING output folder. Otherwise leave blank ("").  

**INPUT FILES**  
- *TRIM_DIRS:*&nbsp;&nbsp;&nbsp;The path(s) to the directory or directories containing the trimmed fastq files to be mapped. If left blank, the workflow will assume the fastq files are in WORKFLOWS_OUTPUTS/DATA_CLEANING/DEMULT_TRIM, which is the path to our DATA_CLEANING workflow output files. To provide several directories, separate them with spaces, e.g.: "/home/user/trim_dir1 /home/user/trim_dir2". Be careful to provide them between quotes.   
- *REFERENCE:*&nbsp;&nbsp;&nbsp;The path to the reference file in fasta format (must end with .fa, .fas or .fasta).  

*If you set CREATE_SUB_BAMS to TRUE, you either have to provide a bed file:*
- *BED:*&nbsp;&nbsp;&nbsp;The path to the bed file listing regions of interest to count reads in. Optional: can be left blank ("").  

*Or to provide ALL of the three following parameters to automatically create a bed file containing the genomic regions with enough coverage in your dataset:*
- *BED_MIN_MEAN_COV:*&nbsp;&nbsp;&nbsp;The minimum mean depth per sample (number of reads per base) to keep a genomic region. Optional: can be left blank ("").
- *BED_MAX_DIST:*&nbsp;&nbsp;&nbsp;The maximum distance between two regions. If several regions are separated by a smaller distance than this, they will be merged into a single one. Optional: can be left blank ("").
- *BED_MIN_LENGTH:*&nbsp;&nbsp;&nbsp;The minimum length to keep a region after merging. Optional: can be left blank ("").

**MAPPING PARAMETERS**  
- *MAPPER:*&nbsp;&nbsp;&nbsp;The name of the mapper you want to use. Currently implemented options are 'bwa-mem2_mem', 'bwa_mem', 'bowtie2' and 'minimap2'.  
- *REMOVE_DUP:*&nbsp;&nbsp;&nbsp;Whether or not to remove duplicates after mapping. They will be marked either way. [TRUE or FALSE]  
- *SEQUENCING_TECHNOLOGY:*&nbsp;&nbsp;&nbsp;The name of the sequencing technology (eg: "ILLUMINA"), which will appear in the reads names after mapping: 'PL:{SEQUENCING_TECHNOLOGY}')  
- *EXTRA_MAPPER_OPTIONS:*&nbsp;&nbsp;&nbsp;Any list of options you would like to pass to the mapper command. Be careful to provide them between quotes. 
- *MAPPING_CPUS_PER_TASK:*&nbsp;&nbsp;&nbsp;The number of CPUs to allocate for each mapping task. Set to 1 if you are not working on a computing cluster. Be careful to never use quotes around this number.  
- *PICARD_MARKDUPLICATES_OPTIONS:*&nbsp;&nbsp;&nbsp;Any list of options you would like to pass to the 'picard MarkDuplicates' command. Be careful to provide them between quotes.  
- *PICARD_MARKDUPLICATES_JAVA_OPTIONS:*&nbsp;&nbsp;&nbsp;Java options to pass to the 'picard MarkDuplicates' command. Eg: "-Xmx4G".  

- *SAMTOOLS_INDEX_OPTIONS:*&nbsp;&nbsp;&nbsp;Any list of options you would like to pass to the 'samtools index' command. Be careful to provide them between quotes. For example, you may need to pass the "-c" option if you need to map your reads to a very big reference file.  

&nbsp;

<ins>Examples of config_ReadsMapping.yml files can be found in the EXAMPLE/.../CONFIG folders</ins>.  

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
When called for the first time, the READS_MAPPING Snakemake workflow will download the tools' packages in a pkgs_dirs folder, and install them in a conda environment that will be stored in a .snakemake/conda folder, in the directory you called the workflow from. Every time you call the workflow from a new directory, the Conda environment will be generated again. To avoid creating the environment multiple times, which can be both time and resource-consuming, you can provide a specific folder where you want Snakemake to store all of its conda environments with the --conda-env-path option of the runSnakemakeWorkflow.sh launcher.  

The pkgs_dirs folder is by default common to your whole system or cluster personnal environment. Conda's standard behaviour is to create it in your home directory, in a .conda folder. If your home space is limited or if you do not have the right to write there from your cluster's nodes, you will need to tell Conda to store its packages somewhere else, thanks to a .condarc file. Place it in your home folder and specify the directory path where you want Conda to store the packages, following this example:  
```
envs_dirs:  
    - /home/username/path/to/appropriate/folder/env  
pkgs_dirs:  
    - /home/username/path/to/appropriate/folder/pkgs  
```




### 5/ Expected outputs  

This workflow will create a "READS_MAPPING" directory in the "WORKFLOWS_OUTPUTS" directory. This directory is structured as follows and contains:  

<img src="https://github.com/BioInfo-GE2POP-BLE/CAPTURE_PIPELINES_SNAKEMAKE/blob/main/readme_img/OutputsTree_ReadsMapping.png" width="600"/>



<ins>Description of the main files:</ins>  

- *bams_list.txt*:&nbsp;&nbsp;&nbsp;File containing the list of paths to bams files
- *reference_chr_size.txt*:&nbsp;&nbsp;&nbsp;File containing the name and size of each chromosome in the genomic reference (just if "CREATE_SUB_BAMS: TRUE")
- *subbams_list.txt*:&nbsp;&nbsp;&nbsp;File containing the list of paths to subbams files (just if "CREATE_SUB_BAMS: TRUE")
- *workflow_info.txt*:&nbsp;&nbsp;&nbsp;File that contains the date and time of the workflow launch, the link to the Github repository and the commit ID

**BAMS directory**  
- *bams files*:&nbsp;&nbsp;&nbsp;One file per sample, reads mapped to the provided reference along with the associated index file (.csi).  

**BAMS/REPORTS directory**  
- *nb_reads_per_sample.tsv*:&nbsp;&nbsp;&nbsp;Number of reads per sample after mapping  
- *DUPLICATES directory*:&nbsp;&nbsp;&nbsp;MarkDuplicates reports (one per sample)  
- *STATS directory*:&nbsp;&nbsp;&nbsp;Samtools stats reports (one per sample)  
- *multiQC_ReadsMapping_Bams_Report.html (and the associated directory)*:&nbsp;&nbsp;&nbsp;Graphic report based on stats reports to visualize the number and percentages of mapped reads  

**EXTRACTED_BAMS/BAMS_ZONES directory** (if a bed file was provided and CREATE_SUB_BAMS was set to TRUE)  
- *bams files*:&nbsp;&nbsp;&nbsp;One file per sample, reads that mapped to the zones provided in the bed file, with coordinates given in the corresponding sub-reference  

**EXTRACTED_BAMS/BAMS_ZONES/REPORTS directory** (if a bed file was provided and CREATE_SUB_BAMS was set to TRUE)  
- *nb_reads_per_sample.tsv*:&nbsp;&nbsp;&nbsp;Number of reads per sample after mapping and extraction  
- *DUPLICATES directory*:&nbsp;&nbsp;&nbsp;MarkDuplicates reports (one per sample)  
- *STATS directory*:&nbsp;&nbsp;&nbsp;Samtools stats reports (one per sample)  
- *multiQC_ReadsMapping_SubBams_Report.html (and the associated directory)*:&nbsp;&nbsp;&nbsp;Graphic report based on stats reports to visualize the number and percentages of mapped reads in sub-bams\*  

**EXTRACTED_BAMS/REFERENCE_ZONES directory** (if a bed file was provided and CREATE_SUB_BAMS was set to TRUE)  
- *zones.bed*:&nbsp;&nbsp;&nbsp;Bed file containing the positions (start - end) of the genomic regions with enough coverage.
- *Reference_zones.fasta*:&nbsp;&nbsp;&nbsp;The sub-reference, corresponding to the extraction of the zones provided in the bed file.  

**ZONES_STATS directory** (if a bed file was provided) 
- *mean_depth_per_zone_per_sample.tsv*:&nbsp;&nbsp;&nbsp;For each zone and each sample, the mean depth per zone (number of reads per base) 

\*⚠ During the extraction step by CrossMap, paired reads are considered as a "proper pair" if both ends are mapped to different strands of the same contig and the distance between them is less than 500 pb. This information is stored in the bams' flags for each read (0x2).



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
- [pandas v1.3.5](https://pandas.pydata.org/)
- [multiqc v1.11](https://github.com/ewels/MultiQC/releases)
- [crossmap v0.6.3](http://crossmap.sourceforge.net/)

These tools are loaded in a CONDA environment from the conda-forge and bioconda channels.


##  List of the snakefile rules
Name, description and tools used for each of the snakemake workflow rules:

| **Rule name**               | **Description**                                                                                                                     | **Tools**                                                                                                                              |
|:---------------------------:|:-----------------------------------------------------------------------------------------------------------------------------------:|:--------------------------------------------------------------------------------------------------------------------------------------:|
| Index_Reference             | Creating the reference index if needed                                                                                              | bwa index // bwa-mem2 index // bowtie2-build // minimap2 -d                                                                            |
| Mapping_PairedEndFastqs     | Mapping the input fastq files onto the reference (paired end reads), sort bams, remove duplicates if needed, create bams index      | bwa mem // bwa-mem2 mem // bowtie2 // minimap2; samtools view; samtools fixmate; picard SortSam; picard MarkDuplicates; samtools index |
| Mapping_SingleEndFastqs     | Mapping the input fastq files onto the reference (single end reads)                                                                 | bwa mem // bwa-mem2 mem // bowtie2 // minimap2; samtools view; samtools fixmate;                                                       |
| MarkDuplicates_Bams         | sort bams and remove duplicates if needed                                                                                           | picard SortSam; picard MarkDuplicates                                                                                                  |
| Index_Bams                  | creating bams index                                                                                                                 | samtools index                                                                                                                         |
| Stats_Bams                  | Computing mapping stats                                                                                                             | samtools stats                                                                                                                         |
| Summarize_BamsReadsCount    | Summarizing reads count in bams from samtools stats output                                                                          |                                                                                                                                        |
| MultiQC_Bams                | Runing MultiQC on samtools stats output                                                                                             | MultiQC                                                                                                                                |
| Create_BamsList             | Writing list of bams files                                                                                                          |                                                                                                                                        |
| Create_BedFile              | Creating bed file containing the positions of zones                                                                                 | samtools merge; bedtools genomecov                                                                                                     |
| CountReadsZones_Bams        | Counting reads in each zone provided in the bed file for each sample                                                                | samtools bedcov                                                                                                                        |
| Create_SubReference         | Generating the reduced reference corresponding to the genomic zones provided in the bed file                                        | samtools faidx                                                                                                                         |
| Create_ChainFile            | Generating a file in chain format corresponding to the genomic zones provided in the bed file, necessary for the Extract_Reads rule | samtools faidx                                                                                                                         |
| Extract_Reads               | Extracting reads that mapped to the zones provided in the bed file, thus creating subbams                                           | crossmap                                                                                                                               |
| MarkDuplicates_Subbams      | sort subbams and remove duplicates if needed                                                                                        | picard SortSam; picard MarkDuplicates                                                                                                  |
| Index_Subbams               | creating subbams index                                                                                                              | samtools index                                                                                                                         |
| Stats_Subbams               | Computing mapping stats for subbams                                                                                                 | samtools stats                                                                                                                         |
| Summarize_SubbamsReadsCount | Summarizing reads count in subbams from samtools stats output                                                                       |                                                                                                                                        |
| MultiQC_Subbams             | Runing MultiQC on samtools stats output                                                                                             | MultiQC                                                                                                                                |
| Create_SubbamsList          | Writing list of subbams files                                                                                                       |                                                                                                                                        |
| Create_RefChrSizeFile       | Creating file with the name and size of each chromosome in the genomic reference                                                    | samtools faidx                                                                                                                         |


![](https://github.com/BioInfo-GE2POP-BLE/CAPTURE_PIPELINES_SNAKEMAKE/blob/main/readme_img/ReadsMapping_Workflow.jpg?raw=true)


