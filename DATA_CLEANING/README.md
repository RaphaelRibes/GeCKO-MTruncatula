# CAPTURE_PIPELINES_SNAKEMAKE : DATA CLEANING

## INTRODUCTION

The DATA_CLEANING workflow allows several steps to be carried out successively in order to obtain cleaned sequences.

This workflow can be used with: 
	- paired_end sequences (noted PE): one way sequencing and in output one fastq.gz file
	- single_end sequences (noted SE): two way sequencing (Read1 et read2) and in output two files: R1.fastq.gz and R2.fastq.gz 

moreover, input sequences can be :
	- multiplexed: sequences of several samples in one (SE) or two (PE) fastq.gz files > all stages of the pipeline are completed 
	- demultiplexed: sequences in one (SE) or two (PE) fastq.gz files per sample > the pipeline automatically starts at step 4


The different stages of the pipeline:
1/  Quality analysis of raw sequences/reads from sequencing ( FASTQC )
2/	Counting the number of reads in one (SE) or two (PE) fastq.gz files
3/	Demultiplex in one (SE) or two (PE) fastq.gz files into individual fastq files based on the barcode or tag specifci to each samples ( CUTADAPT )
4/	Counting the number of reads in one (SE) or two (PE) fastq.gz files by samples
5/	Trimming in one (SE) or two (PE) fastq.gz files by samples to remove adapters sequences, low quality sequences et short sequences ( CUTADAPT )
6/	Counting the number of reads in one (SE) or two (PE) fastq.gz files by samples
7/  Quality analysis of reads in fastq.gz file(s) par samples, after demultiplexing and trimming ( FASTQC ) 
8/ 	Creation of two reports ( MultiQC ) allowing to visualise the impact of the trimming step on the quality of the sequences (observation by samples), and the impact of the whole workflow (merged obervations). 

**Ajouter un shéma récap**


## TOOLS
This workflow uses 3 tools: 
  - cutadapt version 3.5 > https://cutadapt.readthedocs.io/en/v3.5/
  - fastqc version 11.9 > https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
  - multiqc version 1.11 > https://github.com/ewels/MultiQC/releases
 
 These tools are loaded in a CONDA environment via (conda-forge and bioconda).
 
## SUMMARY
 
1/ Import DATA
2/ Clone WORKFLOW directory
3/ Prepare the CONFIG files
4/ Launch the analysis

##  PROCEDURE FOR USING WORKFLOW

### 1/ clone WORKFLOW directory 

Directory to be fully cloned in a workspace/storage of your choice. 
To do this, use the command: 
https://github.com/BioInfo-GE2POP-BLE/CAPTURE_PIPELINES_SNAKEMAKE/tree/main/DATA_CLEANING/WORKFLOW.git** >>> à vérifier**

This directory contains 3 folders:

ENVS
    _Conda.tools.yml_ : file containing the list of tools needed for the conda environment. This file is used in the sanemake workflow in the "conda" parameter of the rules. 

PROFILES
    SGE
        _config.ym_ : file containing the parameters for launching a job on an SGE cluster (qsub, partition, etc). this file is called by the **run??? script  >>> à compléter**
    SLURM
        _config.yml_ : file containing the parameters for launching a job on an SLURM cluster (qsub, partition, etc). this file is called by the **run??? script  >>> à compléter**

SCRIPTS
    _demultiplex_with_cutadapt_PE.sh_ : For Paired End sequencing.  Bash script allowing to demultiplex a pair of fastq files (R1 + R2) into individual paired fastq files according to the barcode/tag, located at the P5 and P7 ends.
    _demultiplex_with_cutadapt_SE.sh_ : For Single End sequencing.  Bash script allowing to demultiplex  fastq file into individual fastq files according to the barcode/tag, located at the P5 ends.
    _trimming_with_cutadapt.sh_PE_ : For Paired End sequencing.  Bash script allowing to trimming a pair of fastq files (R1 + R2 ) to remove adapters sequences, low quality sequences et short sequences. Use of tool CUTADAPT.
    _trimming_with_cutadapt.sh_SE_ : For Single End sequencing. Bash script allowing to trimming a fastq files  to remove adapters sequences, low quality sequences et short sequences. Use of tool CUTADAPT.
    _fastq_read_count.sh_ : bash script allowing to count reads in all fastq.gz files of folder
	
    
### 2/ prepare the CONFIG files

The WORKFLOW DATA_CLEANING needs several pieces of information about the dataset and the analysis parameters. 
All the files containing this information are located in the CONFIG folder which must be placed next to the DATA to be analysed in workspace.

_cluster_config_SGE.json_ : file allowing to configure on a SGE-type computing cluster : the type of partition, the name of the log files, etc. It is possible to define a default version and/or to adapt it to each of the snakmake workflow rules
_cluster_config_SLURM.json_ : file allowing to configure on a SLURM-type computing cluster : the type of queue, the name of the log files, etc. It is possible to define a default version and/or to adapt it to each of the snakmake workflow rules

_config_file.yml_ : 
file containing the information and parameters that will be used by the snakemake workflow DATA_CLEANING
**Détails de tous les paramètres? **

_sample_file.txt_ : 
File containing information about the samples: names and barcode (if demultiplexing required). 

For Paired End sequencing: 
Containing the list of samples names (column 1) , sequences of barcode/tags for read 1 (P5) (column 2) and sequences of barcode/tags for read 2 (P7) (column 3)
 example (no header, tab-separated):
	Tc2208a	AGCGCA	AGCGCA
	Tc2235a	CTCAGC	CTCAGC
	Tc2249a	TAGATC	TAGATC	

For Single End sequencing: 
Containing the list of samples names (column 1) , sequences of barcode/tags (column 2)
example (no header, tab-separated):
	Tc2208a	AGCGCA
	Tc2235a	CTCAGC
	Tc2249a	TAGATC


_adapter_file.txt_ : 
file containing, for each genotype, the sequence of adapters to be deleted when trimming with cutadapt. 

**Ajouter une figure avec les sequences, Read1, read2, ...**


For Paired End sequencing: 
Containing the list of samples names (column 1), sequences of adapter in the direction of read 1 after the sequencing fragment (column 2) and sequences of adapter in the direction of read 2 after the sequencing fragment (column 3)
	example (no header, tab-separated):
        Tc2208a	TGCGCTAGATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCCGCATCTCGTATGCCGTCTTCTGCTTGA	TGCGCTAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTGGCCGTATCATTA
        Tc2235a	GCTGAGAGATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCCGCATCTCGTATGCCGTCTTCTGCTTGA	GCTGAGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTGGCCGTATCATTA
        Tc2249a	GATCTAAGATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCCGCATCTCGTATGCCGTCTTCTGCTTGA	GATCTAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTGGCCGTATCATTA

For Single End sequencing: 
Containing the list of samples names (column 1), sequences of adapter in the direction of read after the sequencing fragment (column 2)
	example (no header, tab-separated):
        Tc2208a	TGCGCTAGATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCCGCATCTCGTATGCCGTCTTCTGCTTGA
        Tc2235a	GCTGAGAGATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCCGCATCTCGTATGCCGTCTTCTGCTTGA
        Tc2249a	GATCTAAGATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCCGCATCTCGTATGCCGTCTTCTGCTTGA


two **example** datasets available: 
- Dataset "DEV" : example of multiplexed data usable with the PairedEnd and SingleEnd version.
- dataset "MEL" : example of already demultiplexed data usable with the PairedEnd version.

Import the chosen folder into your workdir, and run the analysis with the script **run....sh <<< nom à modifier**


### 3/ import DATA

format entrée fastq.gz



### 4/ launch the analysis

.condarc >>> l'utilisation de CONDA pour la gestion des outils (FASTQC, MULTIQC et CUTADAPT) nécessite la création d'un environnement (envs_dirs) et d'un package (pkgs_dirs). Lors du lancement du workflow Snakemake, ils sont créer par défault dans l'espace Homedir qui est soumis en générale à des restrictions d'écritures.
				   Compléter et copier ce fichier à l'endroit où se créer par défault les packages CONDA, pour rediriger la création le l'environnement et des packages CONDA vers l'espace de travail/stockage choisit 
  
** En cours de modification:**
run_DataCleaning_workflow_Slurm.sh >>> bash script to run the snakemake workflow DATA_CLEANING on a SLURM-type computing cluster. /!\ A MODIFIER /!\
											copy this bash script into the workspace next to the data to be processed. 
run_DataCleaning_workflow_SGE.sh >>> bash script to run the snakemake workflow DATA_CLEANING on a SGE-type computing cluster. /!\ A MODIFIER /!\
											copy this bash script into the workspace next to the data to be processed. 




