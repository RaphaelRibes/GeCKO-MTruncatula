# CAPTURE_PIPELINES_SNAKEMAKE : DATA CLEANING

## INTRODUCTION

The DATA_CLEANING workflow allows several steps to be carried out successively in order to obtain cleaned sequences.

This workflow can be used with: 
	- paired_end sequences (noted PE): one way sequencing and in output one fastq.gz file
	- single_end sequences (noted SE): two way sequencing (Read1 et read2) and in output two files: R1.fastq.gz and R2.fastq.gz 

moreover, input sequences can be :
	- multiplexed: sequences of several samples in one (SE) or two (PE) fastq.gz files > all stages of the pipeline are completed 
	- demultiplexed: sequences in one (SE) or two (PE) fastq.gz files per sample > the pipeline automatically starts at step 4

_**EXAMPLE ??? **_

The different stages of the pipeline:
1/  Quality analysis of raw sequences/reads from sequencing ( FASTQC )
2/	Counting the number of reads in one (SE) or two (PE) fastq.gz files
3/	Demultiplex in one (SE) or two (PE) fastq.gz files into individual fastq files based on the barcode or tag specifci to each samples ( CUTADAPT )
4/	Counting the number of reads in one (SE) or two (PE) fastq.gz files by samples
5/	Trimming in one (SE) or two (PE) fastq.gz files by samples to remove adapters sequences, low quality sequences et short sequences ( CUTADAPT )
6/	Counting the number of reads in one (SE) or two (PE) fastq.gz files by samples
7/  Quality analysis of reads in fastq.gz file(s) par samples, after demultiplexing and trimming ( FASTQC ) 
8/ 	Creation of two reports ( MultiQC ) allowing to visualise the impact of the trimming step on the quality of the sequences (observation by samples), and the impact of the whole workflow (merged obervations). 


![DataCleaning_Workflow](C:\Users\ardisson\Documents\Analyses donnees NGS\GitHub\CAPTURE_PIPELINES_SNAKEMAKE\readme_img\DataCleaning_Workflow.jpg)



## TOOLS
This workflow uses 3 tools: 
  - [cutadapt version 3.5 ](https://cutadapt.readthedocs.io/en/v3.5/)
  - [fastqc version 11.9]( https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) 
  - [multiqc version 1.11](https://github.com/ewels/MultiQC/releases)
 
 These tools are loaded in a CONDA environment via (conda-forge and bioconda).
 
## SUMMARY
 
1/ Import DATA
2/ Clone WORKFLOW directory
3/ Prepare the CONFIG files
4/ Launch the analysis

##  PROCEDURE FOR USING WORKFLOW

![DataCleaning_4elements](C:\Users\ardisson\Documents\Analyses donnees NGS\GitHub\CAPTURE_PIPELINES_SNAKEMAKE\readme_img\DataCleaning_4elements.png)

### 1/ Import DATA

the input data are sequences from an Illumina sequencer (Miseq / Hiseq). 
	- paired_end sequences (PE): one way sequencing and in output one fastq.gz file
	- single_end sequences (SE): two way sequencing (Read1 et read2) and in output two files: R1.fastq.gz and R2.fastq.gz 

moreover, input sequences can be :
	- multiplexed: sequences of several samples in one (SE) or two (PE) fastq.gz files.
	- demultiplexed: sequences in one (SE) or two (PE) fastq.gz files per sample.

The sequences are stored in a directory whose path is specified in the config_file.txt


### 2/ clone WORKFLOW directory 

WORKFLOW directory to be fully cloned in a workspace/storage of your choice. 
To do this, use the command: 
git clone git@github.com:BioInfo-GE2POP-BLE/CAPTURE_PIPELINES_SNAKEMAKE/tree/main/DATA_CLEANING/WORKFLOW.git

WORKFLOW directory contains 3 folders:

ENVS
  **  _Conda.tools.yml_** : file containing the list of tools needed for the conda environment. This file is used in the sanemake workflow in the "conda" parameter of the rules. 

PROFILES
    SGE
        **_config.yml_** : file containing the parameters for launching a job on an SGE cluster (qsub, partition, etc). this file is called by the **run??? script  >>> à compléter**
    SLURM
        **_config.yml_** : file containing the parameters for launching a job on an SLURM cluster (qsub, partition, etc). this file is called by the **run??? script  >>> à compléter**

SCRIPTS
    **_demultiplex_with_cutadapt_PE.sh_** : For Paired End sequencing.  Bash script allowing to demultiplex a pair of fastq files (R1 + R2) into individual paired fastq files according to the barcode/tag, located at the P5 and P7 ends.
   ** _demultiplex_with_cutadapt_SE.sh_** : For Single End sequencing.  Bash script allowing to demultiplex  fastq file into individual fastq files according to the barcode/tag, located at the P5 ends.
    **_trimming_with_cutadapt.sh_PE_** : For Paired End sequencing.  Bash script allowing to trimming a pair of fastq files (R1 + R2 ) to remove adapters sequences, low quality sequences et short sequences. Use of tool CUTADAPT.
    **_trimming_with_cutadapt.sh_SE_** : For Single End sequencing. Bash script allowing to trimming a fastq files  to remove adapters sequences, low quality sequences et short sequences. Use of tool CUTADAPT.
    **_fastq_read_count.sh_** : bash script allowing to count reads in all fastq.gz files of folder
	
    
### 2/ prepare the CONFIG files

The WORKFLOW DATA_CLEANING needs several pieces of information about the dataset and the analysis parameters. 
All the files containing this information are located in the CONFIG folder which must be placed next to the DATA to be analysed in workspace.

**_cluster_config_SGE.json_** : file allowing to configure on a SGE-type computing cluster : the type of partition, the name of the log files, etc. It is possible to define a default version and/or to adapt it to each of the snakmake workflow rules
_cluster_config_SLURM.json_ : file allowing to configure on a SLURM-type computing cluster : the type of queue, the name of the log files, etc. It is possible to define a default version and/or to adapt it to each of the snakmake workflow rules

**_config_file.yml_** : 
file containing the information and parameters that will be used by the snakemake workflow DATA_CLEANING

> General variables:
OUTPUTS_DIRNAME: Name of the directory that will contain all the workflow outputs (example: WOKFLOW_OUTPUTS )
> Input files:
_for PairedEnd config.file.txt_                                                                 _for SingleEnd config.file.txt_
FASTQ_R1: path to the folder containing R1.fastq.gz file                      FASTQ: path to the folder containing fastq.gz file
FASTQ_R2: path to the folder containing R2.fastq.gz file
DEMULT_DIR: in the case of an analysis of already demultiplexed data, path to the folder containing the data.
                         in the case of a demultiplexed data analysis, throw the empty field (note "").
BARCODE_FILE:  in case the data has to be demultiplexed using the barcodes specific to each sample, path to the barcode_file.txt. see description below.
                               in the case of an analysis of already demultiplexed data, throw the empty field (note "").                              
ADAPTER_FILE: path to the adapter_file.txt. see description below.
> Demultiplexing parameters 
DEMULT_THREADS: number of threads to be allocated on cluster for the demultiplexing step with CUTADAPT (corresponds to parameter ""--nodes" of CUTADAPT tool)
DEMULT_SUBSTITUTIONS: percentage of substitution by barcode (tag). Example: 1 substitution on a barcode of 8pb, note 0.15. (corresponds to parameter ""--substitutions" of                                                    CUTADAPT tool)
> Trimming parameters
TRIMMING_THREADS: number of threads to be allocated on cluster for the trimming step with CUTADAPT (corresponds to parameter ""--nodes" of CUTADAPT tool)
TRIMMING_QUAL: parameter can be used to trim low-quality ends from reads. exemple:  30 > replacement of nucleotides by N if the quality is lower than Q30 (1 chance out of 30                                 that the base is wrong) (corresponds to parameter ""--quality_cutoff" of CUTADAPT tool)
TRIMMING_MIN_LENGTH: parameter to indicate the minimum size of the sequences to be kept. (corresponds to parameter ""--minimum_length" of CUTADAPT tool)


**_barcode_file.txt_** : 
File containing information about the barcode specific to each sample, if demultiplexing required. 

For Paired End sequencing: 
Containing the list of samples names (column 1) , sequences of barcode for read 1 (P5) (column 2) and sequences of barcode for read 2 (P7) (column 3)
 example (no header, tab-separated):
	Tc2208a	AGCGCA	AGCGCA
	Tc2235a	CTCAGC	CTCAGC
	Tc2249a	TAGATC	TAGATC	

For Single End sequencing: 
Containing the list of samples names (column 1) , sequences of barcode (column 2)
example (no header, tab-separated):
	Tc2208a	AGCGCA
	Tc2235a	CTCAGC
	Tc2249a	TAGATC


**_adapter_file.txt_** : 

Illumina sequencing starts at the first base of the fragment to be sequenced or the barcode. 
The CUTADAPT tool removes adaptor sequences after the fragment to be sequenced which appears when the fragment size is smaller than the sequencing unit. 
adapter_file.txt containing, for each genotype, the sequence of adapters to be deleted when trimming with cutadapt. 
Attention, the sequences of the adapters to be filled in must respect the reading directions 

the sequence of illumina adapters and index i5 / i7 are available: [illumina-adapter-sequences](https://support-docs.illumina.com/SHARE/AdapterSeq/illumina-adapter-sequences.pdf)


Example with Paired End sequencing:

![fragment_structure](C:\Users\ardisson\Documents\Analyses donnees NGS\GitHub\CAPTURE_PIPELINES_SNAKEMAKE\readme_img\DataCleaning_FragmentStructure.jpg )

sequence of the adapter to be deleted on read1 (5'-3'):
[barcode7revcomp]AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC[index7]ATCTCGTATGCCGTCTTCTGCTTG
sequence of the adapter to be deleted on read2 (5'-3'):
[barcode5revcomp]AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT[index5revcomp]GTGTAGATCTCGGTGGTCGCCGTATCATT

_For Paired End sequencing_: 
adapter_file.txt containing the list of samples names (column 1), sequences of adapter in the direction of read 1 after the sequencing fragment (column 2) and sequences of adapter in the direction of read 2 after the sequencing fragment (column 3)
	example (no header, tab-separated):
        Tc2208a	TGCGCTAGATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCCGCATCTCGTATGCCGTCTTCTGCTTGA	TGCGCTAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTGGCCGTATCATTA
        Tc2235a	GCTGAGAGATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCCGCATCTCGTATGCCGTCTTCTGCTTGA	GCTGAGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTGGCCGTATCATTA
        Tc2249a	GATCTAAGATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCCGCATCTCGTATGCCGTCTTCTGCTTGA	GATCTAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTGGCCGTATCATTA

_For Single End sequencing_: 
adapter_file.txt containing the list of samples names (column 1), sequences of adapter in the direction of read after the sequencing fragment (column 2)
	example (no header, tab-separated):
        Tc2208a	TGCGCTAGATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCCGCATCTCGTATGCCGTCTTCTGCTTGA
        Tc2235a	GCTGAGAGATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCCGCATCTCGTATGCCGTCTTCTGCTTGA
        Tc2249a	GATCTAAGATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCCGCATCTCGTATGCCGTCTTCTGCTTGA


### 4/ Run Workflow_DataCleaning

.condarc >>> l'utilisation de CONDA pour la gestion des outils (FASTQC, MULTIQC et CUTADAPT) nécessite la création d'un environnement (envs_dirs) et d'un package (pkgs_dirs). Lors du lancement du workflow Snakemake, ils sont créer par défault dans l'espace Homedir qui est soumis en générale à des restrictions d'écritures.
				   Compléter et copier ce fichier à l'endroit où se créer par défault les packages CONDA, pour rediriger la création le l'environnement et des packages CONDA vers l'espace de travail/stockage choisit 
  
** En cours de modification:**
run_DataCleaning_workflow_Slurm.sh >>> bash script to run the snakemake workflow DATA_CLEANING on a SLURM-type computing cluster. /!\ A MODIFIER /!\
											copy this bash script into the workspace next to the data to be processed. 
run_DataCleaning_workflow_SGE.sh >>> bash script to run the snakemake workflow DATA_CLEANING on a SGE-type computing cluster. /!\ A MODIFIER /!\
											copy this bash script into the workspace next to the data to be processed. 




