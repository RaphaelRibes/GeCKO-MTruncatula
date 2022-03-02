# CAPTURE_PIPELINES_SNAKEMAKE : DATA CLEANING

This DATA_CLEANING workflow generates demultiplexed cleaned sequences from raw sequenced data.

It can be used to process:  
- Single-end sequences (SE), sequenced from only one end of each DNA fragment.  
- Paired-end sequences (PE), sequenced from both ends of each DNA fragment.  

Besides, input sequences can be :
- Multiplexed (sequences of several samples packed in only one (SE) or a pair of (PE) fastq.gz files): all stages of the workflow will be completed
- Demultiplexed (sequences of each sample packed in separate fastq.gz files, one per sample if SE or a pair if PE): the workflow will skip steps 1 to 3



### The DATA_CLEANING workflow's steps  
1) Quality analysis of raw sequences/reads from sequencing (FASTQC)  
2) Counting the number of reads in one (SE) or two (PE) fastq.gz files  
3) Demultiplex in one (SE) or two (PE) fastq.gz files into individual fastq files based on the barcode or tag specifci to each samples (CUTADAPT)  
4) Counting the number of reads in one (SE) or two (PE) fastq.gz files by samples  
5) Trimming in one (SE) or two (PE) fastq.gz files by samples to remove adapters sequences, low quality sequences et short sequences (CUTADAPT)  
6) Counting the number of reads in one (SE) or two (PE) fastq.gz files by samples  
7) Quality analysis of reads in fastq.gz file(s) par samples, after demultiplexing and trimming (FASTQC)  
8) Creation of two reports (MultiQC) allowing to visualise the impact of the trimming step on the quality of the sequences (observation by samples), and the impact of the whole workflow (merged obervations).  

![DataCleaning_Workflow](https://github.com/BioInfo-GE2POP-BLE/CAPTURE_PIPELINES_SNAKEMAKE/blob/main/readme_img/DataCleaning_Workflow.jpg?raw=true)

## Tools
This workflow uses the following tools: 
- [Cutadapt v3.5 ](https://cutadapt.readthedocs.io/en/v3.5/)
- [FastQC v11.9](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) 
- [MultiQC v1.11](https://github.com/ewels/MultiQC/releases)
 
These tools are loaded in a CONDA environment from the conda-forge and bioconda channels.

&nbsp;
## How to use the DATA_CLEANING workflow
 
1) [Prepare your input data](#1-prepare-your-input-data)  
2) [Clone the WORKFLOW directory](#2-clone-the-workflow-directory)  
3) [Prepare the CONFIG files](#3-prepare-the-config-files)  
4) [Launch the analysis](#4-launch-the-analysis)  

&nbsp;

![DataCleaning_4elements](https://github.com/BioInfo-GE2POP-BLE/CAPTURE_PIPELINES_SNAKEMAKE/blob/main/readme_img/DataCleaning_4elements.png?raw=true)


### 1/ Prepare your input data

The input data must be sequences from an Illumina sequencer (Miseq / Hiseq).  

Input sequences can be multiplexed:  
- single-end sequences (SE): you must provide one fastq file named in the format \*name\*.fastq.gz  
- paired-end sequences (PE): you must provide two fastq files named in the format \*name\*_R1.fastq.gz and \*name\*_R2.fastq.gz  

or already demultiplexed:  
- single-end sequences (SE): you must provide one fastq file per sample named in the format \*sample\*.fastq.gz  
- paired-end sequences (PE): you must provide two fastq files per sample named in the format \*sample\*.R1.fastq.gz and \*sample\*.R2.fastq.gz  
The demultiplexed sequences must be placed together in a folder to be specified in the config_file.txt.  


### 2/ Clone the WORKFLOW directory

The DATA_CLEANING WORKFLOW folder must be fully copied or cloned in a workspace/storage of your choice.  
To do this, you can use the command:  
```git clone git@github.com:BioInfo-GE2POP-BLE/CAPTURE_PIPELINES_SNAKEMAKE/tree/main/DATA_CLEANING/WORKFLOW.git```  

    
### 3/ Prepare the config files

The DATA_CLEANING workflow will need information about the dataset and the analysis parameters to perform its different steps.  
These information are provided through two files: *cluster_config_DataCleaning.json* and *config_DataCleaning.yml*.  
If you name them exactly as written above and place them in a folder named 'CONFIG', the bash launching script will detect them automatically. Otherwise, you will have to pass them as arguments with --config and --cluster-config (see [below](#4-launch-the-analysis) for details).

#### *cluster_config_DataCleaning.json file:*
This file will be needed if you run the workflow on a computer cluster and want Snakemake to submit jobs. You <ins>only need to modify the partitions or queues names</ins> to match those of your cluster. The first section of the file gives the default values for the job-scheduler's parameters that Snakemake should use for all its steps (or rules). The following sections correspond to specific Snakemake steps, with new parameters values to overwrite the defaults. If you want to assign a different partition/queue for a specific step that does not yet have its own section, you can create a new section for it, preceded by a comma:  

	"specificStepName" : {
	"q" or "partition"         : "{partitionNameForSpecificStep}"
	}  

Our workflows support SGE and Slurm job-schedulers. <ins>You will find cluster-config files for both in the EXAMPLE/CONFIG folder</ins>.  

&nbsp;

#### *config_DataCleaning.yml file:*  
This file is used to pass all the information and tools parameters that will be used by the DATA_CLEANING workflow. The workflow expects it to contain a specific list of variables and their assigned values, organized in YAML format. Expected variables are:  

**GENERAL VARIABLES**  
*OUTPUTS_DIRNAME:*&nbsp;&nbsp;&nbsp;Name of the directory that the workflow will create, where to store its outputs (example: "WOKFLOW_OUTPUTS")  
*PAIRED_END:*&nbsp;&nbsp;&nbsp;Whether your data is paired-end or single-end [TRUE or FALSE]  

**INPUT FILES**  
*FASTQ:*&nbsp;&nbsp;&nbsp;Path to the raw data fastq.gz file (for single-end AND multiplexed data only, otherwise leave blank: "")  
*FASTQ_R1:*&nbsp;&nbsp;&nbsp;Path to the R1 raw data fastq.gz file, name format must be: name_R1.fastq.gz (for paired-end AND multiplexed data only, otherwise leave blank: "")  
*FASTQ_R2:*&nbsp;&nbsp;&nbsp;Path to the R2 raw data fastq.gz file, name format must be: name_R2.fastq.gz (for paired-end AND multiplexed data only, otherwise leave blank: "")  
*DEMULT_DIR:*&nbsp;&nbsp;&nbsp;Path to the folder containing the demultiplexed input files (for demultiplexed data only, otherwise leave blank: ""). Fastq files in the DEMULT_DIR folder must have the following name format: sampleX.R1.fastq.gz and sampleX.R2.fastq.gz for paired-end data; sampleX.fastq.gz for single-end data.  
*BARCODE_FILE:*&nbsp;&nbsp;&nbsp;Path to the barcode file (for multiplexed data only otherwise leave blank: ""). See description [below](#barcode-file) for more details.  
*ADAPTER_FILE:*&nbsp;&nbsp;&nbsp;Path to the adapter file. See description [below](#adapter-file) for more details.  


**DEMULTIPLEXING PARAMETERS** (mandatory if the raw fastq files are multiplexed, otherwise leave blank: "")  
*DEMULT_THREADS:*&nbsp;&nbsp;&nbsp;Number of threads to be allocated on your cluster for the demultiplexing step with Cutadapt (will be passed to the "--nodes" Cutadapt parameter)  
*DEMULT_SUBSTITUTIONS:*&nbsp;&nbsp;&nbsp;Fraction of authorized substitutions per barcode (tag). Example: to allow 1 substitution for an 8bp barcode, use '0.15'. (will be passed to the "--substitutions" Cutadapt parameter)  


**TRIMMIG PARAMETERS** (mandatory)  
*TRIMMING_THREADS:*&nbsp;&nbsp;&nbsp;Number of threads to be allocated on your cluster for the trimming step with Cutadapt (will be passed to the "--nodes" Cutadapt parameter)  
*TRIMMING_QUAL:*&nbsp;&nbsp;&nbsp;This parameter is used to trim low-quality ends from reads. Example:  If '30': nucleotides with quality score < Q30 (1 chance out of 1000 that the sequenced base is incorrect) will be replaced by N (will be passed to the "--quality_cutoff" Cutadapt parameter)
*TRIMMING_MIN_LENGTH:*&nbsp;&nbsp;&nbsp;parameter to indicate the minimum size of the sequences to be kept, after applying the TRIMMING_QUAL parameter (will be passed to the "--minimum_length" Cutadapt parameter)

&nbsp;

*IN A NUTSHELL*  
**Paired-end vs single-end**  
Expected variables differ for paired-end and single-end data: if PAIRED_END is set to TRUE, both FASTQ_R1 and FASTQ_R2 are expected. If PAIRED_END is set to FALSE, only FASTQ is expected.  

**Demultiplex or multiplexed data**  
Some variables must be left blank (*variable: ""*), depending on whether the raw data is multiplexed or demultiplexed. In case of multiplexed data, DEMULT_DIR must be left blank. In case of demultiplexed data, FASTQ, FASTQ_R1, FASTQ_R2, BARCODE_FILE, DEMULT_THREADS and DEMULT_SUBSTITUTIONS must be left blank.

<ins>Examples of config_DataCleaning.yml files adapted to each case can be found in the EXAMPLE/CONFIG folder</ins>.  

&nbsp;

#### *Barcode file:*  
This file must provide the barcode sequences specific to each sample (only required if the data is multiplexed).  

- For paired-end sequencing:  
Three-column file, specifying samples names (column 1), sequences of barcodes for read 1 (P5) (column 2) and sequences of barcodes for read 2 (P7) (column 3).  
Example (no header, tab-separated):  

	```
	Tc2208a	AGCGCA	AGCGCA  
	Tc2235a	CTCAGC	CTCAGC  
	Tc2249a	TAGATC	TAGATC  
	```  

- For single-end sequencing:  
Two-column file, specifying samples names (column 1), and sequences of barcodes (column 2)  
Example (no header, tab-separated):  
	```
	Tc2208a	AGCGCA  
	Tc2235a	CTCAGC  
	Tc2249a	TAGATC  
	```

&nbsp;

#### *Adapter file:*  
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


### 4/ Launch the analysis

.condarc >>> l'utilisation de CONDA pour la gestion des outils (FASTQC, MULTIQC et CUTADAPT) nécessite la création d'un environnement (envs_dirs) et d'un package (pkgs_dirs). Lors du lancement du workflow Snakemake, ils sont créer par défault dans l'espace Homedir qui est soumis en générale à des restrictions d'écritures.
				   Compléter et copier ce fichier à l'endroit où se créer par défault les packages CONDA, pour rediriger la création le l'environnement et des packages CONDA vers l'espace de travail/stockage choisit 
  
** En cours de modification:**
run_DataCleaning_workflow_Slurm.sh >>> bash script to run the snakemake workflow DATA_CLEANING on a SLURM-type computing cluster. /!\ A MODIFIER /!\
											copy this bash script into the workspace next to the data to be processed. 
run_DataCleaning_workflow_SGE.sh >>> bash script to run the snakemake workflow DATA_CLEANING on a SGE-type computing cluster. /!\ A MODIFIER /!\
											copy this bash script into the workspace next to the data to be processed. 


### 5/ Expected outputs






