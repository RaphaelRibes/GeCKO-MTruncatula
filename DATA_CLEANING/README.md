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
 
1) Prepare your input data 
2) Clone the WORKFLOW directory  
3) Prepare the CONFIG files  
4) Launch the analysis  

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

    
### 3/ Prepare the CONFIG files

The DATA_CLEANING WORKFLOW will need information about the dataset and the analysis parameters to performs its different steps.  
The configuration information is provided through two files: cluster_config_DataCleaning.json and config_DataCleaning.yml.  
If you name them exactly as written above and place them in a folder named 'CONFIG', the bash launching script will detect them automatically. Otherwise, you will have to pass them as arguments with --config and --cluster-config (see [below](####-4/-Run-Workflow_DataCleaning) for details).

**_cluster_config_SGE_DataCleaning.json_** : file allowing to configure on a SGE-type computing cluster : the type of partition, the name of the log files, etc. It is possible to define a default version and/or to adapt it to each of the snakmake workflow rules
_cluster_config_SLURM_DataCleaning.json_ : file allowing to configure on a SLURM-type computing cluster : the type of queue, the name of the log files, etc. It is possible to define a default version and/or to adapt it to each of the snakmake workflow rules

**_config_DataCleaning.yml_** : 
file containing the information and parameters that will be used by the snakemake workflow DATA_CLEANING

> General variables:
OUTPUTS_DIRNAME: Name of the directory that will contain all the workflow outputs (example: WOKFLOW_OUTPUTS )
> Input files:
FASTQ: if the sequencing is in Single_End AND the raw fastq files is multiplexed, if not leave blank: ""
              path to the folder containing fastq.gz file. 
FASTQ_R1: If the sequencing is in Paired_End AND the raw fastq files are  multiplexed, if not leave blank: ""                     
                    path to the folder containing  fastq.gz file Read1. File extension format: name_R1.fastq.gz. 
FASTQ_R2: If the sequencing is in Paired_End AND the raw fastq files are  multiplexed, if not leave blank: ""
                    path to the folder containing  fastq.gz file Read2. File extension format: name_R2.fastq.gz. 
DEMULT_DIR: If the raw fastq files are demultiplexed (sequencing Single_end or paired_end), if not leave blank: ""
                        File extension format: PairedEnd > sampleX.R1.fastq.gz and sampleX.R2.fastq.gz / SingleEnd > sampleX.fastq.gz 
BARCODE_FILE:  If the raw fastq files are multiplexed (sequencing Single_end or paired_end)with the barcodes specific to each sample , if not leave blank: ""
                              path to the barcode_file.txt. see description below.                                                   
ADAPTER_FILE: path to the adapter_file.txt. see description below.
> Demultiplexing parameters : if the raw fastq files is multiplexed, if not leave blank: ""
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


### 5/ Expected outputs






### MORE RANDOM INFO: NO NEED TO READ BELOW
The WORKFLOW directory contains 3 folders:  

ENVS
  **_conda.tools.yml_** : file containing the list of tools needed for the conda environment. This file is used in the sanemake workflow in the "conda" parameter of the rules. 

PROFILES
    SGE
        **_config.yaml_** : file containing the parameters for launching a job on an SGE cluster (qsub, partition, etc). this file is called by the **run??? script  >>> à compléter**
    SLURM
        **_config.yaml_** : file containing the parameters for launching a job on an SLURM cluster (qsub, partition, etc). this file is called by the **run??? script  >>> à compléter**

SCRIPTS
    **_demultiplex_with_cutadapt_PE.sh_** : For Paired End sequencing.  Bash script allowing to demultiplex a pair of fastq files (R1 + R2) into individual paired fastq files according to the barcode/tag, located at the P5 and P7 ends.
   ** _demultiplex_with_cutadapt_SE.sh_** : For Single End sequencing.  Bash script allowing to demultiplex  fastq file into individual fastq files according to the barcode/tag, located at the P5 ends.
    **_trimming_with_cutadapt.sh_PE_** : For Paired End sequencing.  Bash script allowing to trimming a pair of fastq files (R1 + R2 ) to remove adapters sequences, low quality sequences et short sequences. Use of tool CUTADAPT.
    **_trimming_with_cutadapt.sh_SE_** : For Single End sequencing. Bash script allowing to trimming a fastq files  to remove adapters sequences, low quality sequences et short sequences. Use of tool CUTADAPT.
    **_fastq_read_count.sh_** : bash script allowing to count reads in all fastq.gz files of folder
	
