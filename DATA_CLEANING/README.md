# DATA CLEANING

This DATA_CLEANING workflow generates demultiplexed cleaned sequences from raw sequenced data.

It can be used to process:  
- Single-end sequences (SE), sequenced from only one end of each DNA fragment.  
- Paired-end sequences (PE), sequenced from both ends of each DNA fragment.  

Besides, input sequences can be :
- Multiplexed (sequences of several samples packed in only one (SE) or a pair of (PE) fastq.gz files)
- Demultiplexed (sequences of each sample packed in separate fastq.gz files, one per sample if SE or a pair if PE)



### The DATA_CLEANING workflow's steps
Steps 1, 2, 3, 4 are done for mulitplexed data and skipped otherwise. 
1) Quality analysis of raw sequences/reads from sequencing (FASTQC)  
2) Counting the number of reads in one (SE) or two (PE) fastq.gz input files  
3) Demultiplexing one (SE) or two (PE) fastq.gz files into individual fastq files based on the barcode or tag specific to each sample (CUTADAPT)  
4) [Optional] Extracting UMI sequences from reads and storing them in read names (using umi_tools).  
5) Counting the number of reads in one (SE) or two (PE) fastq.gz files per sample  
6) Trimming one (SE) or two (PE) fastq.gz files per sample to remove adapters sequences, low quality sequences and short sequences (CUTADAPT)  
7) Counting the number of reads in one (SE) or two (PE) fastq.gz files per sample  
8) Quality analysis of each sample's reads, after demultiplexing and trimming (FASTQC)  
9) Creation of two reports (MultiQC) allowing to visualise the impact of the trimming step on the quality of the sequences (for each individual sample), and the impact of the whole workflow (for all samples merged together).  


&nbsp;


## QUICK START

Needed files:  
- the full GeCKO/ folder, including the runGeCKO.sh launcher  
- your raw fastq.gz file(s)  
- a DC_CLUSTER_PROFILE folder (in case you work on a cluster) and a config_DataCleaning.yml file  
- a barcode file in case of multiplexed data and an adapter file  

&nbsp;

To easily launch the workflow, use the runGeCKO.sh launcher. For example, to launch the workflow on the MULTIPLEXED_PAIRED_END dataset on a Slurm job-scheduler, run the following command from the EXAMPLE/MULTIPLEXED_PAIRED_END directory:  
```../../../runGeCKO.sh --workflow DataCleaning --config-file CONFIG/config_DataCleaning.yml --cluster-profile CONFIG/DC_CLUSTER_PROFILE_SLURM --jobs 20```  

To launch it on your own data, if you cloned the repository in /home/user and placed your config_DataCleaning.yml file and your DC_CLUSTER_PROFILE folder in a CONFIG folder:  
```WORKFLOW_PATH=/home/user/GeCKO```  
```${WORKFLOW_PATH}/runGeCKO.sh --workflow DataCleaning --cluster-profile CONFIG/DC_CLUSTER_PROFILE --jobs 100```  

&nbsp;
 
![](https://github.com/GE2POP/GeCKO/blob/main/readme_img/DataCleaning_4elements.png?raw=true)

&nbsp;


## How to use the DATA_CLEANING workflow

1) [Clone the GitHub repository](#1-clone-the-github-repository)  
2) [Prepare your input data](#2-prepare-your-input-data)  
3) [Prepare the CONFIG files](#3-prepare-the-config-files)  
4) [Launch the analysis](#4-launch-the-analysis)  
5) [Expected outputs](#5-expected-outputs)


### 1/ Clone the GitHub repository

Follow the procedure described [here](https://github.com/GE2POP/GeCKO/tree/main#installation) to have GeCKO ready to process your data.

### 2/ Prepare your input data

The input data must be sequences from an Illumina sequencer (Miseq/Hiseq/Novaseq/...).  

Input sequences can be multiplexed:  
- single-end sequences (SE): you must provide one fastq file named in the format \*name\*.fastq.gz  
- paired-end sequences (PE): you must provide two fastq files named in the format \*name\*_R1.fastq.gz and \*name\*_R2.fastq.gz  

or already demultiplexed<sup>1</sup>:  
- single-end sequences (SE): you must provide one fastq file per sample named in the format \*sample\*.fastq.gz  
- paired-end sequences (PE): you must provide two fastq files per sample named in the format \*sample\*.R1.fastq.gz and \*sample\*.R2.fastq.gz  
The demultiplexed sequences must be placed together in a folder to be specified in the config_file.txt.  

<sup>1</sup> *either because there was no multiplexing involved, or because the demultiplexing was performed by the sequencer (as it is the case e.g. if your libraries were multiplexed with UDI adapters)*
    
### 3/ Prepare the config files

The DATA_CLEANING workflow will need information about the dataset and the analysis parameters to perform its different steps.  
These information are provided through two files: a *config.yaml* profile file placed in a specific folder, and a *config_DataCleaning.yml* file. For the latter, if you name it exactly as written above and place it in a folder named 'CONFIG', the bash launching script will detect it automatically. Otherwise, you will have to pass it as an argument with ```--config``` (see [below](#4-launch-the-analysis) for details).

#### *A/ The PROFILE config.yaml file:*
If you intend to execute the workflow on a computer cluster and want it to run tasks in parallel, you must provide the ```--cluster-profile``` parameter with a PROFILE folder. This folder should contain a file named 'config.yaml', giving the needed information to properly submit jobs on the cluster you work on. <ins>Examples of this file, adapted to SGE and SLURM job-schedulers, are provided in the CONFIG folder for each DATA_CLEANING example dataset</ins>. Depending on your job-scheduler, pick either the DC_CLUSTER_PROFILE_SGE or the DC_CLUSTER_PROFILE_SLURM folder, and adapt the config.yaml to your needs.  
The yaml file is organized into two parts, but you will only need to modify the first one. In this first part, the first section ('default-resources') provides the default values for the cluster's resources (partitions = 'partition' and memory = 'mem_mb') that the workflow should use for all its steps (or rules). If you want to assign a different partition/queue or memory requirement for a specific step, you can specify it in the second section ('set-resources'), and it will overwrite the defaults. Finally, in the last section ('set-threads') you can provide the number of threads/CPUs needed for each step (default = 1).  

You will find [the list of the steps names](#list-of-the-snakefile-rules) along with what they do and the tools they use at the end of this page.  

&nbsp;

#### *B/ The config_DataCleaning.yml file:*  
This file is used to pass all the information and tools parameters that will be used by the DATA_CLEANING workflow. The workflow expects it to contain a specific list of variables and their assigned values, organized in YAML format. Expected variables are:  

**GENERAL VARIABLES**  
- *PAIRED_END:*&nbsp;&nbsp;&nbsp;Whether your data is paired-end or single-end [TRUE or FALSE]  

**INPUT FILES**  
- *FASTQ:*&nbsp;&nbsp;&nbsp;Path to the raw data fastq.gz file (for single-end AND multiplexed data only, otherwise leave blank: "")  
- *FASTQ_R1:*&nbsp;&nbsp;&nbsp;Path to the R1 raw data fastq.gz file, name format must be: name_R1.fastq.gz (for paired-end AND multiplexed data only, otherwise leave blank: "")  
- *FASTQ_R2:*&nbsp;&nbsp;&nbsp;Path to the R2 raw data fastq.gz file, name format must be: name_R2.fastq.gz (for paired-end AND multiplexed data only, otherwise leave blank: "")  
- *DEMULT_DIR:*&nbsp;&nbsp;&nbsp;Path to the folder containing the demultiplexed input files (for demultiplexed data only, otherwise leave blank: ""). Fastq files in the DEMULT_DIR folder must have the following name format: sampleX.R1.fastq.gz and sampleX.R2.fastq.gz for paired-end data; sampleX.fastq.gz for single-end data.  
- *BARCODE_FILE:*&nbsp;&nbsp;&nbsp;Path to the barcode file (for multiplexed data only otherwise leave blank: ""). See description [below](#barcode-file) for more details.  
- *ADAPTER_FILE:*&nbsp;&nbsp;&nbsp;Path to the adapter file. See description [below](#adapter-file) for more details.  


**DEMULTIPLEXING PARAMETERS** (mandatory if the raw fastq files are multiplexed, otherwise leave blank: "")  
- *DEMULT_SUBSTITUTIONS:*&nbsp;&nbsp;&nbsp;Fraction of authorized substitutions per barcode (tag). Example: to allow 1 substitution for an 8bp barcode, use '0.15'. (will be passed to the "--substitutions" Cutadapt parameter)    
- *CUTADAPT_DEMULT_EXTRA_OPTIONS:*&nbsp;&nbsp;&nbsp;Any list of options or parameters you would like to pass to the cutadapt command. Be careful to provide them between quotes. The '--pair-adapters' option is automatically added if PAIRED_END is set to TRUE.

**UMI EXTRACTION PARAMETERS**
- *UMI:*&nbsp;&nbsp;&nbsp;Wether or not UMI sequences should be extracted from reads. Set to TRUE if UMIs were incorporated during library construction. This option is currently only supported for demultiplexed data. [TRUE or FALSE]
- *UMITOOLS_EXTRACT_OPTIONS:*&nbsp;&nbsp;&nbsp;Any list of options or parameters you would like to pass to the '[umi-tools extract](https://umi-tools.readthedocs.io/en/latest/reference/extract.html)' command. Be careful to provide them between quotes. For example, if you expect the UMI sequences to be 8 bp long at the 5' end of your R1 reads, you should use: "--extract-method=string --bc-pattern=NNNNNNNN". See description [below](#extracting-umi-sequences) for more details.

**TRIMMING PARAMETERS** (mandatory)    
- *TRIMMING_QUAL:*&nbsp;&nbsp;&nbsp;This parameter is used to trim low-quality ends from reads. Example:  If '30': nucleotides with quality score < Q30 (1 chance out of 1000 that the sequenced base is incorrect) will be replaced by N (will be passed to the "--quality_cutoff" Cutadapt parameter)  
- *TRIMMING_MIN_LENGTH:*&nbsp;&nbsp;&nbsp;parameter to indicate the minimum size of the sequences to be kept, after applying the TRIMMING_QUAL parameter (will be passed to the "--minimum_length" Cutadapt parameter)  
- *CUTADAPT_TRIMMING_EXTRA_OPTIONS:*&nbsp;&nbsp;&nbsp;Any list of options or parameters you would like to pass to the 'cutadapt --action=trim' command. Be careful to provide them between quotes.

&nbsp;

*IN A NUTSHELL*  
**Paired-end vs single-end**  
Expected variables differ for paired-end and single-end data: if PAIRED_END is set to TRUE, both FASTQ_R1 and FASTQ_R2 are expected. If PAIRED_END is set to FALSE, only FASTQ is expected.  

**Demultiplex or multiplexed data**  
Some variables must be left blank (*variable: ""*), depending on whether the raw data is multiplexed or demultiplexed. In case of multiplexed data, DEMULT_DIR must be left blank. In case of demultiplexed data, FASTQ, FASTQ_R1, FASTQ_R2, BARCODE_FILE, DEMULT_CORES and DEMULT_SUBSTITUTIONS must be left blank.

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

This file is used for the trimming step, which removes unwanted technical sequences remaining in the reads. These sequences appear when biological fragments are shorter than expected: in such cases, the whole biological fragment is sequenced from one end to the other, and the sequencer reaches the technical adapter sequence on the 3' end of the fragment, adding its sequence (or part of it) at the end of the read. This technical part of the sequenced read then need to be removed, which is what Cutadapt does in our workflow.  

The adapter file must provide, for each sample, the expected technical sequences that may have been accidentally sequenced. They must be given in the 5'-3' reading direction, as they are expected to be sequenced following the biological sequence part of the read.  

Illumina adapters and index i5/i7 sequences are available in [this document](https://support-docs.illumina.com/SHARE/AdapterSeq/illumina-adapter-sequences.pdf).


**Example for paired-end sequencing:**
> ⚠ *The following figure is a generic representation of Illumina adapters. Some elements may differ or not be present in your data depending on your libraries preparation protocol. To properly trim your reads, you will need to know what type of adapter was used to prepare your libraries, and provide adapters sequences for each sample accordingly.*  

![fragment_structure](https://github.com/GE2POP/GeCKO/blob/main/readme_img/DataCleaning_FragmentStructure.jpg)


Technical sequence to look for and remove at the 3' end of R1 reads (given in 5'-3' reading direction):  

	[barcode7revcomp]AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC[index7]ATCTCGTATGCCGTCTTCTGCTTG  
	
Technical sequence to look for and remove at the 3' end of R2 reads (given in 5'-3' reading direction):  

	[barcode5revcomp]AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT[index5revcomp]GTGTAGATCTCGGTGGTCGCCGTATCATT  


- _For paired-end sequencing_:  
Three-column file specifying the samples names (column 1), technical sequences to look for and remove in R1 (column 2) and technical sequences to look for and remove in R2 (column 3)  
	Example (no header, tab-separated):  

	```
	Tc2208a	TGCGCTAGATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCCGCATCTCGTATGCCGTCTTCTGCTTGA	TGCGCTAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTGGCCGTATCATTA  
	Tc2235a	GCTGAGAGATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCCGCATCTCGTATGCCGTCTTCTGCTTGA	GCTGAGAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTGGCCGTATCATTA  
	Tc2249a	GATCTAAGATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCCGCATCTCGTATGCCGTCTTCTGCTTGA	GATCTAAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTGGCCGTATCATTA  
	```
	
- _For single-end sequencing_:  
Two-column file specifying the samples names (column 1) and technical sequences to look for and remove in reads (column 2)  
	Example (no header, tab-separated):  

	```
	Tc2208a	TGCGCTAGATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCCGCATCTCGTATGCCGTCTTCTGCTTGA  
	Tc2235a	GCTGAGAGATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCCGCATCTCGTATGCCGTCTTCTGCTTGA  
	Tc2249a	GATCTAAGATCGGAAGAGCACACGTCTGAACTCCAGTCACGTCCGCATCTCGTATGCCGTCTTCTGCTTGA  
	```

&nbsp;


> 📌 **_TIP_**  
> The Cutadapt software will assume 3' adapters to be ligated to the 3' end of the reads, and 5' adapters to be ligated to the 5' end of the reads. When such a sequence shows up in a read, not only is the adapter sequence trimmed, but also the sequence following it (if there is any).  
> This means that it is sufficient to provide only the beginning of the technical sequences to be removed for each sample.  
> For example, in case of UDI Truseq adapters where barcodes (b5 and b7) are absent from the construction, you can provide the first part of the adapters (instead of the whole sequences), which is the same for all samples:  
> 	```
>	Tc2208a	AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC	AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT  
>	Tc2235a	AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC	AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT  
>	Tc2249a	AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC	AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT  
>	```
>

&nbsp;


#### *Extracting UMI sequences:* 
Unique Molecular Identifiers (UMIs) are short sequences added to each DNA fragment during library construction. They help identify and distinguish between original molecules and duplicates, which can arise from PCR amplification. If UMIs were incorporated during your library construction, they need to be extracted from the reads prior to further processing steps such as read mapping. This can be done with this workflow by setting ```UMI: TRUE```. In this case, ```umi-tools```'s ```extract``` command will be used to extract the UMI sequences and move them to the read name for easy tracking. For more details about the command's options to pass to ```UMITOOLS_EXTRACT_OPTIONS``` in the config_DataCleaning.yml file, see [here](https://umi-tools.readthedocs.io/en/latest/reference/extract.html).

**Note:** At this stage, <ins>duplicate reads are not yet removed</ins>. Deduplication will be handled after mapping, using the UMI information stored in the read names. If you proceed with the ReadMapping workflow, be sure to set REMOVE_DUP_UMI: TRUE to enable duplicate removal based on UMIs.

⚠ **Warning:** Currently, ```UMI: TRUE``` can only be used with <ins>demultiplexed data</ins>. If UMI is set to TRUE for multiplexed data, UMI sequences <ins>will not</ins> be extracted.

&nbsp;

### 4/ Launch the analysis

**Environment**  
You can run this workflow on a computer or on a computer cluster. You will need Snakemake and Mamba to be available. If you chose to [create the GeCKO_env conda environment with runGeCKO.sh](https://github.com/GE2POP/GeCKO/tree/main#recommended-method), you first need to activate it:  
```conda activate GeCKO_env```  

**Launching**  
To launch the DATA_CLEANING workflow, assuming you placed your config_DataCleaning.yml and DC_CLUSTER_PROFILE folder in a CONFIG folder, use the launching script runGeCKO.sh with the option --workflow DataCleaning:  
```WORKFLOW_PATH=/home/user/GeCKO```  
```${WORKFLOW_PATH}/runGeCKO.sh --workflow DataCleaning --cluster-profile CONFIG/DC_CLUSTER_PROFILE --jobs 100``` 
 
For more help on how to use the launcher, see GeCKO's general [README](https://github.com/GE2POP/GeCKO/tree/main#quick-start), or run:  
```${WORKFLOW_PATH}/runGeCKO.sh --help```  

### 5/ Expected outputs  

This workflow will create a "DATA_CLEANING" directory in the "WORKFLOWS_OUTPUTS" directory. This directory is structured as follows and contains:  

<img src="https://github.com/GE2POP/GeCKO/blob/main/readme_img/OutputsTree_DataCleaning.png" width="600"/>

&nbsp;

<ins>Description of the main files:</ins>  

- *multiQC_DataCleaning_Report.html (and the associated directory)*:&nbsp;&nbsp;&nbsp;Graphic report based on raw data, (demultiplexing) and trimming fastQC reports to visualize the impact of DATA_CLEANING.  
- *workflow_info.txt*:&nbsp;&nbsp;&nbsp;File containing the date and time of the workflow launch, the link to the Github repository, the corresponding commit ID, and a copy of the config files provided by the user

**RAWDATA/REPORTS directory**  
- *FASTQC directory*:&nbsp;&nbsp;&nbsp;Graphics representations (R1 /R2) of the quality of the raw reads before cleaning (fastQC)  
- *Reads_Count_RawData.txt*:&nbsp;&nbsp;&nbsp;Number of reads per initial fastq.gz file (illumina sequencer)  

**DEMULT (or DEMULT_UMI) directory** (if the data is multiplexed using barcodes)  
- *fastq.gz files*:&nbsp;&nbsp;&nbsp;One pair per sample after demultiplexing and optionnal UMI extraction (sample1.R1.fastq.gz + sample1.R2.fastq.gz) and unassigned reads in "unknown" files.  

**DEMULT/REPORTS directory**  
- *Reads_Count_Demult.txt*:&nbsp;&nbsp;&nbsp;Number of reads per sample after demultiplexing  
- *CUTADAPT_INFOS directory*:&nbsp;&nbsp;&nbsp;Cutadapt demultiplexing informations, for each barcode  
- *FASTQC directory*:&nbsp;&nbsp;&nbsp;One html file (and the associated .zip file) per sample, to visualize the quality of reads after demultiplexing (fastQC)  
- *multiQC_Demult_Report.html (and the associated directory)*:&nbsp;&nbsp;&nbsp;Graphic report based on fastQC and cutadapt reports to visualize the quality of reads after demultiplexing  


**DEMULT_TRIM directory**  
- *fastq.gz files*:&nbsp;&nbsp;&nbsp;One pair per sample after trimming (sample1_trimmed.R1.fastq.gz + sample1_trimmed.R2.fastq.gz)  

**DEMULT_TRIM/REPORTS directory**  
- *Reads_Count_DemultTrim.txt*:&nbsp;&nbsp;&nbsp;Number of reads per sample after trimming  
- *CUTADAPT_INFOS directory*:&nbsp;&nbsp;&nbsp;Cutadapt trimming reports (one per sample)  
- *FASTQC directory*:&nbsp;&nbsp;&nbsp;One html file (and the associated .zip file) per sample, to visualize the quality of reads after trimming (fastQC)  
- *multiQC_Trimming_Report.html (and the associated directory)*:&nbsp;&nbsp;&nbsp;Graphic report based on fastQC and cutadapt reports to visualize the quality of reads after trimming  
	 



## Tools
This workflow uses the following tools: 
- [Cutadapt v3.5](https://cutadapt.readthedocs.io/en/v3.5/)
- [UMI-tools v1.1.5](https://umi-tools.readthedocs.io/en/latest/index.html)
- [FastQC v11.9](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) 
- [MultiQC v1.11](https://github.com/ewels/MultiQC/releases)
 
These tools are automatically downloaded from the conda-forge and bioconda channels and installed in a Conda environment by Snakemake with Mamba.

##  List of the snakefile rules
Name, description and tools used for each of the snakemake workflow rules:

| **Rule name**              | **Description**                                                                  | **Tools** |
|:--------------------------:|:--------------------------------------------------------------------------------:|:---------:|
| Fastqc_RawFastqs           | Runing FastQC on raw fastq files                                                 | FastQC    |
| CountReads_RawFastqs       | Counting reads in raw fastq files                                                |           |
| Demultiplex_RawFastqs      | Demultiplexing raw fastq files                                                   | Cutadapt  |
| ExtractUMI_DemultFastqs    | Extracting UMI sequences from reads in demultiplexed fastq files and storing them in read names    | umi_tools extract  |
| CountReads_DemultFastqs    | Counting reads in demultiplexed fastq files                                      |           |
| Fastqc_DemultFastqs        | Runing FastQC on demultiplexed fastq files                                       | FastQC    |
| MultiQC_DemultFastqs       | Runing MultiQC on demultiplexed fastq files                                      | MultiQC   |
| Concatenate_DemultFastqs   | Concatenating all demultiplexed fastq files to compute global quality statistics |           |
| Fastqc_ConcatDemultFastqs  | Runing FastQC on concatenated demultiplexed fastq files                          | FastQC    |
| Trimming_DemultFastqs      | Trimming demultiplexed fastq files                                               | Cutadapt  |
| CountReads_TrimmedFastqs   | Counting reads in trimmed fastq files                                            |           |
| Fastqc_TrimmedFastqs       | Runing FastQC on trimmed fastq files                                             | FastQC    |
| MultiQC_TrimmedFastqs      | Runing MultiQC on trimmed fastq files                                            | MultiQC   |
| Concatenate_TrimmedFastqs  | Concatenating all trimmed fastq files to compute global quality statistics       |           |
| Fastqc_ConcatTrimmedFastqs | Runing FastQC on concatenated trimmed fastq files                                | FastQC    |
| MultiQC_Global             | Runing MultiQC on raw and concatenated trimmed fastq files                       | MultiQC   |



<!-- ![](https://github.com/GE2POP/GeCKO/blob/main/readme_img/DataCleaning_Workflow.jpg?raw=true) -->

