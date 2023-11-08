# READ MAPPING

This READ_MAPPING workflow generates bams files from demultiplexed cleaned sequences.

It can be used to process:  
- Single-end sequences (SE), sequenced from only one end of each DNA fragment.  
- Paired-end sequences (PE), sequenced from both ends of each DNA fragment.  

### The READ_MAPPING workflow's steps
1) An index of the provided reference is created if it does not exist yet.
2) Reads are mapped to the reference, the resulting bams are sorted, the duplicates are removed if needed, the reads are filtered as specified by the user, and the final bams are indexed.
3) Total and mapped reads are counted at each step, and a MultiQC report is created, showing the reads numbers and quality after the mapping step.
4) [Optional] For each genomic region provided in a bed file (or automatically determined by the workflow), reads are counted in each sample.
5) [Optional] A sub-reference corresponding to the bed regions is produced. The reads that mapped to these regions are extracted and re-mapped to the sub-reference. The resulting sub-bams are sorted, filtered as specified by the user, and the final sub-bams are indexed.
6) [Optional] Total and mapped reads are counted before and after the second filtering step, and a MultiQC report is created, showing the reads numbers and quality after remapping the reads to the sub-reference.



## QUICK START

Needed files:  
- the full GeCKO/ folder, including the runGeCKO.sh launcher  
- your demultiplexed and trimmed fastq.gz files
- a reference file in fasta format to map your reads to
- (optional) a bed file listing genomic regions of interest 
- the cluster_config_ReadMapping.yml (in case you work on a cluster) and config_ReadMapping.yml files  


&nbsp;

To easily launch the workflow, use our runGeCKO.sh launcher. For example, to launch the workflow on our example PAIRED_END DEMULT_TRIM dataset on a Slurm job-scheduler, run the following command from the EXAMPLE/PAIRED_END directory:  
```../../../runGeCKO.sh --workflow ReadMapping --workflow-path ../../../../GeCKO --config-file CONFIG/config_ReadMapping.yml --cluster-config CONFIG/cluster_config_ReadMapping_SLURM.yml --jobs 20 --job-scheduler SLURM```  

To launch it on your own data, if you cloned the repository in /home/user and placed your config_ReadMapping.yml and cluster_config_ReadMapping.yml files in a CONFIG folder:  
```WORKFLOW_PATH=/home/user/GeCKO```  
```${WORKFLOW_PATH}/runGeCKO.sh --workflow ReadMapping --workflow-path ${WORKFLOW_PATH} --jobs 50 --job-scheduler SLURM``` 


&nbsp;

![](https://github.com/GE2POP/GeCKO/blob/main/readme_img/ReadMapping\_4elements.png)

&nbsp;


## How to use the READ_MAPPING workflow
 
1) [Clone the GitHub repository](#1-clone-the-github-repository)  
2) [Prepare your input data](#2-prepare-your-input-data) 
3) [Prepare the CONFIG files](#3-prepare-the-config-files)  
4) [Launch the analysis](#4-launch-the-analysis)  
5) [Expected outputs](#5-expected-outputs)

### 1/ Clone the GitHub repository

Follow the procedure described [here](https://github.com/GE2POP/GeCKO/tree/main#readme) to have GeCKO ready to process your data.

### 2/ Prepare your input data

The input data must be sequences from an Illumina sequencer (Miseq/Hiseq/Novaseq/...). 

Input sequences can be:  
- single-end sequences (SE): you must provide fastq files named in the format \*name\*.fastq.gz  
- paired-end sequences (PE): you must provide pairs of fastq files named in the format \*name\*.R1.fastq.gz and \*name\*.R2.fastq.gz  


### 3/ Prepare the config files

The READ_MAPPING workflow will need information about the dataset and the analysis parameters to perform its different steps.  
These information are provided through two files: *cluster_config_ReadMapping.yml* and *config_ReadMapping.yml*.  
If you name them exactly as written above and place them in a folder named 'CONFIG', the bash launching script will detect them automatically. Otherwise, you will have to pass them as arguments with --config and --cluster-config (see [below](#4-launch-the-analysis) for details).

#### *A/ The cluster_config_ReadMapping.yml file:*
This file will be needed if you run the workflow on a computer cluster and want Snakemake to submit jobs. You will <ins>only need to modify two things: the partitions or queues names</ins> to match those of your cluster, and <ins>the memory to be requested for each submitted job</ins>. The first section of the file gives the default values for the job-scheduler's parameters that Snakemake should use for all its steps (or rules). The following sections correspond to specific Snakemake steps, with new parameters values to overwrite the defaults. If you want to assign a different partition/queue or memory requirement for a specific step that does not yet have its own section, you can create a new section for it:  

	specificStepName:
    	q or partition: {partition name for specificStep}
    	mem-per-cpu or h_vmem: {needed memory for each job submitted in specificStep}

You will find [the list of the steps names](#list-of-the-snakefile-rules) along with what they do and the tools they use at the end of this page.  
Our workflows support SGE and Slurm job-schedulers. <ins>You will find cluster-config files for both in the EXAMPLE/.../CONFIG folders</ins>.  


#### *B/ The config_ReadMapping.yml file:*  
This file is used to pass all the information and tools parameters that will be used by the READ_MAPPING workflow. The workflow expects it to contain a specific list of variables and their assigned values, organized in YAML format. Expected variables are:  

**GENERAL VARIABLES**  
- *PAIRED_END:*&nbsp;&nbsp;&nbsp;Whether your data is paired-end or single-end. [TRUE or FALSE]  
- *CREATE_SUB_BAMS:*&nbsp;&nbsp;&nbsp;Whether to extract reads from regions of interest (listed in bed file) and to create the corresponding sub-bams. Cannot be set to TRUE if the BED variable is left blank. [TRUE or FALSE]  
- *MAPPING_SUBFOLDER:*&nbsp;&nbsp;&nbsp;If you want to separate results from different mapping parameters (different reference, mapping options...), provide a name for an extra folder to create in the READ_MAPPING output folder. Otherwise leave blank ("").  

**INPUT FILES**  
- *TRIM_DIRS:*&nbsp;&nbsp;&nbsp;The path(s) to the directory or directories containing the trimmed fastq files to be mapped. If left blank, the workflow will assume the fastq files are in WORKFLOWS_OUTPUTS/DATA_CLEANING/DEMULT_TRIM, which is the path to our DATA_CLEANING workflow output files. To provide several directories, separate them with spaces, e.g.: "/home/user/trim_dir1 /home/user/trim_dir2". Be careful to provide them between quotes.   
- *REFERENCE:*&nbsp;&nbsp;&nbsp;The path to the reference file in fasta format (must end with .fa, .fas or .fasta).  

*If you set CREATE_SUB_BAMS to TRUE, you either have to provide a bed file:*
- *BED:*&nbsp;&nbsp;&nbsp;The path to the bed file listing regions of interest to count reads in. Optional: can be left blank ("").  

*Or to provide ALL of the three following parameters to automatically create a bed file containing the genomic regions with enough coverage in your dataset:*
- *BED_MIN_MEAN_COV:*&nbsp;&nbsp;&nbsp;The minimum mean depth per sample (number of reads per base) to keep a genomic region. Optional: can be left blank ("").
- *BED_MIN_DIST:*&nbsp;&nbsp;&nbsp;The minimum distance allowed between two regions. If several regions are separated by a smaller distance than this, they will be merged into a single one. Optional: can be left blank ("").
- *BED_MIN_LENGTH:*&nbsp;&nbsp;&nbsp;The minimum length to keep a region after merging. Optional: can be left blank ("").

**MAPPING PARAMETERS**  
- *MAPPER:*&nbsp;&nbsp;&nbsp;The name of the mapper you want to use. Currently implemented options are 'bwa-mem2_mem', 'bwa_mem', 'bowtie2' and 'minimap2'.  
- *REMOVE_DUP:*&nbsp;&nbsp;&nbsp;Whether or not to remove duplicates after the mapping step. They will be marked either way. [TRUE or FALSE]  
- *SEQUENCING_TECHNOLOGY:*&nbsp;&nbsp;&nbsp;The name of the sequencing technology (eg: "ILLUMINA"), which will appear in the reads names after mapping: 'PL:{SEQUENCING_TECHNOLOGY}')  
- *EXTRA_MAPPER_OPTIONS:*&nbsp;&nbsp;&nbsp;Any list of options you would like to pass to the mapper command. Be careful to provide them between quotes. 
- *MAPPING_CPUS_PER_TASK:*&nbsp;&nbsp;&nbsp;The number of CPUs to allocate for each mapping task. Set to 1 if you are not working on a computing cluster. Be careful to never use quotes around this number.  
- *PICARD_MARKDUPLICATES_OPTIONS:*&nbsp;&nbsp;&nbsp;Any list of options you would like to pass to the 'picard MarkDuplicates' command. Be careful to provide them between quotes.  
- *PICARD_MARKDUPLICATES_JAVA_OPTIONS:*&nbsp;&nbsp;&nbsp;Java options to pass to the 'picard MarkDuplicates' command. Eg: "-Xmx4G".  
- *SAMTOOLS_INDEX_OPTIONS:*&nbsp;&nbsp;&nbsp;Any list of options you would like to pass to the 'samtools index' command. Be careful to provide them between quotes. For example, you may need to pass the "-c" option if you need to map your reads to a very big reference file.  

**BAM FILTERING**
- *SAMTOOLS_VIEW_FILTERS1:*&nbsp;&nbsp;&nbsp;Any list of filters you would like to pass to the 'samtools view' command after the mapping step. Eg for CREATE_SUB_BAMS set to TRUE: "-F 256 -F 2048 -f 2" to discard unproperly paired, non primary and supplementary reads before the extraction step. Eg for CREATE_SUB_BAMS set to FALSE: "-q 30" to discard alignments with a mapping quality (MAPQ score) lower than 30.   
- *SAMTOOLS_VIEW_FILTERS2:*&nbsp;&nbsp;&nbsp;Any list of filters you would like to pass to the 'samtools view' command after the remapping step. Eg: "-q 30" to discard alignments with a mapping quality (MAPQ score) lower than 30. If CREATE_SUB_BAMS is set to FALSE, leave blank (the command will not be called anyway).

See the samtools view documentation [here](https://www.htslib.org/doc/samtools-view.html) to read more about the -F, -f, -q options.

&nbsp;

<ins>Examples of config_ReadMapping.yml files can be found in the EXAMPLE/.../CONFIG folders</ins>.  

&nbsp;

#### *Filling in the config file if you want to extract and remap reads from specific zones (CREATE_SUB_BAMS set to TRUE)*
- Bed file:  
<ul>
In case you <ins>provide your own bed file</ins>, we strongly recommend merging proximal regions. This will prevent the potential issue of initially properly paired reads being remapped into distinct zones, which would result in them being erroneously flagged as improperly paired by the mapping software during the remapping step.  

Please note that any overlapping regions present in your BED file will be automatically merged into a single contiguous region. The resultant BED file (user_clean.bed), featuring the merged zones, will be available in the WORKFLOWS_OUTPUTS/READ_MAPPING/*/EXTRACTED_BAMS/REFERENCE_zones/ output folder.
</ul>
<ul>
Should you prefer to have the workflow <ins>automatically identify regions of interest</ins>, it will compute the read depth across the genome for all samples, and only retain zones that exceed the BED_MIN_MEAN_COV mean depth threshold. Subsequently, regions within a distance less than BED_MIN_DIST will be merged together, and any resulting regions shorter than BED_MIN_LENGTH will be excluded.  

We advise setting the BED_MIN_DIST parameter to a value exceeding the difference between the insert size and twice the read sequencing length. For instance, with DNA fragments averaging 500bp and both R1 and R2 reads being 150bp, a BED_MIN_DIST value greater than 200 is recommended.
</ul>

- Filtering:  
Bam filtering options are available after both the initial mapping (SAMTOOLS_VIEW_FILTERS1) and remapping (SAMTOOLS_VIEW_FILTERS2) steps. To ensure the integrity of read extraction and remapping, it is best to remove unproperly paired reads, as well as non-primary and supplementary alignments, after the first mapping. This precaution helps prevent the misclassification of improperly paired reads as singletons if only one read of a pair is preserved. Furthermore, in the absence of the primary read (in case it mapped out of the extracted zones), secondary or supplementary alignments could be incorrectly designated as primary during the remapping step. Such misclassifications can lead to the erroneous interpretation of mapping quality, resulting in inaccuracies in downstream analysis and the potential overestimation of certain reads' reliability. To filter out unproperly paired, non primary and supplementary reads, set SAMTOOLS_VIEW_FILTERS1 to "-F 256 -F2048 -f2" (see [here](https://broadinstitute.github.io/picard/explain-flags.html) to understand sam flags and [here](https://www.htslib.org/doc/samtools-view.html) for more information on samtools view's -F and -f options).

### 4/ Launch the analysis

**Environment**  
You can run this workflow on a computer or on a computer cluster. You will need Snakemake and Conda to be available.

**Launching**  
To launch the READ_MAPPING workflow, use the launching script runGeCKO.sh with the option --workflow ReadMapping:  
```WORKFLOW_PATH=/home/user/GeCKO```  
```${WORKFLOW_PATH}/runGeCKO.sh --workflow ReadMapping --workflow-path ${WORKFLOW_PATH} --jobs 50 --job-scheduler SLURM``` 

For more help on how to use the launcher, see GeCKO's general [README](https://github.com/GE2POP/GeCKO/tree/main#quick-start), or run:  
```${WORKFLOW_PATH}/runGeCKO.sh --help --workflow-path ${WORKFLOW_PATH}```  



### 5/ Expected outputs  

This workflow will create a "READ_MAPPING" directory in the "WORKFLOWS_OUTPUTS" directory. This directory is structured as follows and contains:  

<img src="https://github.com/GE2POP/GeCKO/blob/main/readme_img/OutputsTree_ReadMapping.png" width="600"/>



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
- *multiQC_ReadMapping_Bams_Report.html (and the associated directory)*:&nbsp;&nbsp;&nbsp;Graphic report based on stats reports to visualize the number and percentages of mapped reads  

**EXTRACTED_BAMS/BAMS_ZONES directory** (if a bed file was provided and CREATE_SUB_BAMS was set to TRUE)  
- *bams files*:&nbsp;&nbsp;&nbsp;One file per sample, reads that mapped to the zones provided in the bed file, with coordinates given in the corresponding sub-reference  

**EXTRACTED_BAMS/BAMS_ZONES/REPORTS directory** (if a bed file was provided and CREATE_SUB_BAMS was set to TRUE)  
- *nb_reads_per_sample.tsv*:&nbsp;&nbsp;&nbsp;Number of reads per sample after mapping and extraction  
- *DUPLICATES directory*:&nbsp;&nbsp;&nbsp;MarkDuplicates reports (one per sample)  
- *STATS directory*:&nbsp;&nbsp;&nbsp;Samtools stats reports (one per sample)  
- *multiQC_ReadMapping_SubBams_Report.html (and the associated directory)*:&nbsp;&nbsp;&nbsp;Graphic report based on stats reports to visualize the number and percentages of mapped reads in sub-bams\*  

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
- [multiqc v1.11](https://github.com/ewels/MultiQC/releases)
- [crossmap v0.6.3](http://crossmap.sourceforge.net/)
- [bedtools v2.30.0](https://github.com/arq5x/bedtools2)

These tools are loaded in a CONDA environment from the conda-forge and bioconda channels.


##  List of the snakefile rules
Name, description and tools used for each of the snakemake workflow rules:

| **Rule name**               | **Description**                                                                                                                     | **Tools**                                                                                                                              |
|:---------------------------:|:-----------------------------------------------------------------------------------------------------------------------------------:|:--------------------------------------------------------------------------------------------------------------------------------------:|
| Index_Reference             | Creating the reference index if needed                                                                                              | bwa index // bwa-mem2 index // bowtie2-build // minimap2 -d                                                                            |
| Mapping_PairedEndFastqs     | Mapping the input fastq files onto the reference (paired end reads)                                                                 | bwa mem // bwa-mem2 mem // bowtie2 // minimap2; samtools view, sort, fixmate;                                                          |
| Mapping_SingleEndFastqs     | Mapping the input fastq files onto the reference (single end reads)                                                                 | bwa mem // bwa-mem2 mem // bowtie2 // minimap2; samtools view, sort, fixmate;                                                          |
| MarkDuplicates_Bams         | Sorting bams and removing duplicates if needed                                                                                      | picard SortSam; picard MarkDuplicates                                                                                                  |
| Index_Bams                  | Creating bams index                                                                                                                 | samtools index                                                                                                                         |
| Stats_Bams                  | Computing mapping stats                                                                                                             | samtools stats                                                                                                                         |
| Summarize_BamsReadsCount    | Summarizing reads count in bams from samtools stats output                                                                          |                                                                                                                                        |
| MultiQC_Bams                | Runing MultiQC on samtools stats output                                                                                             | MultiQC                                                                                                                                |
| Create_BamsList             | Writing list of bams files                                                                                                          |                                                                                                                                        |
| Create_BedFile              | Creating bed file containing the positions of zones                                                                                 | samtools merge; bedtools genomecov                                                                                                     |
| CountReadsZones_Bams        | Counting reads in each zone provided in the bed file for each sample                                                                | samtools bedcov                                                                                                                        |
| Create_SubReference         | Generating the reduced reference corresponding to the genomic zones provided in the bed file                                        | samtools faidx                                                                                                                         |
| Create_ChainFile            | Generating a file in chain format corresponding to the genomic zones provided in the bed file, necessary for the Extract_Reads rule | samtools faidx                                                                                                                         |
| Extract_Reads               | Extracting reads that mapped to the zones provided in the bed file, thus creating subbams                                           | crossmap                                                                                                                               |
| MarkDuplicates_Subbams      | Sorting subbams and removing duplicates if needed                                                                                   | picard SortSam; picard MarkDuplicates                                                                                                  |
| Index_Subbams               | Creating subbams index                                                                                                              | samtools index                                                                                                                         |
| Stats_Subbams               | Computing mapping stats for subbams                                                                                                 | samtools stats                                                                                                                         |
| Summarize_SubbamsReadsCount | Summarizing reads count in subbams from samtools stats output                                                                       |                                                                                                                                        |
| MultiQC_Subbams             | Runing MultiQC on samtools stats output                                                                                             | MultiQC                                                                                                                                |
| Create_SubbamsList          | Writing list of subbams files                                                                                                       |                                                                                                                                        |
| Create_RefChrSizeFile       | Creating file with the name and size of each chromosome in the genomic reference                                                    | samtools faidx                                                                                                                         |


![](https://github.com/GE2POP/GeCKO/blob/main/readme_img/ReadMapping_Workflow.jpg?raw=true)


