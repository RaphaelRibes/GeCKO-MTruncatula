### Project overview  
This repository provides several workflows to clean sequenced reads, map them to the reference of your choice and call SNPs. They are especially appropriate for data sequenced after target enrichment capture, but can be used for other types of sequenced reads.

You will find the different workflows in the corresponding folders:
- DATA_CLEANING
- READS_MAPPING
- SNP_CALLING  

These workflows rely on [snakemake](https://snakemake.readthedocs.io/en/stable/), which ensures reproducible and scalable data analysis and make pipelines installation straightforward since the only requirement is to have snakemake and [conda](https://docs.conda.io/en/latest/) available on your environment. Each pipeline produces an html report (generating thanks to [multiQC](https://multiqc.info/)) which summarize key information for this step. Finally, we provide a bash launcher called runSnakemakeWorkflow.sh that can be used to easily run these different workflows on your own data and HPC environments with minimal effort. SGE and Slurm job schedulers are currently supported.

To execute one of the workflows, follow the steps:  

1) Clone or copy the WORKFLOW folder (and all its contents) corresponding to the workflow you want to run.  
e.g. to clone all our workflows folders:  
```git clone git@github.com:BioInfo-GE2POP-BLE/CAPTURE_PIPELINES_SNAKEMAKE.git```  

2) Copy the appropriate config and cluster_config files and adapt them to your data and cluster.  
For more information on this step, see the more detailed README placed in each pipeline folder.  
It is advised to place these files in a CONFIG folder in your working directory, under the names config_WorkflowName.yml and cluster_config_WorkflowName.json  

3) Use the launcher script runSnakemakeWorkflow.sh to run the workflow.  

&nbsp;
### Installation  
- Make the script executable with:  
```chmod 755 runSnakemakeWorkflow.sh```  

- If the script was open on a Windows system and you will execute it on a Linux system, you may need to remove windows carriage returns ('\r') with:  
```dos2unix runSnakemakeWorkflow.sh```  
or  
```sed -i 's/\r$//g' runSnakemakeWorkflow.sh ; sed -i 's/\r/\n/g' runSnakemakeWorkflow.sh```  

- Make sure Snakemake and Conda are available to your working environment.  
Either install them on your computer, or if you are working on a cluster, you may need to 'module load' them, or to 'conda activate' them, depending on your cluster's software management policy.  
    - For cluster using module environment. You can add the 'module load' lines inside runSnakemakeWorkflow.sh: you will find a dedicated zone "WRITE YOUR MODULE LOADS HERE" at the top of the script. It is advised to precede it with 'module purge' to avoid potential conflicts with previously loaded modules. To find out the exact name of the needed modules, use the 'module avail' command. The modules will be loaded every time you execute the script.  
    - For cluster using Conda environment. Conda will likely be readily available, and you will only need to conda activate Snakemake. To find out the precise name of the snakemake module, use the 'conda info --envs' command. You may need to call conda activate outside of the script itself.  

&nbsp;
### Using the pipelines

The following section describes the different pipeline actions and parameters. It what follows /PATH/TO/SMK_DIR refer to the directory containing the snakemake workflow you wish to use, e.g something like /home/vranwez/CAPTURE_PIPELINES_SNAKEMAKE/.

&nbsp;
#### QUICK STARTER:  
There are only to mandatory option, one specifying where the snakemake pipeline folder is, the other is the name of the pipeline you want to run. So to demultiplex and trim your reads simply type:

```./runSnakemakeWorkflow.sh --workflow-path /PATH/TO/SMK_DIR --workflow DataCleaning```

on my computer this will look like:

```./runSnakemakeWorkflow.sh --workflow-path /home/vranwez/CAPTURE_PIPELINES_SNAKEMAKE/ --workflow DataCleaning```

The information regarding the fastq files, read index etc. are, by default, retrieved from the config file /PATH/TO/SMK_DIR/Pipeline/CONFIG/cluster_config_Pipeline.json  (/home/vranwez/CAPTURE_PIPELINES_SNAKEMAKE/DataCleaning/CONFIG/config_DataCleaning.json in this exemple). The same folder also contain the file (cluster_config_DataCleaning.json) used by default to provide specific cluster information (e.g. job queue names) related to this pipeline. 

To use the full resource of my HPC environment (slurm), it thus suffices to adapt this cluster config file and to type the following command

```./runSnakemakeWorkflow.sh --workflow-path /home/vranwez/CAPTURE_PIPELINES_SNAKEMAKE/ --workflow DataCleaning --job-scheduler SLURM```  


#### POSSIBLE ACTIONS:  

**--help**&nbsp;&nbsp;&nbsp;*print the help*  
```./runSnakemakeWorkflow.sh --workflow-path /PATH/TO/SMK_DIR --help```  

**--dryrun**&nbsp;&nbsp;&nbsp;*only dryrun the workflow (and detect potential errors) without actually running it*  
```./runSnakemakeWorkflow.sh --workflow-path /PATH/TO/SMK_DIR --workflow WorkflowName --dryrun```  

**--unlock**&nbsp;&nbsp;&nbsp;*unlock the folder*  
This must be used in case a previous Snakemake run was ended abruptly, and Snakemake returns an error such as:  

    Error: Directory cannot be locked. Please make sure that no other Snakemake process is trying to create the same files in the following directory:  
    ...  
    If you are sure that no other instances of snakemake are running on this directory, the remaining lock was likely caused by a kill signal or a power loss. It can be removed with the --unlock argument.  

It can be combined with other actions. For example, to unlock and then run the workflow:  
```./runSnakemakeWorkflow.sh --workflow-path /PATH/TO/SMK_DIR --workflow DataCleaning --unlock --run-with-conda```  

**--report**&nbsp;&nbsp;&nbsp;*write an html report of the workflow's last run*  
```./runSnakemakeWorkflow.sh --workflow-path /PATH/TO/SMK_DIR --workflow DataCleaning --report```  

**--diagram**&nbsp;&nbsp;&nbsp;*write an svg diagram of the workflow*  
```./runSnakemakeWorkflow.sh --workflow-path /PATH/TO/SMK_DIR --workflow DataCleaning --diagram```  

**--run-with-conda**&nbsp;&nbsp;&nbsp;*actually run the workflow*  
Conda is used to manage the tools used in the different steps, and a conda environment will be built every time you run a new workflow in a new folder.  
```./runSnakemakeWorkflow.sh --workflow-path /PATH/TO/SMK_DIR --workflow DataCleaning --run-with-conda```  

**--conda-create-envs-only**&nbsp;&nbsp;&nbsp;*only create the workflow's conda environment without running the workflow*  
```./runSnakemakeWorkflow.sh --workflow-path /PATH/TO/SMK_DIR --workflow DataCleaning --conda-create-envs-only```  


&nbsp;
#### MANDATORY PARAMETERS FOR ALL ACTIONS (except --workflow for --help):  
**--workflow-path [...]**&nbsp;&nbsp;&nbsp;*path to the directory containing the workflow's snakefile (.smk)*  
If the directory was cloned from GitHub, it should end with /WORKFLOW)  

**--workflow [...]**&nbsp;&nbsp;&nbsp;*name of the workflow you want to run*  
Current existing options are 'DataCleaning' and 'ReadsMapping'  

&nbsp;
#### CONFIGURATION PARAMETERS:  
**--job-scheduler [...]**&nbsp;&nbsp;&nbsp;*name of the job scheduler that is installed on your cluster*  
Current supported options are 'SLURM' and 'SGE'. If omitted, the workflow will run without submitting any job and any parallelization.  

**--cluster-config [...]**&nbsp;&nbsp;&nbsp;*path to the cluster config file*  
If omitted, this script will look for a cluster_config_WorkflowName.json file (eg: cluster_config_DataCleaning.json) in a CONFIG/ folder in the directory it was executed from. This argument can also be absent if the job-scheduler is not specified (no jobs submitted).  

**--config-file [...]**&nbsp;&nbsp;&nbsp;*path to the workflow's config file*  
If omitted, this script will look for a config_WorkflowName.json file (eg: config_DataCleaning.json) in a CONFIG/ folder in the directory it was executed from.  

&nbsp;
#### EXTRA OPTIONS FOR --run-with-conda:  
**--jobs [int]**&nbsp;&nbsp;&nbsp;*maximum number of jobs that can be run in parallel (default: 1)*  

**--printshellcmds**&nbsp;&nbsp;&nbsp;*print the shell commands run by Snakemake for each step*  

**--latency-wait [int]**&nbsp;&nbsp;&nbsp;*number of seconds Snakemake will wait after the end of a task to look for the expected output files (default: 20)*  
You should increase it if Snakemake returns an error such as:  

    MissingOutputException in ... of ...:  
    Job Missing files after 20 seconds:  
    ...  
    This might be due to filesystem latency. If that is the case, consider to increase the wait time with --latency-wait.
    
**--forceall**&nbsp;&nbsp;&nbsp;*run all the workflow's steps, even if output files already exist (and overwrite them)*  
Without this argument, Snakemake's default behavior is to only run the steps for which output files are missing.  

&nbsp;
#### ANY OTHER OPTIONS:  
**--extra-snakemake-options ["..."]**&nbsp;&nbsp;&nbsp;*any list of other Snakemake options that you may want to pass to the Snakemake command*  
Be careful to provide them between quotes. For an exhaustive list of Snakemake options see https://snakemake.readthedocs.io/en/stable/index.html.  

