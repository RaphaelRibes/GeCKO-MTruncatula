### Project overview  
This repository provides several workflows to clean sequenced reads, map them to the reference of your choice, call SNPs and filter VCF file with polymorphism data. They are especially appropriate for data sequenced after target enrichment capture, but can be used for other types of sequenced reads.

You will find the different workflows in the corresponding folders:
- DATA_CLEANING
- READ_MAPPING
- VARIANT_CALLING 
- VCF_FILTERING 

Each corresponding folder has a specific README to explain the use of the workflow, its various options and settings, as well as an EXAMPLE subfolder with a small dataset and appropriate config files to test it.

We provide a bash launcher called runGeCKO.sh that can be used to easily run these different workflows on your own data and HPC environments with minimal effort. Each workflow will produce an html report (generating thanks to [multiQC](https://multiqc.info/)) summarizing key information for this step.

To execute one of the workflows, follow the steps:  


1) If it is the first time you clone a GitHub repository, you will first need to [generate your SSH key](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent#generating-a-new-ssh-key) and [add it to your GitHub account](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account)

2) Clone this repository in your environment:   
```git clone git@github.com:GE2POP/GeCKO.git```  
This will create a GeCKO directory, for example /home/user/GeCKO.

3) Copy the appropriate config and cluster_config files and adapt them to your data and cluster.  
For more information on this step, see the more detailed README placed in each workflow folder.  

4) Use the launcher script runGeCKO.sh to run the workflow.  

&nbsp;
### Environment  
These workflows rely on snakemake, which ensures reproducible and scalable data analysis and make worflows installation straightforward since the only requirement is to have snakemake and conda available on your environment. 

#### Snakemake
We chose to utilize [Snakemake](https://snakemake.readthedocs.io/en/stable/) as our workflow management system because of its many advantages. First, Snakemake ensures reproducibility by natively supporting conda, allowing to encapsulate software environments and dependencies for accurate result replication. In addition, it allows for the optimization of resource utilization and performance by intelligently distributing jobs across multiple cores or nodes, enabling efficient handling of complex pipelines. Finally, Snakemake's portability allows workflows to be executed on different systems, facilitating collaboration and flexibility in scientific projects.
We made dedicated efforts to extensively document the usage of our workflows, ensuring they are user-friendly and accessible to non-bioinformatician biologists. Our aim was to empower them to benefit from the numerous advantages of Snakemake without requiring any coding expertise on their part.

#### Conda 
The workflow will download and find the [tools it needs](#tools) through [Conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html), which means you do not need to have them installed in your working environment behorehand.  
When called for the first time, the DATA_CLEANING Snakemake workflow will download the tools' packages in a pkgs_dirs folder, and install them in a conda environment that will be stored in a .snakemake/conda folder, in the directory you called the workflow from. Every time you call the workflow from a new directory, the Conda environment will be generated again. To avoid creating the environment multiple times, which can be both time and resource-consuming, you can provide a specific folder where you want Snakemake to store all of its conda environments with the --conda-env-path option of the runGeCKO.sh launcher.  

The pkgs_dirs folder is by default common to your whole system or cluster personnal environment. Conda's standard behaviour is to create it in your home directory, in a .conda folder. If your home space is limited or if you do not have the right to write there from your cluster's nodes, you will need to tell Conda to store its packages somewhere else, thanks to a .condarc file. Place it in your home folder and specify the directory path where you want Conda to store the packages, following this example:  
```
envs_dirs:  
    - /path/to/appropriate/directory/env  
pkgs_dirs:  
    - /path/to/appropriate/directory/pkgs  
```
Alternatively, you can use a symbolic link to have the .conda folder in your home directory point to another folder where you would rather have Conda store its files.  
To do so:  
```
mv   /home/user/.conda   /path/to/appropriate/directory/.conda
ln -nfs /path/to/appropriate/directory/.conda   /home/user/.conda
```

Sometimes you will also need to create a symbolic link for your .cache folder:  
```
mv   /home/user/.cache   /path/to/appropriate/directory/.cache
ln -nfs /path/to/appropriate/directory/.cache   /home/user/.cache
```

#### Cluster's environment
SGE and Slurm job schedulers are currently supported.



&nbsp;
### Installation  

- If the script was open on a Windows system and you will execute it on a Linux system, you may need to remove windows carriage returns ('\r') with:  
```dos2unix runGeCKO.sh```  
or  
```sed -i 's/\r$//g' runGeCKO.sh ; sed -i 's/\r/\n/g' runGeCKO.sh```  

- Make sure Snakemake and Conda are available to your working environment.  
Either install them on your computer, or if you are working on a cluster, you may need to 'module load' them, or to 'conda activate' them, depending on your cluster's software management policy.  
    - For clusters using module environment, you can add the 'module load' lines inside runGeCKO.sh: you will find a dedicated zone "WRITE YOUR MODULE LOADS HERE" at the top of the script. It is advised to precede it with 'module purge' to avoid potential conflicts with previously loaded modules. To find out the exact name of the needed modules, use the 'module avail' command. The modules will be loaded every time you execute the script.  
    - For clusters using Conda environment, Conda will likely be readily available, and you will only need to conda activate Snakemake. To find out the precise name of the snakemake environment, use the 'conda info --envs' command. You may need to call conda activate outside of the script itself.  

&nbsp;
### Using the workflow

The following section describes the different workflow actions and parameters. The /home/user/GeCKO path refers to the directory you cloned from GitHub.

&nbsp;
#### QUICK START:  
There are only two mandatory options: one specifying the WORKFLOW directory, and another to provide the name of the workflow you want to run. So to demultiplex and trim your reads simply type:

```bash runGeCKO.sh --workflow-path /home/user/GeCKO --workflow DataCleaning```


The information regarding the fastq files, read index etc. are, by default, retrieved from the config file CONFIG/config_WorkflowName.yml. The same folder can also contain the cluster_config_WorkflowName.yml file used by default to provide specific cluster information (e.g. job queue names) related to this workflow.

To use the full resource of my HPC environment (Slurm), and allow up to 100 submitted jobs at the same time, it thus suffices to adapt this cluster config file and to type the following command:  

```bash runGeCKO.sh --workflow-path /home/jgirodolle/GeCKO --workflow DataCleaning --job-scheduler SLURM --jobs 100```  

&nbsp;

#### POSSIBLE ACTIONS:  
The launcher's default behavior is to run the workflow, but other actions can be called instead:

**--help**&nbsp;&nbsp;&nbsp;*print the help*  
```bash runGeCKO.sh --workflow-path /PATH/TO/GeCKO --help```  

**--dryrun**&nbsp;&nbsp;&nbsp;*only dryrun the workflow (and detect potential errors) without actually running it*  
```bash runGeCKO.sh --workflow-path /PATH/TO/GeCKO --workflow WorkflowName --dryrun```  

**--report**&nbsp;&nbsp;&nbsp;*write an html report of the workflow's last run*  
```bash runGeCKO.sh --workflow-path /PATH/TO/GeCKO --workflow WorkflowName --report```  

**--diagram**&nbsp;&nbsp;&nbsp;*write an svg diagram of the workflow*  
```bash runGeCKO.sh --workflow-path /PATH/TO/GeCKO --workflow WorkflowName --diagram```  

&nbsp;

#### MANDATORY PARAMETERS FOR ALL ACTIONS (except --workflow for --help):  
**--workflow-path [...]**&nbsp;&nbsp;&nbsp;*the path to the directory you cloned from GitHub*  
If the directory was cloned from GitHub, it should end with /GeCKO  

**--workflow [...]**&nbsp;&nbsp;&nbsp;*name of the workflow you want to run*  
Current existing options are 'DataCleaning', 'ReadMapping', 'VariantCalling' and 'VcfFiltering'  

&nbsp;

#### CONFIGURATION PARAMETERS:  
**--job-scheduler [...]**&nbsp;&nbsp;&nbsp;*name of the job scheduler that is installed on your cluster*  
Current supported options are 'SLURM' and 'SGE'. If omitted, the workflow will run without submitting any job and any parallelization.  

**--cluster-config [...]**&nbsp;&nbsp;&nbsp;*path to the cluster config file*  
If omitted, this script will look for a cluster_config_WorkflowName.yml file (eg: cluster_config_DataCleaning.yml) in a CONFIG/ folder in the directory it was executed from. This argument can also be absent if the job-scheduler is not specified (no jobs submitted).  

**--config-file [...]**&nbsp;&nbsp;&nbsp;*path to the workflow's config file*  
If omitted, this script will look for a config_WorkflowName.yml file (eg: config_DataCleaning.yml) in a CONFIG/ folder in the directory it was executed from.  

&nbsp;

#### USEFUL EXTRA OPTIONS:  
**--jobs [int]**&nbsp;&nbsp;&nbsp;*maximum number of jobs that can be run in parallel (default: 1)*  

**--latency-wait [int]**&nbsp;&nbsp;&nbsp;*number of seconds Snakemake will wait after the end of a task to look for the expected output files (default: 20)*  
You should increase it if Snakemake returns an error such as:  

    MissingOutputException in ... of ...:  
    Job Missing files after 20 seconds:  
    ...  
    This might be due to filesystem latency. If that is the case, consider to increase the wait time with --latency-wait.
    
**--forceall**&nbsp;&nbsp;&nbsp;*run all the workflow's steps, even if output files already exist (and overwrite them)*  
Without this argument, Snakemake's default behavior is to only run the steps for which output files are missing.  

**--conda-env-path**&nbsp;&nbsp;&nbsp;*the path to a folder of your choice where conda will build the workflow's environment*  
If specified, the environment will only be built once in this folder (the first time you run the workflow). Otherwise, a new environment will be built in every new folder the workflow is launched from.

&nbsp;

#### ANY OTHER OPTIONS:  
**--extra-snakemake-options ["..."]**&nbsp;&nbsp;&nbsp;*any list of other Snakemake options that you may want to pass to the Snakemake command*  
Be careful to provide them between quotes. For an exhaustive list of Snakemake options see https://snakemake.readthedocs.io/en/stable/index.html.  




