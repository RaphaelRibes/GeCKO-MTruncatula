## Project overview  
This repository provides several workflows to clean sequenced reads, map them to the reference of your choice, call variants and filter the resulting VCF file. They are especially appropriate for data sequenced after target enrichment capture, but can be used for other types of sequenced reads.

You will find the different workflows in the corresponding folders:
- DATA_CLEANING
- READ_MAPPING
- VARIANT_CALLING 
- VCF_FILTERING 

Each corresponding folder has a specific README to explain the use of the workflow, its various options and settings, as well as an EXAMPLE subfolder with a small dataset and appropriate config files to test it.

We provide a bash launcher called runGeCKO.sh that can be used to easily run these different workflows on your own data and HPC environments with minimal effort. Each workflow will produce an html report (generating thanks to [multiQC](https://multiqc.info/)) summarizing key information for this step.



&nbsp;
## Installation  

1) If it is the first time you clone a GitHub repository, you will first need to [generate your SSH key](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent#generating-a-new-ssh-key) and [add it to your GitHub account](https://docs.github.com/en/authentication/connecting-to-github-with-ssh/adding-a-new-ssh-key-to-your-github-account)

2) Clone this repository in your environment:   
```git clone git@github.com:GE2POP/GeCKO.git```  
This will create a GeCKO directory, for example /home/user/GeCKO.

3) Make the script executable with:  
```chmod u+x runGeCKO.sh```

4) GeCKO will need [Snakemake](#snakemake) and [Singularity](#singularity) to run its workflows.  
Make sure Snakemake (v7) and Singularity (≥v2.4.1) are available to your working environment. Either install them on your computer, or if you are working on a cluster, you may need to ```module load``` them, or to ```conda activate``` them, depending on your cluster's software management policy.  
⚠ Please be aware that GeCKO may encounter issues with certain versions of Snakemake v7 (it was successfully launched with versions 7.15.1 and 7.32.4), and is not yet compatible with Snakemake v8.  
    - For clusters using module environment, you can add the ```module load``` lines in runGeCKO.sh: you will find a dedicated zone "WRITE YOUR MODULE LOAD HERE" at the top of the script. It is advised to precede it with ```module purge``` to avoid potential conflicts with previously loaded modules. To find out the exact name of the Snakemake and Singularity modules on your cluster, use the ```module avail``` command. The modules will be loaded every time you execute the script.  
    - For clusters using Conda environment, you will need to conda activate Snakemake and Singularity. To find out the precise name of the environments, use the ```conda info --envs``` command. You may need to call conda activate outside of the script itself.  


&nbsp;
## Running a workflow
*The /home/user/GeCKO path refers to the directory you cloned from GitHub.*

#### Prepare your data and config files 
Copy the appropriate config and profile files and adapt them to your data and cluster.  
For more information on this step, see the more detailed README placed in each workflow folder.  

#### The runGeCKO launcher  
You can use the launcher script runGeCKO.sh to run the workflow of your choice.  
For example, the following command will demultiplex and trim your reads with the DataCleaning workflow, using the full resource of your SLURM HPC environment with up to 100 submitted jobs at the same time:  

```/home/user/GeCKO/runGeCKO.sh --workflow DataCleaning --config-file CONFIG/config_DataCleaning.yml --cluster-profile CONFIG/DC_CLUSTER_PROFILE_SLURM --jobs 100``` 

The information regarding the fastq files, read index etc. will be retrieved from the config_DataCleaning.yml config file, while the config.yaml found in the CONFIG/DC_CLUSTER_PROFILE_SLURM folder will provide information specific to your cluster (e.g. job queue names, SLURM or SGE job-sheduler, etc).

&nbsp;
> ##### POSSIBLE ACTIONS:  
> The launcher's default behavior is to run the workflow, but other actions can be called instead:
>
> **--help**&nbsp;&nbsp;&nbsp;*print the help*  
> ```./runGeCKO.sh --help```
> 
> **--dryrun**&nbsp;&nbsp;&nbsp;*only dryrun the workflow (and detect potential errors) without actually running it*  
> ```./runGeCKO.sh --workflow WorkflowName --dryrun```  
>
> **--report**&nbsp;&nbsp;&nbsp;*write an html report of the workflow's last run*  
> ```./runGeCKO.sh --workflow WorkflowName --report```  
>
> **--diagram**&nbsp;&nbsp;&nbsp;*write an svg diagram of the workflow*  
> ```./runGeCKO.sh --workflow WorkflowName --diagram```  
>
> ##### MANDATORY PARAMETER FOR ALL ACTIONS except --help:  
> **--workflow [...]**&nbsp;&nbsp;&nbsp;*name of the workflow you want to run*  
> Current existing options are 'DataCleaning', 'ReadMapping', 'VariantCalling' and 'VcfFiltering'  
>
> ##### CONFIGURATION PARAMETERS:  
> **--cluster-profile [...]**&nbsp;&nbsp;&nbsp;*path to the cluster profile folder*  
> This folder must contain a yaml configuration file named 'config.yaml'. This file's first part provides information specific to your cluster (queues/partitions names, and needed cpus and memory for each task). The second part determines how jobs should be submitted depending on your cluster's job-scheduler. Please use the GeCKO's templates provided for SGE and SLURM job-shedulers in the example section of each workflow.  
> If omitted, the workflow will run without submitting any job and any parallelization. 
> 
> **--config-file [...]**&nbsp;&nbsp;&nbsp;*path to the workflow's config file*  
> If omitted, this script will look for a config_WorkflowName.yml file (eg: config_DataCleaning.yml) in a CONFIG/ folder in the directory it was executed from.  
>
> ##### USEFUL EXTRA OPTIONS:  
> **--jobs [int]**&nbsp;&nbsp;&nbsp;*maximum number of jobs that can be run in parallel (default: 1)*  
>
> **--latency-wait [int]**&nbsp;&nbsp;&nbsp;*number of seconds Snakemake will wait after the end of a task to look for the expected output files (default: 20)*  
> You should increase it if Snakemake returns an error such as:  
>    ```
>    MissingOutputException in ... of ...:  
>    Job Missing files after 20 seconds:  
>    ...  
>    This might be due to filesystem latency. If that is the case, consider to increase the wait time with --latency-wait.
>    ```
>    
> **--forceall**&nbsp;&nbsp;&nbsp;*run all the workflow's steps, even if output files already exist (and overwrite them)*  
> Without this argument, Snakemake's default behavior is to only run the steps for which output files are missing.  
>
> ##### ANY OTHER OPTIONS:  
> **--extra-snakemake-options ["..."]**&nbsp;&nbsp;&nbsp;*any list of other Snakemake options that you may want to pass to the Snakemake command*  
> Be careful to provide them between quotes. For an exhaustive list of Snakemake options see https://snakemake.readthedocs.io/en/stable/index.html.  

#### Outputs
Each of the four workflows generates its own set of files, organized in a WORKFLOW_OUTPUTS folder. One of their common features, however, is the creation of a ‘workflow_info.txt’ file. This file records the date and time of each run, the configuration files used, and the GitHub commit ID of the workflow version you’re working with. It also tracks when each output file was created. If you run the workflow again, i.e. to produce some missing files, this file will be updated with the new run's details, while keeping the previous run’s information intact.

&nbsp;
## Environment  
These workflows rely on snakemake, which ensures reproducible and scalable data analysis and make worflows installation straightforward since the only requirement is to have snakemake and conda available on your environment. 

#### Snakemake
We chose to utilize [Snakemake](https://snakemake.readthedocs.io/en/stable/) as our workflow manager. It ensures reproducibility by natively supporting conda, enabling the encapsulation of software environments and dependencies for accurate result replication. It also optimizes resource utilization and performance in a cluster computing environment by intelligently distributing jobs across multiple cores or nodes. Additionally, Snakemake's portability allows workflows to be executed on different systems, promoting collaboration and flexibility in scientific projects.   
We made dedicated efforts to extensively document the usage of our workflows, ensuring they are user-friendly and accessible to non-bioinformatician biologists. Our aim was to empower them to benefit from the numerous advantages of Snakemake without requiring any coding expertise on their part.

#### Singularity
##### Running workflows with Singularity
Our launcher executes GeCKO's workflows using [Singularity](https://docs.sylabs.io/guides/latest/user-guide/) (or [Apptainer](https://apptainer.org/docs/user/latest/), Singularity's open-source successor) via Snakemake's ```--use-singularity```
option. Singularity/Apptainer encapsulates software dependencies and environments within a single "image" file (.sif), enabling the creation containers for running applications in a consistent and isolated environment. This ensures reproducibility and portability without requiring you to pre-install the software required by GeCKO on your system.  

##### Automatic image download
On the first launch, our launcher automatically downloads the required image (GeCKO.sif, ~2GB) from [Sylabs Cloud](https://cloud.sylabs.io/library/ge2pop_gecko/gecko/gecko) and stores it in GeCKO/launcher_files/container. Once downloaded, all workflow rules are executed in a container created from this image.  



#### Cluster's environment
SGE and Slurm job schedulers are currently supported.

