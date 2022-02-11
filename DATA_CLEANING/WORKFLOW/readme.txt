To run these snakemake pipelines you need a configuration files (config.yml) that provides information regarding your data localization, software choices and options as well as your computing environnement. Each examples of our EXAMPLE folder come with its configuration files (for SGE and SLURM). You have to copy one of these templates in this folder and, adapt it to your needs, before being able to run any of these snakemake pipelines.

Once you got your config.yml file, simply runs the pipeline of your choice using the ./snakemake.sh command followed by the name of the pipeline to be run e.g., to demultiplex and trim unpaired reads:
./snakemake.sh DataCleaning_unpaired.smk 
