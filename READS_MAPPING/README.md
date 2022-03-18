# READS MAPPING

This READS_MAPPING workflow generates bams files from demultiplexed cleaned sequences.

It can be used to process:  
- Single-end sequences (SE), sequenced from only one end of each DNA fragment.  
- Paired-end sequences (PE), sequenced from both ends of each DNA fragment.  

### The READS_MAPPING workflow's steps
1) An index of the provided reference is created if it does not exist yet.
2) Reads are mapped to the reference, the resulting bams are sorted, the duplicates are removed if needed, and the final bams are indexed.
3) Bams reads are counted.
4) [Optionnal] For each genomic region provided in a bed file, reads are counted in each sample. A heatmap representing this data is generated.
5) [Optionnal] The reads that mapped to these regions are extracted and sub-bams are created. A corresponding sub-reference is also produced.
6) Two MultiQC reports are created, showing the reads numbers and quality after mapping, both before and after extracting reads from regions of interest.



## QUICK START

To easily launch the workflow, use our runSnakemakeWorkflow.sh launcher:  
```./runSnakemakeWorkflow.sh --workflow ReadsMapping --workflow-path PATH/TO/CAPTURE_SNAKEMAKE_WORKFLOWS```  

Needed files:  
- the full CAPTURE_SNAKEMAKE_WORKFLOWS/ folder  
- the runSnakemakeWorkflow.sh launcher  
- your demultiplexed and trimmed fastq.gz files
- a reference in fasta format to map your reads unto
- the cluster_config_ReadsMapping.yml (in case you work on a cluster) and config_ReadsMapping.yml files in a CONFIG folder  
- a bed file listing genomic regions of interest

&nbsp;

For example, if you need to launch the workflow on our ... dataset on a Slurm job-scheduler, run the following command from the EXAMPLE/... directory:  
```./runSnakemakeWorkflow.sh --workflow ReadsMapping --workflow-path /home/jogirodolle/save/CAPTURE_PIPELINES_SNAKEMAKE --config-file CONFIG/config_DataCleaning.yml --cluster-config CONFIG/cluster_config_Slurm_DataCleaning.json --jobs 20 --job-scheduler SLURM```  


&nbsp;
