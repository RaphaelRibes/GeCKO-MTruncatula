#!/bin/bash

### Variables help
cluster_config_help=$(grep -e "^--cluster-config " ${WORKFLOW_PATH}/scripts/launcher_help.txt)


## Jobs ##
if [[ "$JOBS" = 1 ]] ; then
  echo -e "\nWARNING: The workflow will be run with only 1 job allowed at a time. The tasks will not be parallelized. You can increase the maximum number of jobs allowed to be run in parallel with the --jobs option.\n"
fi

## Job scheduler ##
if [[ -z "$JOB_SCHEDULER" ]] ; then
  CLUSTER_CONFIG=""
  echo -e "\nWARNING: No job scheduler was specified. The pipeline will be run without submitting any job and any parallelization.\n"
elif [[ "$JOB_SCHEDULER" != "SLURM" && "$JOB_SCHEDULER" != "SGE" ]] ; then
  echo -e "\nERROR: The provided job scheduler is not implemented. Implemented job schedulers are 'SLURM' and 'SGE'. To run the pipeline without any job scheduler, skip the --job-scheduler option."
  echo -e "\nExiting.\n"
  exit 1
else
  PROFILE="--profile ${workflow_folder}/PROFILES/$JOB_SCHEDULER"
fi


## CLUSTER_CONFIG ##

if [[ (-z "$CLUSTER_CONFIG" || "$CLUSTER_CONFIG" = --* || "$CLUSTER_CONFIG" = -*) && ! -z "$JOB_SCHEDULER" ]] ; then
  if [[ -f "CONFIG/cluster_config_${WORKFLOW}.json" ]] ; then
    CLUSTER_CONFIG=CONFIG/cluster_config_${WORKFLOW}.json
    echo -e "\nINFO: The cluster_config_${WORKFLOW}.json file was automatically found in ${HERE}/CONFIG and will be used as the cluster_config file for this workflow."
    echo -e "If you would rather use another cluster_config file, please provide it with --cluster-config\n"
	else
    echo -e "\nERROR: You provided a job scheduler, but your cluster config file cannot be found. Please provide your cluster config file with --cluster-config or place it in a CONFIG folder under the name cluster_config_${WORKFLOW}.json."
		echo "As a reminder:"
		echo $cluster_config_help
		echo -e "\nExiting.\n"
    exit 1
  fi
fi


if [[ ! -z "$CLUSTER_CONFIG" && ! -f "$CLUSTER_CONFIG" ]] ; then
    echo -e "\nERROR: the file given in the --cluster-config parameter is not valid. Please make sure the file exists and the path is correctly written."
		echo "As a reminder:"
		echo $cluster_config_help
		echo -e "\nExiting.\n"
    exit 1
fi

if [[ ! -z "$CLUSTER_CONFIG" ]] ; then
  nb_carriage_returns=$(grep -c $'\r' $CLUSTER_CONFIG)
  if [[ "$nb_carriage_returns" -gt 0 ]] ; then
    echo "Removing windows carriage returns in ${CLUSTER_CONFIG}..."
    sed -i 's/\r$//g' $CLUSTER_CONFIG
    sed -i 's/\r/\n/g' $CLUSTER_CONFIG
  fi
  CLUSTER_CONFIG=$(absolutePath $CLUSTER_CONFIG)
  CLUSTER_CONFIG_CMD="--cluster-config $CLUSTER_CONFIG"
fi


## Singularity ##
# Comment or uncomment the singularity line in the smk
#... singularity: "docker://condaforge/mambaforge"

## Check if Conda module is available
if ! [[ -x "$(command -v conda)" ]] ; then
  echo -e "\nERROR: Conda is not available. You must install it, or make it available to your working environment (eg: module load it)."
  echo "As a reminder:"
  awk '/^- Make sure Snakemake and Conda/,/^$/' ${WORKFLOW_PATH}/scripts/launcher_help.txt
  echo -e "\nExiting.\n"
  exit 1
fi



## Create logs folder ##
mkdir -p Logs_${WORKFLOW}Workflow
