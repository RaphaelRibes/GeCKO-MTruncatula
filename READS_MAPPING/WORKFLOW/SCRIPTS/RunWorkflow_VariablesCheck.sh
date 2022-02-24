#!/bin/bash


if [[ "$JOBS" = 1 ]] ; then
  echo -e "\nWARNING: The workflow will be run with only 1 job allowed at a time. The tasks will not be parallelized. You can increase the maximum number of jobs allowed to be run in parallel with the --jobs option.\n"
fi

if [[ -z "$JOB_SCHEDULER" ]] ; then
  echo -e "\nWARNING: No job scheduler was provided. The pipeline will be run without submitting any job and any parallelization.\n"
elif [[ "$JOB_SCHEDULER" != "SLURM" && "$JOB_SCHEDULER" != "SGE" ]] ; then
  echo -e "\nERROR: The provided job scheduler is not implemented. Implemented job schedulers are 'SLURM' and 'SGE'. To run the pipeline without any job scheduler, skip the --job-scheduler option."
  echo -e "\nExiting.\n"
  exit 1
else
  PROFILE="--profile ${WORKFLOW_PATH}/PROFILES/$JOB_SCHEDULER"
fi

## Singularity:
# Comment or uncomment the singularity line in the smk
#... singularity: "docker://condaforge/mambaforge"


# Create logs folder
mkdir -p Logs_${WORKFLOW}Workflow
