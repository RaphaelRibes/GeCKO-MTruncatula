#!/bin/bash


### v WRITE YOUR MODULE LOADS HERE v ###

module purge
module load snakemake/5.13.0
module load anaconda/python3.8
#   module load singularity/3.6.3


### ^ WRITE YOUR MODULE LOADS HERE ^ ###


### DEFAULT ACTION VALUES
UNLOCK="FALSE"
CONDA_CREATE_ENV_ONLY="FALSE"
DRYRUN="FALSE"
DIAGRAM="FALSE"
REPORT="FALSE"
USE_CONDA="FALSE"
USE_CONDA_AND_SINGULARITY="FALSE"


### ARGUMENTS
POSITIONAL=()
while [[ $# -gt 0 ]]
do
  key="$1"

  case $key in
    --help)
    HELP="TRUE"
    shift # past argument
    ;;
    --workflow)
    WORKFLOW="$2"
    shift
    shift
    ;;
    --workflow-path)
    WORKFLOW_PATH="$2"
    shift
    shift
    ;;
    --cluster-config)
    CLUSTER_CONFIG="$2"
    shift
    shift
    ;;
    --config-file)
    CONFIG="$2"
    shift
    shift
    ;;
    --use-conda)
    USE_CONDA="TRUE"
    shift
    ;;
    --use-conda-and-singularity)
    USE_CONDA_AND_SINGULARITY="TRUE"
    shift
    ;;
    --singularity-args)
    SINGULARITY_ARGS="--singularity-args \"$2\""
    shift
    shift
    ;;
    --jobs)
    JOBS="$2"
    shift
    shift
    ;;
    --printshellcmds)
    PRINTSHELLCMDS="--printshellcmds"
    shift
    ;;
    --forceall)
    FORCEALL="--forceall"
    shift
    ;;
    --latency-wait)
    LATENCY_WAIT="$2"
    shift
    shift
    ;;
    --dryrun)
    DRYRUN="TRUE"
    shift
    ;;
    --conda-create-envs-only)
    CONDA_CREATE_ENV_ONLY="TRUE"
    shift
    ;;
    --unlock)
    UNLOCK="TRUE"
    shift
    ;;
    --report)
    REPORT="TRUE"
    REPORT_NAME="$2"
    shift
    shift
    ;;
    --diagram)
    DIAGRAM="TRUE"
    DIAGRAM_NAME="$2"
    shift
    shift
    ;;
    --job-scheduler)
    JOB_SCHEDULER="$2"
    shift
    shift
    ;;
    *)
    POSITIONAL+=("$1")
    shift
    ;;
  esac
done
set -- "${POSITIONAL[@]}"



# --------------------------------------------------------------------------------------------------------------#


### Print help
if [ "${HELP}" = "TRUE" ] ; then
  echo "The help can't be displayed yet."
  #cat ...
  exit 0
fi

# --------------------------------------------------------------------------------------------------------------#



### WORKFLOW_PATH (needed for next step)
if [[ -z "$WORKFLOW_PATH" || "$WORKFLOW_PATH" = --* || "$WORKFLOW_PATH" = -* ]] ; then
	echo -e "\nERROR: the --workflow-path parameter is missing, please include it in your command."
	echo "As a reminder:"
	echo $workflow_path_help
	echo -e "\nExiting.\n"
	exit 1
fi
if [ ! -d "$WORKFLOW_PATH" ] ; then
	echo -e "\nERROR: the path given in the --workflow-path parameter is not valid. Please make sure the directory exists and the path is correctly written."
	echo "As a reminder:"
	echo $workflow_path_help
	echo -e "\nExiting.\n"
	exit 1
fi
if [[ ! "$WORKFLOW_PATH" = /* ]] ; then
  WORKFLOW_PATH=$(readlink -f $WORKFLOW_PATH) ;
fi


### Make scripts executable
if [[ ! -d "${WORKFLOW_PATH}/SCRIPTS/" ]] ; then
  echo -e "\nERROR: No SCRIPTS/ folder was found in the provided workflow path (${WORKFLOW_PATH}). Please clone or copy the whole WORKFLOW folder from GitHub: https://github.com/BioInfo-GE2POP-BLE/CAPTURE_PIPELINES_SNAKEMAKE containing all sub-directories."
  echo -e "\nExiting.\n"
  exit 1
fi

for script in $(ls "${WORKFLOW_PATH}/SCRIPTS/") ; do
  if [[ ! -x "${WORKFLOW_PATH}/SCRIPTS/${script}" ]] ; then
    echo "Making $script executable..."
    chmod 755 "${WORKFLOW_PATH}/SCRIPTS/${script}"
  fi
done



### Check variables and paths
source "${WORKFLOW_PATH}/SCRIPTS/GenericVariablesCheck.sh"



### Check if snakemake is available
if ! [[ -x "$(command -v snakemake)" ]] ; then
  echo -e "\nERROR: Snakemake is not available. You must install it, or make it available to your working environment (eg: module load it or activate it with conda)."
  echo -e "\nExiting.\n"
  exit 1
fi



### Comment or uncomment the singularity line in the smk
#... singularity: "docker://condaforge/mambaforge"







### RUN APPROPRIATE SNAKEMAKE COMMANDS ###

## UNLOCK ##
if [ "${UNLOCK}" = "TRUE" ] ; then
  snakemake --snakefile ${WORKFLOW_PATH}/${WORKFLOW_SMK} --unlock --cluster-config ${CLUSTER_CONFIG} --configfile ${CONFIG}
fi


## CREATE CONDA ENVIRONMENT ##
if [ "${CONDA_CREATE_ENV_ONLY}" = "TRUE" ] ; then
  snakemake --snakefile ${WORKFLOW_PATH}/${WORKFLOW_SMK} $PRINTSHELLCMDS --use-conda --conda-create-envs-only --jobs $JOBS --cluster-config ${CLUSTER_CONFIG} --configfile ${CONFIG}
  exit 0
fi


## DRYRUN ##
if [ "${DRYRUN}" = "TRUE" ] ; then
  source "${WORKFLOW_PATH}/SCRIPTS/${WORKFLOW}_VariablesCheck.sh"
  snakemake --snakefile ${WORKFLOW_PATH}/${WORKFLOW_SMK} $PRINTSHELLCMDS --dryrun --dag --forceall --cluster-config ${CLUSTER_CONFIG} --configfile ${CONFIG}
  exit 0
fi


## DIAGRAM ##
if [ "${DIAGRAM}" = "TRUE" ] ; then
  snakemake --snakefile ${WORKFLOW_PATH}/${WORKFLOW_SMK} $PRINTSHELLCMDS --dryrun --dag --forceall --cluster-config ${CLUSTER_CONFIG} --configfile ${CONFIG} | dot -Tsvg > $DIAGRAM_NAME
  exit 0
fi


## REPORT ##
if [ "${REPORT}" = "TRUE" ] ; then
  snakemake --snakefile ${WORKFLOW_PATH}/${WORKFLOW_SMK} $PRINTSHELLCMDS --report $REPORT_NAME --cluster-config ${CLUSTER_CONFIG} --configfile ${CONFIG}
  exit 0
fi


## RUN WITH CONDA ##
if [ "${USE_CONDA}" = "TRUE" ] ; then
  source "${WORKFLOW_PATH}/SCRIPTS/RunWorkflow_VariablesCheck.sh"
  source "${WORKFLOW_PATH}/SCRIPTS/${WORKFLOW}_VariablesCheck.sh"
  snakemake --snakefile ${WORKFLOW_PATH}/${WORKFLOW_SMK} $PRINTSHELLCMDS $FORCEALL --latency-wait $LATENCY_WAIT --jobs $JOBS --use-conda --cluster-config ${CLUSTER_CONFIG} --configfile ${CONFIG} ${PROFILE}
  exit 0
fi


## RUN WITH CONDA + SINGULARITY ##
if [ "${USE_CONDA_AND_SINGULARITY}" = "TRUE" ] ; then

  echo "This does not work yet sorry :("
  exit 0

  source "${WORKFLOW_PATH}/SCRIPTS/RunWorkflow_VariablesCheck.sh"
  source "${WORKFLOW_PATH}/SCRIPTS/${WORKFLOW}_VariablesCheck.sh"
  snakemake --snakefile ${WORKFLOW_PATH}/${WORKFLOW_SMK} $PRINTSHELLCMDS $FORCEALL --latency-wait $LATENCY_WAIT --jobs $JOBS --use-conda --use-singularity ${SINGULARITY_ARGS} --cluster-config ${CLUSTER_CONFIG} --configfile ${CONFIG} ${PROFILE}
fi


## NONE OF THE ABOVE ##
if ( [ "${UNLOCK}" = "FALSE" ] && [ "${CONDA_CREATE_ENV_ONLY}" = "FALSE" ] && [ "${DRYRUN}" = "FALSE" ] && [ "${DIAGRAM}" = "FALSE" ] && [ "${REPORT}" = "FALSE" ] && [ "${USE_CONDA}" = "FALSE" ] && [ "${USE_CONDA_AND_SINGULARITY}" = "FALSE" ] ) ; then
  echo "Nothing to be done. To print the help run ./runSnakemakeWorkflow.sh --help"
  exit 0
fi


#HELP :

# Run with conda -> ok sur Migale
#./runSnakemakeWorkflow.sh --workflow DataCleaning --workflow-path /home/jogirodolle/save/CAPTURE_PIPELINES_SNAKEMAKE/DATA_CLEANING/WORKFLOW \
#   [--cluster-config CONFIG/cluster_config.json] [--config-file CONFIG/config.yml] --use-conda --job-scheduler SGE \
#   [--printshellcmds] [--forceall] [--latency-wait 60] [--jobs 200]
#>> qsub -N RunSMK -V -b yes -q short.q -cwd "conda activate snakemake-6.9.1; ./runSnakemakeWorkflow.sh --workflow DataCleaning --workflow-path /home/jogirodolle/save/CAPTURE_PIPELINES_SNAKEMAKE/DATA_CLEANING/WORKFLOW --cluster-config /home/jogirodolle/work/DATA_TEST/MEL_R19/CONFIG/cluster_config_SGE.json --config-file /home/jogirodolle/work/DATA_TEST/MEL_R19/CONFIG/config_SGE.yml --jobs 200 --use-conda --job-scheduler SGE"

# Run with conda + singularity
#./runSnakemakeWorkflow.sh --workflow DataCleaning --workflow-path /home/jogirodolle/save/CAPTURE_PIPELINES_SNAKEMAKE/DATA_CLEANING/WORKFLOW \
#   --cluster-config CONFIG/cluster_config.json --config-file CONFIG/config.yml --jobs 200 --use-conda --use-singularity --job-scheduler SGE \
#   [--printshellcmds] [--forceall] [--latency-wait 60]
#>> qsub -N RunSMK -V -b yes -q short.q -cwd "conda activate snakemake-6.9.1; ./runSnakemakeWorkflow.sh --workflow DataCleaning --workflow-path /home/jogirodolle/save/CAPTURE_PIPELINES_SNAKEMAKE/DATA_CLEANING/WORKFLOW --cluster-config /home/jogirodolle/work/DATA_TEST/MEL_R19/CONFIG/cluster_config_SGE.json --config-file /home/jogirodolle/work/DATA_TEST/MEL_R19/CONFIG/config_SGE.yml --jobs 200 --use-conda-and-singularity --job-scheduler SGE"

# Dryrun the workflow -> ok sur Migale
#./runSnakemakeWorkflow.sh --workflow DataCleaning --workflow-path /home/jogirodolle/save/CAPTURE_PIPELINES_SNAKEMAKE/DATA_CLEANING/WORKFLOW \
#   --cluster-config CONFIG/cluster_config.json --config-file CONFIG/config.yml --dryrun
#>> conda activate snakemake-6.9.1; ./runSnakemakeWorkflow.sh --workflow DataCleaning --workflow-path /home/jogirodolle/save/CAPTURE_PIPELINES_SNAKEMAKE/DATA_CLEANING/WORKFLOW --cluster-config /home/jogirodolle/work/DATA_TEST/MEL_R19/CONFIG/cluster_config_SGE.json --config-file /home/jogirodolle/work/DATA_TEST/MEL_R19/CONFIG/config_SGE.yml --dryrun

# Run to create the conda environment only -> ok sur Migale
#./runSnakemakeWorkflow.sh --workflow DataCleaning --workflow-path /home/jogirodolle/save/CAPTURE_PIPELINES_SNAKEMAKE/DATA_CLEANING/WORKFLOW \
#   --cluster-config CONFIG/cluster_config.json --config-file CONFIG/config.yml --conda-create-envs-only
#>> qsub -N RunSMK -V -b yes -q short.q -cwd "conda activate snakemake-6.9.1; ./runSnakemakeWorkflow.sh --workflow DataCleaning --workflow-path /home/jogirodolle/save/CAPTURE_PIPELINES_SNAKEMAKE/DATA_CLEANING/WORKFLOW --cluster-config /home/jogirodolle/work/DATA_TEST/MEL_R19/CONFIG/cluster_config_SGE.json --config-file /home/jogirodolle/work/DATA_TEST/MEL_R19/CONFIG/config_SGE.yml --conda-create-envs-only"

# Unlock -> ok sur Migale
#./runSnakemakeWorkflow.sh --workflow DataCleaning --workflow-path /home/jogirodolle/save/CAPTURE_PIPELINES_SNAKEMAKE/DATA_CLEANING/WORKFLOW \
#   --cluster-config CONFIG/cluster_config.json --config-file CONFIG/config.yml --unlock
#>> conda activate snakemake-6.9.1; ./runSnakemakeWorkflow.sh --workflow DataCleaning --workflow-path /home/jogirodolle/save/CAPTURE_PIPELINES_SNAKEMAKE/DATA_CLEANING/WORKFLOW --cluster-config /home/jogirodolle/work/DATA_TEST/MEL_R19/CONFIG/cluster_config_SGE.json --config-file /home/jogirodolle/work/DATA_TEST/MEL_R19/CONFIG/config_SGE.yml --unlock

# Run to write the workflow's report -> ok sur Migale
#./runSnakemakeWorkflow.sh --workflow DataCleaning --workflow-path /home/jogirodolle/save/CAPTURE_PIPELINES_SNAKEMAKE/DATA_CLEANING/WORKFLOW \
#   --cluster-config CONFIG/cluster_config.json --config-file CONFIG/config.yml --report Pipeline_Report.html
#>> conda activate snakemake-6.9.1; ./runSnakemakeWorkflow.sh --workflow DataCleaning --workflow-path /home/jogirodolle/save/CAPTURE_PIPELINES_SNAKEMAKE/DATA_CLEANING/WORKFLOW --cluster-config /home/jogirodolle/work/DATA_TEST/MEL_R19/CONFIG/cluster_config_SGE.json --config-file /home/jogirodolle/work/DATA_TEST/MEL_R19/CONFIG/config_SGE.yml --report Pipeline_Report.html

# Run to write diagram -> ok sur Migale
#./runSnakemakeWorkflow.sh --workflow DataCleaning --workflow-path /home/jogirodolle/save/CAPTURE_PIPELINES_SNAKEMAKE/DATA_CLEANING/WORKFLOW \
#   --cluster-config CONFIG/cluster_config.json --config-file CONFIG/config.yml --diagram Pipeline_Diagram.svg
#>> conda activate snakemake-6.9.1; ./runSnakemakeWorkflow.sh --workflow DataCleaning --workflow-path /home/jogirodolle/save/CAPTURE_PIPELINES_SNAKEMAKE/DATA_CLEANING/WORKFLOW --cluster-config /home/jogirodolle/work/DATA_TEST/MEL_R19/CONFIG/cluster_config_SGE.json --config-file /home/jogirodolle/work/DATA_TEST/MEL_R19/CONFIG/config_SGE.yml --diagram Pipeline_Diagram.svg
