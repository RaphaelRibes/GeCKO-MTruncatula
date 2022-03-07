#!/bin/bash


### v WRITE YOUR MODULE LOADS HERE v ###

# module purge
# module load snakemake/5.13.0
# module load anaconda/python3.8
#   (module load singularity/3.6.3)


### ^ WRITE YOUR MODULE LOADS HERE ^ ###

### DEFAULT OPTIONS
JOBS=1
LATENCY_WAIT=20

### DEFAULT ACTION VALUES
HELP="FALSE"
UNLOCK="FALSE"
CONDA_CREATE_ENV_ONLY="FALSE"
DRYRUN="FALSE"
DIAGRAM="FALSE"
REPORT="FALSE"
USE_CONDA="FALSE"
#USE_CONDA_AND_SINGULARITY="FALSE"


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
    --run-with-conda)
    USE_CONDA="TRUE"
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
    --extra-snakemake-options)
    EXTRA_SNAKEMAKE_OPTIONS="$2"
    shift
    shift
    ;;
#    --use-conda-and-singularity)
#    USE_CONDA_AND_SINGULARITY="TRUE"
#    shift
#    ;;
#    --singularity-args)
#    SINGULARITY_ARGS="--singularity-args \"$2\""
#    shift
#    shift
#    ;;
    *)
    POSITIONAL+=("$1")
    shift
    ;;
  esac
done
set -- "${POSITIONAL[@]}"


# --------------------------------------------------------------------------------------------------------------#

### Functions

absolutePath () {
  if [[ ! -z "$1" && ! "$1" = /* ]] ; then
    fileOrFolder_absolutePath=$(readlink -f $1) ;
    echo ${fileOrFolder_absolutePath%/}
  else
    echo ${1%/}
  fi
}

# --------------------------------------------------------------------------------------------------------------#

### WORKFLOW_PATH (needed for next steps)
WORKFLOW_PATH=$(absolutePath $WORKFLOW_PATH)

if [[ -z "$WORKFLOW_PATH" || "$WORKFLOW_PATH" = --* || "$WORKFLOW_PATH" = -* ]] ; then
	echo -e "\nERROR: the --workflow-path parameter is missing, please include it in your command."
	echo "As a reminder:"
	echo "--workflow-path [...]: the path to the directory containing the workflow's Snakemake (.smk) file. If the directory was cloned from GitHub, it should end with /WORKFLOW)"
	echo -e "\nExiting.\n"
	exit 1
fi
if [ ! -d "$WORKFLOW_PATH" ] ; then
	echo -e "\nERROR: the path given in the --workflow-path parameter is not valid. Please make sure the directory exists and the path is correctly written."
	echo "As a reminder:"
  echo "--workflow-path [...]: the path to the directory containing the workflow's Snakemake (.smk) file. If the directory was cloned from GitHub, it should end with /WORKFLOW)"
	echo -e "\nExiting.\n"
	exit 1
fi

# --------------------------------------------------------------------------------------------------------------#

## NO ACTION PROVIDED ##
if ( [ "${UNLOCK}" = "FALSE" ] && [ "${CONDA_CREATE_ENV_ONLY}" = "FALSE" ] && [ "${DRYRUN}" = "FALSE" ] && [ "${DIAGRAM}" = "FALSE" ] && [ "${REPORT}" = "FALSE" ] && [ "${USE_CONDA}" = "FALSE" ] && [ "${HELP}" = "FALSE" ] ) ; then
  echo -e "\nYou did not specify any action to perform. To print the help run ./runSnakemakeWorkflow.sh --help --workflow-path PATH/TO/SMK_DIR\n"
  exit 0
fi

### Check if folder exists
if [[ ! -d "${WORKFLOW_PATH}/SCRIPTS/" ]] ; then
  echo -e "\nERROR: No SCRIPTS/ folder was found in the provided workflow path (${WORKFLOW_PATH}). Please clone or copy the whole WORKFLOW folder from GitHub: https://github.com/BioInfo-GE2POP-BLE/CAPTURE_PIPELINES_SNAKEMAKE containing all sub-directories."
  echo -e "\nExiting.\n"
  exit 1
fi


### Remove Windows \r from help file
nb_carriage_returns=$(grep -c $'\r' ${WORKFLOW_PATH}/SCRIPTS/launcher_help.txt)
if [[ "$nb_carriage_returns" -gt 0 ]] ; then
  sed -i 's/\r$//g' ${WORKFLOW_PATH}/SCRIPTS/launcher_help.txt
  sed -i 's/\r/\n/g' ${WORKFLOW_PATH}/SCRIPTS/launcher_help.txt

fi


### Print the help
if [ "${HELP}" = "TRUE" ] ; then
  cat ${WORKFLOW_PATH}/SCRIPTS/launcher_help.txt
  exit 0
fi


### Make scripts executable
for script in $(ls "${WORKFLOW_PATH}/SCRIPTS/") ; do
  if [[ -f "${WORKFLOW_PATH}/SCRIPTS/${script}" ]] ; then
    nb_carriage_returns=$(grep -c $'\r' ${WORKFLOW_PATH}/SCRIPTS/${script})
    if [[ "$nb_carriage_returns" -gt 0 ]] ; then
      echo "Removing windows carriage returns in ${script}..."
      sed -i 's/\r$//g' ${WORKFLOW_PATH}/SCRIPTS/$script
      sed -i 's/\r/\n/g' ${WORKFLOW_PATH}/SCRIPTS/$script
    fi
    if [[ ! -x "${WORKFLOW_PATH}/SCRIPTS/${script}" ]] ; then
      echo "Making $script executable..."
      chmod 755 "${WORKFLOW_PATH}/SCRIPTS/${script}"
    fi
  fi
done


### Check variables and paths
source "${WORKFLOW_PATH}/SCRIPTS/launcher_basicCheck.sh"




### RUN APPROPRIATE SNAKEMAKE COMMANDS ###

## UNLOCK ##
if [ "${UNLOCK}" = "TRUE" ] ; then
  echo -e "\nCalling Snakemake:"
  echo -e "snakemake --snakefile ${WORKFLOW_PATH}/${WORKFLOW_SMK} --unlock --configfile ${CONFIG} ${EXTRA_SNAKEMAKE_OPTIONS}\n"
  snakemake --snakefile ${WORKFLOW_PATH}/${WORKFLOW_SMK} --unlock --configfile ${CONFIG} ${EXTRA_SNAKEMAKE_OPTIONS}
fi


## CREATE CONDA ENVIRONMENT ##
if [ "${CONDA_CREATE_ENV_ONLY}" = "TRUE" ] ; then
  echo -e "\nCalling Snakemake:"
  echo -e "snakemake --snakefile ${WORKFLOW_PATH}/${WORKFLOW_SMK} $PRINTSHELLCMDS --use-conda --conda-create-envs-only --jobs $JOBS --configfile ${CONFIG} ${EXTRA_SNAKEMAKE_OPTIONS}\n"
  snakemake --snakefile ${WORKFLOW_PATH}/${WORKFLOW_SMK} $PRINTSHELLCMDS --use-conda --conda-create-envs-only --jobs $JOBS --configfile ${CONFIG} ${EXTRA_SNAKEMAKE_OPTIONS}
  exit 0
fi


## DRYRUN ##
if [ "${DRYRUN}" = "TRUE" ] ; then
  source "${WORKFLOW_PATH}/SCRIPTS/launcher_${WORKFLOW}Check.sh"
  echo -e "\nCalling Snakemake:"
  echo -e "snakemake --snakefile ${WORKFLOW_PATH}/${WORKFLOW_SMK} $PRINTSHELLCMDS --dryrun --dag --forceall --configfile ${CONFIG} ${EXTRA_SNAKEMAKE_OPTIONS}\n"
  snakemake --snakefile ${WORKFLOW_PATH}/${WORKFLOW_SMK} $PRINTSHELLCMDS --dryrun --dag --forceall --configfile ${CONFIG} ${EXTRA_SNAKEMAKE_OPTIONS}
  exit 0
fi


## DIAGRAM ##
if [ "${DIAGRAM}" = "TRUE" ] ; then
  echo -e "\nCalling Snakemake:"
  echo -e "snakemake --snakefile ${WORKFLOW_PATH}/${WORKFLOW_SMK} $PRINTSHELLCMDS --dryrun --dag --forceall --configfile ${CONFIG} ${EXTRA_SNAKEMAKE_OPTIONS} | dot -Tsvg > $DIAGRAM_NAME\n"
  snakemake --snakefile ${WORKFLOW_PATH}/${WORKFLOW_SMK} $PRINTSHELLCMDS --dryrun --dag --forceall --configfile ${CONFIG} ${EXTRA_SNAKEMAKE_OPTIONS} | dot -Tsvg > $DIAGRAM_NAME
  exit 0
fi


## REPORT ##
if [ "${REPORT}" = "TRUE" ] ; then
  echo -e "\nCalling Snakemake:"
  echo -e "snakemake --snakefile ${WORKFLOW_PATH}/${WORKFLOW_SMK} $PRINTSHELLCMDS --report $REPORT_NAME --configfile ${CONFIG} ${EXTRA_SNAKEMAKE_OPTIONS}\n"
  snakemake --snakefile ${WORKFLOW_PATH}/${WORKFLOW_SMK} $PRINTSHELLCMDS --report $REPORT_NAME --configfile ${CONFIG} ${EXTRA_SNAKEMAKE_OPTIONS}
  exit 0
fi


## RUN WITH CONDA ##
if [ "${USE_CONDA}" = "TRUE" ] ; then
  source "${WORKFLOW_PATH}/SCRIPTS/launcher_allWorkflowsCheck.sh"
  source "${WORKFLOW_PATH}/SCRIPTS/launcher_${WORKFLOW}Check.sh"
  echo -e "\nCalling Snakemake:"
  echo -e "snakemake --snakefile ${WORKFLOW_PATH}/${WORKFLOW_SMK} $PRINTSHELLCMDS $FORCEALL --latency-wait $LATENCY_WAIT --jobs $JOBS --use-conda ${CLUSTER_CONFIG_CMD} --configfile ${CONFIG} ${PROFILE} ${EXTRA_SNAKEMAKE_OPTIONS}\n"
  snakemake --snakefile ${WORKFLOW_PATH}/${WORKFLOW_SMK} $PRINTSHELLCMDS $FORCEALL --latency-wait $LATENCY_WAIT --jobs $JOBS --use-conda ${CLUSTER_CONFIG_CMD} --configfile ${CONFIG} ${PROFILE} ${EXTRA_SNAKEMAKE_OPTIONS}
  exit 0
fi


## RUN WITH CONDA + SINGULARITY ##
#if [ "${USE_CONDA_AND_SINGULARITY}" = "TRUE" ] ; then

#  echo "This does not work yet sorry :("
#  exit 0

#  source "${WORKFLOW_PATH}/SCRIPTS/launcher_allWorkflowsCheck.sh"
#  source "${WORKFLOW_PATH}/SCRIPTS/launcher_${WORKFLOW}Check.sh"
#  snakemake --snakefile ${WORKFLOW_PATH}/${WORKFLOW_SMK} $PRINTSHELLCMDS $FORCEALL $LATENCY_WAIT $JOBS --use-conda --use-singularity ${SINGULARITY_ARGS} --cluster-config ${CLUSTER_CONFIG} --configfile ${CONFIG} ${PROFILE}
#fi
