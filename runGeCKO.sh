#!/bin/bash


### v WRITE YOUR MODULE LOADS HERE v ###

# module purge
# module load snakemake/5.13.0
# module load anaconda/python3.8


### ^ WRITE YOUR MODULE LOADS HERE ^ ###



### DEFAULT OPTIONS
JOBS=1
LATENCY_WAIT=20

### DEFAULT ACTION VALUES
HELP="FALSE"
DRYRUN="FALSE"
DIAGRAM="FALSE"
REPORT="FALSE"
USE_CONDA="TRUE"
CONDA_ENV_PATH=""
CONDA_ENV_PATH_CMD=""

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
    --jobs)
    JOBS="$2"
    shift
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
    --conda-env-path)
    CONDA_ENV_PATH="$2"
    shift
    shift
    ;;
    --extra-snakemake-options)
    EXTRA_SNAKEMAKE_OPTIONS="$2"
    shift
    shift
    ;;
    -*)
    echo -e "\nWARNING: $1 option is unknown and will be ignored.\n"
    POSITIONAL+=("$1")
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

if [[ -z "$WORKFLOW_PATH" || "$WORKFLOW_PATH" = --* || "$WORKFLOW_PATH" = -* ]] ; then
	echo -e "\nERROR: the --workflow-path parameter is missing, please include it in your command."
	echo "As a reminder:"
	echo "--workflow-path [...]: the path to the directory you cloned from GitHub, ending with /GeCKO"
	echo -e "\nExiting.\n"
	exit 1
fi
if [ ! -d "$WORKFLOW_PATH" ] ; then
	echo -e "\nERROR: the path given in the --workflow-path parameter is not valid. Please make sure the directory exists and the path is correctly written."
	echo "As a reminder:"
  echo "--workflow-path [...]: the path to the directory you cloned from GitHub, ending with /GeCKO"
	echo -e "\nExiting.\n"
	exit 1
fi

WORKFLOW_PATH=$(absolutePath $WORKFLOW_PATH)

# --------------------------------------------------------------------------------------------------------------#


### Check if folder exists
if [[ ! -d "${WORKFLOW_PATH}/scripts/" ]] ; then
  echo -e "\nERROR: No scripts/ folder was found in the provided workflow path (${WORKFLOW_PATH}). Please clone or copy the whole repository from GitHub: https://github.com/GE2POP/GeCKO containing all sub-directories."
  echo -e "\nExiting.\n"
  exit 1
fi


### Remove Windows \r from help file
nb_carriage_returns=$(grep -c $'\r' ${WORKFLOW_PATH}/scripts/launcher_help.txt)
if [[ "$nb_carriage_returns" -gt 0 ]] ; then
  sed -i 's/\r$//g' ${WORKFLOW_PATH}/scripts/launcher_help.txt
  sed -i 's/\r/\n/g' ${WORKFLOW_PATH}/scripts/launcher_help.txt

fi


### Print the help
if [ "${HELP}" = "TRUE" ] ; then
  cat ${WORKFLOW_PATH}/scripts/launcher_help.txt
  exit 0
fi


### Make scripts executable
for script in $(ls "${WORKFLOW_PATH}/scripts/") ; do
  if [[ -f "${WORKFLOW_PATH}/scripts/${script}" ]] ; then
    nb_carriage_returns=$(grep -c $'\r' ${WORKFLOW_PATH}/scripts/${script})
    if [[ "$nb_carriage_returns" -gt 0 ]] ; then
      echo "Removing windows carriage returns in ${script}..."
      sed -i 's/\r$//g' ${WORKFLOW_PATH}/scripts/$script
      sed -i 's/\r/\n/g' ${WORKFLOW_PATH}/scripts/$script
    fi
    if [[ ! -x "${WORKFLOW_PATH}/scripts/${script}" ]] ; then
      echo "Making $script executable..."
      chmod 755 "${WORKFLOW_PATH}/scripts/${script}"
    fi
  fi
done


### Check variables and paths
source "${WORKFLOW_PATH}/scripts/launcher_allWorkflowsCheck.sh"
if [[ -f "${WORKFLOW_PATH}/scripts/launcher_${WORKFLOW}Check.sh" ]] ; then
  source "${WORKFLOW_PATH}/scripts/launcher_${WORKFLOW}Check.sh"
fi


### RUN APPROPRIATE SNAKEMAKE COMMANDS ###
# Always unlock in case the folder is locked
snakemake --snakefile ${workflow_folder}/${WORKFLOW_SMK} --jobs $JOBS --unlock --configfile ${CONFIG}


## DRYRUN ##
if [ "${DRYRUN}" = "TRUE" ] ; then
  snakemake_command="snakemake --snakefile ${workflow_folder}/${WORKFLOW_SMK} --printshellcmds --dryrun --dag --forceall --configfile ${CONFIG} ${EXTRA_SNAKEMAKE_OPTIONS}"
  echo -e "\nCalling Snakemake:"
  echo -e $snakemake_command"\n"
  $snakemake_command
  exit 0
fi


## DIAGRAM ##
if [ "${DIAGRAM}" = "TRUE" ] ; then
  snakemake_command="snakemake --snakefile ${workflow_folder}/${WORKFLOW_SMK} --printshellcmds --dryrun --dag --forceall --configfile ${CONFIG} ${EXTRA_SNAKEMAKE_OPTIONS} | dot -Tsvg > $DIAGRAM_NAME"
  echo -e "\nCalling Snakemake:"
  echo -e $snakemake_command"\n"
  $snakemake_command
  exit 0
fi


## REPORT ##
if [ "${REPORT}" = "TRUE" ] ; then
  snakemake_command="snakemake --snakefile ${workflow_folder}/${WORKFLOW_SMK} --printshellcmds --report $REPORT_NAME --configfile ${CONFIG} ${EXTRA_SNAKEMAKE_OPTIONS}"
  echo -e "\nCalling Snakemake:"
  echo -e $snakemake_command"\n"
  $snakemake_command
  exit 0
fi


## RUN WITH CONDA ##
if [ "${USE_CONDA}" = "TRUE" ] ; then
  mkdir -p Logs_${WORKFLOW}Workflow

  snakemake_command="snakemake --snakefile ${workflow_folder}/${WORKFLOW_SMK} --printshellcmds $FORCEALL --latency-wait $LATENCY_WAIT --jobs $JOBS --use-conda ${CLUSTER_CONFIG_CMD} --configfile ${CONFIG} ${PROFILE} ${CONDA_ENV_PATH_CMD} ${EXTRA_SNAKEMAKE_OPTIONS}"
  echo -e "\nCalling Snakemake:"
  echo -e $snakemake_command"\n"
  $snakemake_command
  exit 0
fi
