#!/usr/bin/env bash

set -e -o pipefail


### v WRITE YOUR MODULE LOADS HERE v ###

module purge
source /home/girodollej/.bashrc_new
conda activate snakemake7.32.4_mamba_python
#module load snakemake/7.15.1-conda
#module load snakemake/8.4.2-conda
#conda activate snakemake7.15.2_mamba_python
#conda activate snakemake8.4.2_mamba
#conda activate snakemake8.4.8_mamba
#conda activate snakemake8.5.3_mamba
#module load anaconda/python3.8


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
    #--workflow-path)
    #WORKFLOW_PATH="$2"
    #shift
    #shift
    #;;
    --cluster-profile)
    CLUSTER_PROFILE="$2"
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
    fileOrFolder_absolutePath=$(readlink -f "$1") ;
    echo ${fileOrFolder_absolutePath%/}
  else
    echo ${1%/}
  fi
}

# --------------------------------------------------------------------------------------------------------------#

### GeCKO_path (needed for next steps)
GeCKO_path=$(dirname $(absolutePath "$0"))

#if [[ -z "$WORKFLOW_PATH" || "$WORKFLOW_PATH" = --* || "$WORKFLOW_PATH" = -* ]] ; then
#	echo -e "\nERROR: the --workflow-path parameter is missing, please include it in your command."
#	echo "As a reminder:"
#	echo "--workflow-path [...]: the path to the directory you cloned from GitHub, ending with /GeCKO"
#	echo -e "\nExiting.\n"
#	exit 1
#fi
#if [ ! -d "$WORKFLOW_PATH" ] ; then
#	echo -e "\nERROR: the path given in the --workflow-path parameter is not valid. Please make sure the directory exists and the path is correctly written."
#	echo "As a reminder:"
#  echo "--workflow-path [...]: the path to the directory you cloned from GitHub, ending with /GeCKO"
#	echo -e "\nExiting.\n"
#	exit 1
#fi

#WORKFLOW_PATH=$(absolutePath $WORKFLOW_PATH)



# --------------------------------------------------------------------------------------------------------------#


### Check if folder exists
if [[ ! -d "${GeCKO_path}/launcher_files/" ]] ; then
  echo -e "\nERROR: No launcher_files/ folder was found in the provided workflow path (${GeCKO_path}). Please clone or copy the whole repository from GitHub: https://github.com/GE2POP/GeCKO containing all sub-directories."
  echo -e "\nExiting.\n"
  exit 1
fi


### Remove Windows \r from help file

if grep -q $'\r' ${GeCKO_path}/launcher_files/launcher_help.txt; then
  sed -i 's/\r$//g' ${GeCKO_path}/launcher_files/launcher_help.txt
  sed -i 's/\r/\n/g' ${GeCKO_path}/launcher_files/launcher_help.txt
fi



### Print the help
if [ "${HELP}" = "TRUE" ] ; then
  cat ${GeCKO_path}/launcher_files/launcher_help.txt
  exit 0
fi


### Make scripts executable
for script in $(ls "${GeCKO_path}/launcher_files/") ; do
  if [[ -f "${GeCKO_path}/launcher_files/${script}" ]] ; then
    if grep -q $'\r' ${GeCKO_path}/launcher_files/${script}; then
      echo "Removing windows carriage returns in ${script}..."
      sed -i 's/\r$//g' ${GeCKO_path}/launcher_files/$script
      sed -i 's/\r/\n/g' ${GeCKO_path}/launcher_files/$script
    fi
    if [[ ! -x "${GeCKO_path}/launcher_files/${script}" ]] ; then
      echo "Making $script executable..."
      chmod 755 "${GeCKO_path}/launcher_files/${script}"
    fi
  fi
done


### Check variables and paths
source "${GeCKO_path}/launcher_files/launcher_allWorkflowsCheck.sh"
if [[ -f "${GeCKO_path}/launcher_files/launcher_${WORKFLOW}Check.sh" ]] ; then
  source "${GeCKO_path}/launcher_files/launcher_${WORKFLOW}Check.sh"
fi


### RUN APPROPRIATE SNAKEMAKE COMMANDS ###
# Always unlock in case the folder is locked
snakemake --snakefile ${workflow_path}/${WORKFLOW_SMK} --jobs $JOBS --unlock --configfile ${CONFIG}


## DRYRUN ##
if [ "${DRYRUN}" = "TRUE" ] ; then
  snakemake_command="snakemake --snakefile ${workflow_path}/${WORKFLOW_SMK} --printshellcmds --dryrun --dag --forceall --configfile ${CONFIG} ${EXTRA_SNAKEMAKE_OPTIONS}"
  echo -e "\nCalling Snakemake:"
  echo -e $snakemake_command"\n"
  $snakemake_command
  exit 0
fi


## DIAGRAM ##
if [ "${DIAGRAM}" = "TRUE" ] ; then
  snakemake_command="snakemake --snakefile ${workflow_path}/${WORKFLOW_SMK} --printshellcmds --dryrun --dag --forceall --configfile ${CONFIG} ${EXTRA_SNAKEMAKE_OPTIONS} | dot -Tsvg > $DIAGRAM_NAME"
  echo -e "\nCalling Snakemake:"
  echo -e $snakemake_command"\n"
  $snakemake_command
  exit 0
fi


## REPORT ##
if [ "${REPORT}" = "TRUE" ] ; then
  snakemake_command="snakemake --snakefile ${workflow_path}/${WORKFLOW_SMK} --printshellcmds --report $REPORT_NAME --configfile ${CONFIG} ${EXTRA_SNAKEMAKE_OPTIONS}"
  echo -e "\nCalling Snakemake:"
  echo -e $snakemake_command"\n"
  $snakemake_command
  exit 0
fi


## RUN WITH CONDA ##
if [ "${USE_CONDA}" = "TRUE" ] ; then
  snakemake_command="snakemake --snakefile ${workflow_path}/${WORKFLOW_SMK} --printshellcmds $FORCEALL --latency-wait $LATENCY_WAIT --jobs $JOBS --use-conda --configfile ${CONFIG} ${PROFILE} --config configfile_name=${CONFIG} clusterprofile_name=${PROFILE_FILE} ${CONDA_ENV_PATH_CMD} ${EXTRA_SNAKEMAKE_OPTIONS}"
  echo -e "\nCalling Snakemake:"
  echo -e $snakemake_command"\n"
  $snakemake_command
  exit 0
fi
