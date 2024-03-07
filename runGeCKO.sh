#!/usr/bin/env bash


### v WRITE YOUR MODULE LOAD OR CONDA ACTIVATE HERE v ###

#module purge
#conda activate GeCKO_env

### ^ WRITE YOUR MODULE LOAD OR CONDA ACTIVATE HERE ^ ###



### DEFAULT OPTIONS
JOBS=1
LATENCY_WAIT=20

### DEFAULT ACTION VALUES
HELP="FALSE"
DRYRUN="FALSE"
DIAGRAM="FALSE"
REPORT="FALSE"
CREATE_GECKO_CONDA_ENV="FALSE"
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
    --create-gecko-env)
    CREATE_GECKO_CONDA_ENV="TRUE"
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


# --------------------------------------------------------------------------------------------------------------#


### Check if folder exists
if [[ ! -d "${GeCKO_path}/launcher_files/" ]] ; then
  echo -e "\nERROR: No launcher_files/ folder was found in the provided workflow path (${GeCKO_path}). Please clone or copy the whole repository from GitHub: https://github.com/GE2POP/GeCKO containing all sub-directories."
  echo -e "\nExiting.\n"
  exit 1
fi


### Remove CR and make scripts executable
for file in $(ls "${GeCKO_path}/launcher_files/") ; do
  if [[ -f "${GeCKO_path}/launcher_files/${file}" ]] ; then
    if grep -q $'\r' ${GeCKO_path}/launcher_files/${file}; then
      echo "Removing windows carriage returns in ${file}..."
      sed -i 's/\r$//g' ${GeCKO_path}/launcher_files/$file
      sed -i 's/\r/\n/g' ${GeCKO_path}/launcher_files/$file
    fi
    if [[ ${file} == *.sh && ! -x "${GeCKO_path}/launcher_files/${file}" ]] ; then
      echo "Making $file executable..."
      chmod 755 "${GeCKO_path}/launcher_files/${file}"
    fi
  fi
done


### Print the help
if [ "${HELP}" = "TRUE" ] ; then
  cat ${GeCKO_path}/launcher_files/launcher_help.txt
  exit 0
fi


### Create GeCKO conda environment
if [ "${CREATE_GECKO_CONDA_ENV}" = "TRUE" ] ; then
  if ! command -v conda &> /dev/null ; then
    echo -e "\nERROR: Conda is not available. Please install it, or make it available to your working environment (eg: module load it)."
    exit 1
  else
    echo "Creating GeCKO conda environment..."
    conda env create -n GeCKO_env -f ${GeCKO_path}/launcher_files/GeCKO_env.yaml
    exit 0
  fi
fi


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
