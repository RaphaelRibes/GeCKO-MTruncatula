#!/usr/bin/env bash
set -eo pipefail

SingularityImageVersion=1.2.0

### v WRITE YOUR MODULE LOAD OR CONDA ACTIVATE HERE v ###

#module purge
#module load snakemake/7.32.4-conda
#module load singularity/3.6.3

### ^ WRITE YOUR MODULE LOAD OR CONDA ACTIVATE HERE ^ ###



### DEFAULT OPTIONS
JOBS=1
LATENCY_WAIT=20

### DEFAULT ACTION VALUES
HELP="FALSE"
DRYRUN="FALSE"
DIAGRAM="FALSE"
REPORT="FALSE"

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

printAbsolutePath () {
  if [[ ! -z "$1" ]] ; then
    fileOrFolder_printAbsolutePath=$(readlink -f "$1") ;
    echo ${fileOrFolder_printAbsolutePath}
  else
    echo $1
  fi
}

# --------------------------------------------------------------------------------------------------------------#


### Paths

GeCKO_path=$(dirname $(printAbsolutePath "$0"))
checks_path="${GeCKO_path}/launcher_files"


# --------------------------------------------------------------------------------------------------------------#

### Functions

source "${checks_path}/lib.sh"


# --------------------------------------------------------------------------------------------------------------#


### Check if launcher_files folder exists
checkMissingDir $checks_path "GeCKOdir"


### Check if Singularity/Apptainer and Snakemake are available
isAvailable "Snakemake" "snakemake"
isAvailable "Singularity/Apptainer" "singularity"


### Download the singularity container if it can't be found
dlImageSylabs "library://ge2pop_gecko/gecko/gecko:${SingularityImageVersion}" "${checks_path}/singularity_image/GeCKO.sif"


### Make scripts executable
for script in $(ls ${checks_path}/*.sh) ; do
  makeExecutable $script
done



### Print the help
if [ "${HELP}" = "TRUE" ] ; then
  cat ${checks_path}/launcher_help.txt
  exit 0
fi


### Check variables and paths
source "${checks_path}/launcher_allWorkflowsCheck.sh"
if [[ -f "${checks_path}/launcher_${WORKFLOW}Check.sh" ]] ; then
  source "${checks_path}/launcher_${WORKFLOW}Check.sh"
fi



### Unlock in case the folder is locked
snakemake --snakefile ${workflow_path}/${WORKFLOW_SMK} --jobs $JOBS --unlock --configfile ${CONFIG}


### RUN APPROPRIATE SNAKEMAKE COMMANDS ###

## DRYRUN ##
if [ "${DRYRUN}" = "TRUE" ] ; then
  snakemake_command="snakemake --snakefile ${workflow_path}/${WORKFLOW_SMK} --printshellcmds --dryrun --dag --forceall --configfile ${CONFIG} ${EXTRA_SNAKEMAKE_OPTIONS}"
  echo -e "\nCalling Snakemake:"
  echo -e $snakemake_command"\n"
  eval $snakemake_command
  exit 0
fi


## DIAGRAM ##
if [ "${DIAGRAM}" = "TRUE" ] ; then
  snakemake_command="snakemake --snakefile ${workflow_path}/${WORKFLOW_SMK} --printshellcmds --dryrun --dag --forceall --configfile ${CONFIG} ${EXTRA_SNAKEMAKE_OPTIONS} | dot -Tsvg > $DIAGRAM_NAME"
  echo -e "\nCalling Snakemake:"
  echo -e $snakemake_command"\n"
  eval $snakemake_command
  exit 0
fi


## REPORT ##
if [ "${REPORT}" = "TRUE" ] ; then
  snakemake_command="snakemake --snakefile ${workflow_path}/${WORKFLOW_SMK} --printshellcmds --report $REPORT_NAME --configfile ${CONFIG} ${EXTRA_SNAKEMAKE_OPTIONS}"
  echo -e "\nCalling Snakemake:"
  echo -e $snakemake_command"\n"
  eval $snakemake_command
  exit 0
fi


## RUN ##
snakemake_command="snakemake --snakefile ${workflow_path}/${WORKFLOW_SMK} --printshellcmds $FORCEALL --latency-wait $LATENCY_WAIT --jobs $JOBS --use-singularity --configfile ${CONFIG} ${PROFILE} --config configfile_name=${CONFIG} clusterprofile_name=${PROFILE_FILE} ${EXTRA_SNAKEMAKE_OPTIONS} --singularity-args \"--bind ${GeCKO_path} --bind $(pwd)\""
echo -e "\nCalling Snakemake:"
echo -e $snakemake_command"\n"
eval $snakemake_command
exit 0
