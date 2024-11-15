#!/bin/bash


HERE=$PWD

### 0/ Variables help
workflow_help=$(grep -e "^--workflow " ${GeCKO_path}/launcher_files/launcher_help.txt)
config_file_help=$(grep -e "^--config-file " ${GeCKO_path}/launcher_files/launcher_help.txt)
cluster_profile_help=$(grep -e "^--cluster-profile " ${GeCKO_path}/launcher_files/launcher_help.txt)


### 1/ WORKFLOW FOLDER AND ITS CONTENTS ###

if [[ -z "$WORKFLOW" || "$WORKFLOW" = --* ]] ; then
	echo -e "\nERROR: the --workflow parameter is missing, please include it in your command."
	echo "As a reminder:"
	echo $workflow_help
	echo -e "\nExiting.\n"
	exit 1
fi

workflow_folder_name=$(grep -w $WORKFLOW ${GeCKO_path}/launcher_files/workflows_list.tsv | cut -f2)
if [[ -z $workflow_folder_name ]] ; then
	echo -e "\nERROR: The workflow name provided with --workflow is unknown."
	echo "As a reminder:"
	echo $workflow_help
	echo -e "\nExiting.\n"
	exit 1
fi

if [[ ! -d "${GeCKO_path}/${workflow_folder_name}/" ]] ; then
  echo -e "\nERROR: No ${workflow_folder_name}/ folder was found in the provided workflow path (${GeCKO_path}). Please clone or copy the whole repository from GitHub: https://github.com/GE2POP/GeCKO containing all sub-directories."
  echo -e "\nExiting.\n"
  exit 1
fi


workflow_path="${GeCKO_path}/${workflow_folder_name}/WORKFLOW"
if [[ ! -d "${GeCKO_path}/${workflow_folder_name}/WORKFLOW" ]] ; then
  echo -e "\nERROR: No ${workflow_folder_name}/WORKFLOW folder was found in the provided workflow path (${GeCKO_path}). Please clone or copy the whole repository from GitHub: https://github.com/GE2POP/GeCKO containing all sub-directories."
  echo -e "\nExiting.\n"
  exit 1
fi

workflow_scripts_folder="${workflow_path}/SCRIPTS"
workflow_model_files_folder="${workflow_path}/SCRIPTS/model_files"
if [[ -d ${workflow_scripts_folder} ]] ; then
	for script in $(ls "${workflow_scripts_folder}") ; do
	  if [[ -f "${workflow_scripts_folder}/${script}" ]] ; then
	    if grep -q $'\r' ${workflow_scripts_folder}/${script}; then
	      echo "Removing windows carriage returns in ${script}..."
	      sed -i 's/\r$//g' ${workflow_scripts_folder}/$script
	      sed -i 's/\r/\n/g' ${workflow_scripts_folder}/$script
	    fi
	    if [[ ! -x "${workflow_scripts_folder}/${script}" ]] ; then
	      echo "Making $script executable..."
	      chmod 755 "${workflow_scripts_folder}/${script}"
	    fi
	  fi
	done
	if [[ -d ${workflow_model_files_folder} ]] ; then
		for config_file in $(ls ${workflow_model_files_folder}) ; do
      if grep -q $'\r' ${workflow_model_files_folder}/${config_file}; then
	      		echo "Removing windows carriage returns in ${config_file}..."
	      		sed -i 's/\r$//g' ${workflow_model_files_folder}/$config_file
	      		sed -i 's/\r/\n/g' ${workflow_model_files_folder}/$config_file
	    	fi
		done
	fi
fi




### 2/ CONFIG_FILE ###

CONFIG=$(absolutePath $CONFIG)

# In case it was not provided
if [[ -z "$CONFIG" || "$CONFIG" = --* || "$CONFIG" = -* ]] ; then
  if [[ -f "CONFIG/config_${WORKFLOW}.yml" ]] ; then
    CONFIG=CONFIG/config_${WORKFLOW}.yml
		CONFIG=$(absolutePath $CONFIG)
    echo -e "\nINFO: The config_${WORKFLOW}.yml file was automatically found in ${HERE}/CONFIG and will be used as the config file for this workflow."
    echo -e "If you would rather use another config file, please provide it with --config-file\n"
  else
    echo -e "\nERROR: The workflow config file cannot be found. Please provide your config file with --config-file or place it in a CONFIG folder under the name config_${WORKFLOW}.yml."
		echo "As a reminder:"
		echo $config_file_help
		echo -e "\nExiting.\n"
    exit 1
  fi # In case it was provided
elif [[ ! -f "$CONFIG" ]] ; then
    echo -e "\nERROR: the file given in the --config-file parameter is not valid. Please make sure the file exists and the path is correctly written."
		echo "As a reminder:"
		echo $config_file_help
		echo -e "\nExiting.\n"
    exit 1
fi

if grep -q '  .*' $CONFIG; then
	echo "Removing extra spaces in ${CONFIG}..."
  sed -i 's/  */ /g' $CONFIG
fi


if grep -q $'\r' $CONFIG; then
	echo "Removing windows carriage returns in ${CONFIG}..."
  sed -i 's/\r$//g' $CONFIG
  sed -i 's/\r/\n/g' $CONFIG
fi


# Format
nb_col_config_file=$(grep -v '^#' $CONFIG | sed 's/#.*$//' | grep -v '^[[:space:]]*$' | awk -F ': ' '{print NF}' | sort | uniq)
nb_tabs_config_file=$(grep -v '^#' $CONFIG | sed 's/#.*$//' | grep -v '^[[:space:]]*$' | grep -c $'\t')
if [[ "$nb_col_config_file" != 2 || $nb_tabs_config_file -gt 0 ]] ; then
  echo -e "\nERROR: Your config file (${CONFIG}) does not seem to be properly formated."
  echo "As a reminder:"
  echo "Variables must be specified in the yaml format, with one variable per row, followed by a ':' and a space before the assigned value, for example : VARIABLE_NAME: value."
  echo "Some variables can be left empty, for example: VARIABLE_NAME: \"\" "
  echo -e "\nExiting.\n"
  exit 1
fi


### 3/ WORKFLOW_SMK
WORKFLOW_SMK="${WORKFLOW}.smk"

# if DataCleaning: is it paired or single
if [[ "$WORKFLOW" = "DataCleaning" ]] ; then
	PAIRED_END=$(grep "^PAIRED_END" $CONFIG | cut -f2 -d ' ')
	if [[ "$PAIRED_END" == "TRUE" || "$PAIRED_END" == "True" || "$PAIRED_END" == "true" || "$PAIRED_END" == "T" ]] ; then
		WORKFLOW_SMK="${WORKFLOW}_PairedEnd.smk"
	elif [[ "$PAIRED_END" == "FALSE" || "$PAIRED_END" == "False" || "$PAIRED_END" == "false" || "$PAIRED_END" == "F" ]] ; then
		WORKFLOW_SMK="${WORKFLOW}_SingleEnd.smk"
	else
		echo -e "\nERROR: The PAIRED_END variable is either missing or incorrect in your config file (${CONFIG}). Please set it to TRUE or FALSE."
		echo "As a reminder:"
		echo "PAIRED_END: set to TRUE in case of paired end data (R1 + R2), to FALSE in case of single end data."
		echo -e "\nExiting.\n"
		exit 1
	fi
fi

if [[ ! -f "${workflow_path}/${WORKFLOW_SMK}" ]] ; then
	echo -e "\nERROR: The workflow name provided with --workflow is incorrect."
	echo "As a reminder:"
	echo $workflow_help
	echo -e "\nExiting.\n"
	exit 1
fi


if grep -q $'\r' ${workflow_path}/$WORKFLOW_SMK; then
	echo "Removing windows carriage returns in ${WORKFLOW_SMK}..."
	sed -i 's/\r$//g' ${workflow_path}/$WORKFLOW_SMK
	sed -i 's/\r/\n/g' ${workflow_path}/$WORKFLOW_SMK
fi

### 4/ Check if Snakemake module is available
if ! command -v snakemake &> /dev/null; then
  echo -e "\nERROR: Snakemake is not available. You must install it, or make it available to your working environment (eg: module load it or activate it with conda)."
  echo -e "\nExiting.\n"
  exit 1
fi



### 5/ CLUSTER_PROFILE and jobs


if [[ ! -z "$CLUSTER_PROFILE" ]] ; then
	if [[ "$JOBS" = 1 ]] ; then
	  echo -e "\nWARNING: The workflow will be run with only 1 job allowed at a time. The tasks will not be parallelized. You can increase the maximum number of jobs allowed to be run in parallel with the --jobs option.\n"
	fi
  if [[ -d "$CLUSTER_PROFILE" ]] ; then
    if [[ -f "${CLUSTER_PROFILE}/config.yaml" ]] ; then
      if grep -q $'\r' ${CLUSTER_PROFILE}/config.yaml; then
        echo "Removing windows carriage returns in ${CLUSTER_PROFILE}/config.yaml..."
        sed -i 's/\r$//g' ${CLUSTER_PROFILE}/config.yaml
        sed -i 's/\r/\n/g' ${CLUSTER_PROFILE}/config.yaml
      fi
      CLUSTER_PROFILE=$(absolutePath $CLUSTER_PROFILE)
      PROFILE="--profile $CLUSTER_PROFILE"
      PROFILE_FILE=${CLUSTER_PROFILE}/config.yaml
    else
      echo -e "\nERROR: $CLUSTER_PROFILE must contain a config.yaml file. Please make sure it exists and it is properly named."
		  echo "As a reminder:"
		  echo $cluster_profile_help
		  echo -e "\nExiting.\n"
      exit 1
    fi
  else
    echo -e "\nERROR: the path given in the --cluster-profile parameter is not valid. Please make sure the folder exists and the path is correctly written."
		echo "As a reminder:"
		echo $cluster_profile_help
		echo -e "\nExiting.\n"
    exit 1
  fi
else
  PROFILE=""
  PROFILE_FILE="NULL"
  echo -e "\nWARNING: No cluster profile was provided. The pipeline will be run without submitting any job and any parallelization.\n"
fi


### 6/ Check if Singularity is available
if ! command -v singularity &> /dev/null; then
  echo -e "\nERROR: Singularity is not available. You must install it, or make it available to your working environment (eg: module load it)."
  echo -e "\nExiting.\n"
  exit 1
fi


