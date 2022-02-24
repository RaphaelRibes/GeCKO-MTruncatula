#!/bin/bash

HERE=$PWD

### 0/ Variables help
workflow_path_help=""
workflow_help=""
cluster_config_help=""
config_file_help=""



### 1/ WORKFLOW  ###

if [[ -z "$WORKFLOW" || "$WORKFLOW" = --* || "$WORKFLOW_PATH" = -* ]] ; then
	echo -e "\nERROR: the --workflow parameter is missing, please include it in your command."
	echo "As a reminder:"
	echo $workflow_help
	echo -e "\nExiting.\n"
	exit 1
fi

available_workflows=$(ls $WORKFLOW_PATH/*smk | xargs -n 1 basename | sed 's/.smk//' | sed 's/_PairedEnd//' | sed 's/_SingleEnd//' | sort | uniq)
if [[ -z $(echo $available_workflows | grep -w $WORKFLOW) ]] ; then
	echo -e "\nERROR: The workflow name provided with --workflow is incorrect."
	echo "As a reminder:"
	echo $workflow_help
	echo -e "\nExiting.\n"
	exit 1
fi


### 2/ CLUSTER_CONFIG and CONFIG_FILE ###

## CLUSTER_CONFIG ##
# In case it was not provided
if [[ -z "$CLUSTER_CONFIG" || "$CLUSTER_CONFIG" = --* || "$CLUSTER_CONFIG" = -* ]] ; then
  if [[ -f "CONFIG/cluster_config_${WORKFLOW}.json" ]] ; then
    CLUSTER_CONFIG=CONFIG/cluster_config_${WORKFLOW}.json
    echo -e "\nThe cluster_config_${WORKFLOW}.json file was automatically found in ${HERE}/CONFIG and will be used as the cluster_config file for this workflow."
    echo -e "If you would rather use another cluster_config file, please provide it with --cluster-config\n"
  else
    echo -e "\nERROR: The cluster config file cannot be found. Please provide your cluster config file with --cluster-config or place it in a CONFIG folder under the name cluster_config_${WORKFLOW}.json."
		echo "As a reminder:"
		echo $cluster_config_help
		echo -e "\nExiting.\n"
    exit 1
  fi
else # In case it was provided
  if [[ ! -f "$CLUSTER_CONFIG" ]] ; then
    echo -e "\nERROR: the file given in the --cluster-config parameter is not valid. Please make sure the file exists and the path is correctly written."
		echo "As a reminder:"
		echo $cluster_config_help
		echo -e "\nExiting.\n"
    exit 1
  fi
fi
# Absolute path
if [[ ! "$CLUSTER_CONFIG" = /* ]] ; then
  CLUSTER_CONFIG=$(readlink -f $CLUSTER_CONFIG) ;
fi
# Remove potential \r from file
sed -i 's/\r$//g' $CLUSTER_CONFIG
sed -i 's/\r/\n/g' $CLUSTER_CONFIG

## CONFIG ##
# In case it was not provided
if [[ -z "$CONFIG" || "$CONFIG" = --* || "$CONFIG" = -* ]] ; then
  if [[ -f "CONFIG/config_${WORKFLOW}.yml" ]] ; then
    CONFIG=CONFIG/config_${WORKFLOW}.yml
    echo -e "\nThe config_${WORKFLOW}.yml file was automatically found in ${HERE}/CONFIG and will be used as the config file for this workflow."
    echo -e "If you would rather use another config file, please provide it with --config-file\n"
  else
    echo -e "\nERROR: The workflow config file cannot be found. Please provide your config file with --config-file or place it in a CONFIG folder under the name config_${WORKFLOW}.yml."
		echo "As a reminder:"
		echo $config_file_help
		echo -e "\nExiting.\n"
    exit 1
  fi
else # In case it was provided
  if [[ ! -f "$CONFIG" ]] ; then
    echo -e "\nERROR: the file given in the --config-file parameter is not valid. Please make sure the file exists and the path is correctly written."
		echo "As a reminder:"
		echo $config_file_help
		echo -e "\nExiting.\n"
    exit 1
  fi
fi
# Absolute path
if [[ ! "$CONFIG" = /* ]] ; then
  CONFIG=$(readlink -f $CONFIG) ;
fi
# Remove potential \r from file
sed -i 's/\r$//g' $CONFIG
sed -i 's/\r/\n/g' $CONFIG





### 3/ WORKFLOW_SMK
WORKFLOW_SMK="${WORKFLOW}.smk"

# if DataCleaning: is it paired or single
if [[ "$WORKFLOW" = "DataCleaning" ]] ; then
	PAIRED_END=$(grep "^PAIRED_END" $CONFIG | cut -f2 -d ' ')
	if [[ "$PAIRED_END" = "TRUE" || "$PAIRED_END" == "True" || "$PAIRED_END" == "true" || "$PAIRED_END" == "T" ]] ; then
		WORKFLOW_SMK="${WORKFLOW}_PairedEnd.smk"
	elif [[ "$PAIRED_END" = "FALSE" || "$PAIRED_END" == "False" || "$PAIRED_END" == "false" || "$PAIRED_END" == "F" ]] ; then
		WORKFLOW_SMK="${WORKFLOW}_SingleEnd.smk"
	else
		echo -e "\nERROR: The PAIRED_END variable is either missing or incorrect in your config file (${CONFIG}). Please set it to TRUE or FALSE."
		echo "As a reminder:"
		echo "PAIRED_END: set to TRUE in case of paired end data (R1 + R2), to FALSE in case of single end data."
		echo -e "\nExiting.\n"
		exit 1
	fi
fi

# if the smk file doesn't exist

if [[ ! -f "${WORKFLOW_PATH}/${WORKFLOW_SMK}" ]] ; then
	echo -e "\nERROR: The workflow name provided with --workflow is incorrect."
	echo "As a reminder:"
	echo $workflow_help
	echo -e "\nExiting.\n"
	exit 1
fi
