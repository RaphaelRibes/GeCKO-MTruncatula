#!/bin/bash

HERE=$PWD

### 0/ Variables help
workflow_help=$(grep -e "^--workflow " ${WORKFLOW_PATH}/SCRIPTS/launcher_help.txt)
config_file_help=$(grep -e "^--config-file " ${WORKFLOW_PATH}/SCRIPTS/launcher_help.txt)


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

nb_spaces=$(grep -c '  .*' $CONFIG)
if [[ "$nb_spaces" -gt 0 ]] ; then
	echo "Removing extra spaces in ${CONFIG}..."
  sed -i 's/  */ /g' $CONFIG
fi


nb_carriage_returns=$(grep -c $'\r' $CONFIG)
if [[ "$nb_carriage_returns" -gt 0 ]] ; then
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

if [[ ! -f "${WORKFLOW_PATH}/${WORKFLOW_SMK}" ]] ; then
	echo -e "\nERROR: The workflow name provided with --workflow is incorrect."
	echo "As a reminder:"
	echo $workflow_help
	echo -e "\nExiting.\n"
	exit 1
fi

nb_carriage_returns=$(grep -c $'\r' ${WORKFLOW_PATH}/$WORKFLOW_SMK)
if [[ "$nb_carriage_returns" -gt 0 ]] ; then
	echo "Removing windows carriage returns in ${WORKFLOW_SMK}..."
	sed -i 's/\r$//g' ${WORKFLOW_PATH}/$WORKFLOW_SMK
	sed -i 's/\r/\n/g' ${WORKFLOW_PATH}/$WORKFLOW_SMK
fi

### 4/ Check if Snakemake module is available
if ! [[ -x "$(command -v snakemake)" ]] ; then
  echo -e "\nERROR: Snakemake is not available. You must install it, or make it available to your working environment (eg: module load it or activate it with conda)."
  echo "As a reminder:"
  awk '/^- Make sure Snakemake and Conda/,/^$/' ${WORKFLOW_PATH}/SCRIPTS/launcher_help.txt
  echo -e "\nExiting.\n"
  exit 1
fi
