#!/bin/bash

# Model config file
model_config="${workflow_folder}/SCRIPTS/model_files/config_VariantsCalling.yml"


# Config file variables
given_variables=$(grep -v '^#' $CONFIG | grep -v '^$' | cut -f1 -d ' ')
expected_variables=$(grep -v '^#' $model_config | grep -v '^$' | cut -f1 -d ' ')

for var in $expected_variables ; do
  var_in_config=$(grep "^$var " $CONFIG)
  if [[ -z "$var_in_config" ]] ; then
    echo -e "\nERROR: The expected variable $var was not found in your config file (${CONFIG}). Please make sure to include it."
    echo -e "\nList of expected variables (some can be left empty but must appear in the file nonetheless) :\n${expected_variables}"
    echo -e "\nAs a reminder:"
    echo "All expected variables must be specified in YAML format, with one variable per row followed by a ':' and a space before its assigned value. For example: VARIABLE_NAME: value."
    echo "To leave a variable empty, use quotes: VARIABLE_NAME: \"\" "
    echo -e "\nExiting.\n"
    exit 1
  fi
done


                          ### CONFIG FILE VARIABLES VALUES ###


## Input data (BAMS_LIST, REFERENCE) ##
BAMS_LIST=$(grep "^BAMS_LIST:" $CONFIG | sed 's/#.*$//' | cut -d ' ' -f2 | sed 's/"//g')
REFERENCE=$(grep "^REFERENCE:" $CONFIG | sed 's/#.*$//' | cut -d ' ' -f2 | sed 's/"//g')

# **BAMS_LIST**
if [[ -z "$BAMS_LIST" ]] ; then
  echo -e "\nERROR: The BAMS_LIST variable is missing in your config file (${CONFIG}). Please provide it."
  echo "As a reminder:"
  echo "BAMS_LIST: The path to the directory containing the mapped bam files and the index file in .bam.bai format."
  echo -e "\nExiting.\n"
  exit 1
fi
if [[ ! -f "$BAMS_LIST" ]] ; then
  echo -e "\nERROR: The BAMS_LIST provided in the config file (${BARCODE_FILE}) does not exist. Please make sure the file exists and its path is correct."
  echo -e "\nExiting.\n"
  exit 1
fi

# **REFERENCE**
if [[ -z "$REFERENCE" ]] ; then
  echo -e "\nERROR: The REFERENCE variable is missing in your config file (${CONFIG}). Please provide the reference file used to map your reads."
  echo -e "\nExiting.\n"
  exit 1
elif [[ ! -f "$REFERENCE" ]] ; then
  echo -e "\nERROR: The REFERENCE file (${REFERENCE}) provided in your config file (${CONFIG}) does not exist. Please make sure the file exists and the path is correctly written."
  echo -e "\nExiting.\n"
  exit 1
fi
if [[ "$REFERENCE" != *.fasta && "$REFERENCE" != *.fas && "$REFERENCE" != *.fa ]] ; then
  echo -e "\nERROR: The REFERENCE file (${REFERENCE}) provided in your config file (${CONFIG}) does not have a proper fasta extension. Please make sure the file is a fasta file and ends with '.fa', '.fas', or '.fasta'."
  echo -e "\nExiting.\n"
  exit 1
fi


# **Variant Calling parameters**
GATK_HAPLOTYPE_CALLER_CPUS_PER_TASK=$(grep "^GATK_HAPLOTYPE_CALLER_CPUS_PER_TASK:" $CONFIG | sed 's/#.*$//' | cut -d ' ' -f2 | sed 's/"//g')
GATK_GENOMICS_DB_IMPORT_CPUS_PER_TASK=$(grep "^GATK_GENOMICS_DB_IMPORT_CPUS_PER_TASK:" $CONFIG | sed 's/#.*$//' | cut -d ' ' -f2 | sed 's/"//g')

if [[ -z "$GATK_HAPLOTYPE_CALLER_CPUS_PER_TASK" ]] ; then
  echo -e "\nERROR: You must provide the GATK_HAPLOTYPE_CALLER_CPUS_PER_TASK in the config_file."
  echo -e "\nExiting.\n"
  exit 1
fi

if [[ -z "$GATK_GENOMICS_DB_IMPORT_CPUS_PER_TASK" ]] ; then
  echo -e "\nERROR: You must provide the GATK_GENOMICS_DB_IMPORT_CPUS_PER_TASK in the config_file."
  echo -e "\nExiting.\n"
  exit 1
fi
