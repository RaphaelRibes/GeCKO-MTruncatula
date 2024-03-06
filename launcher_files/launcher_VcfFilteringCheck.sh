#!/bin/bash

# Model config file
model_config="${workflow_path}/SCRIPTS/model_files/config_VcfFiltering.yml"


# Config file variables
given_variables=$(grep -v '^#' $CONFIG | grep -v '^$' | cut -f1 -d ' ' || true)
expected_variables=$(grep -v '^#' $model_config | grep -v '^$' | cut -f1 -d ' ' || true)

for var in $expected_variables ; do
  var_in_config=$(grep "^$var " $CONFIG || true)
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


## INPUT FILES
VCF_FILE=$(grep "^VCF_FILE:" $CONFIG | sed 's/#.*$//' | cut -d ' ' -f2 | sed 's/"//g' || true)

if [[ -z "$VCF_FILE" ]] ; then
  echo -e "\nERROR: The VCF_FILE variable is missing in your config file (${CONFIG}). Please provide the vcf file produced by your variant caller."
  echo -e "\nExiting.\n"
  exit 1
elif [[ ! -f "$VCF_FILE" ]] ; then
  echo -e "\nERROR: The VCF_FILE file (${VCF_FILE}) provided in your config file (${CONFIG}) does not exist. Please make sure the file exists and the path is correctly written."
  echo -e "\nExiting.\n"
  exit 1
fi
if [[ "$VCF_FILE" != *.vcf.gz ]] ; then
  echo -e "\nERROR: The VCF_FILE file (${VCF_FILE}) provided in your config file (${CONFIG}) does not have a proper vcf extension. Please make sure the file is a zipped vcf file and ends with '.vcf.gz'."
  echo -e "\nExiting.\n"
  exit 1
fi


### VCF FILTERING PARAMETERS ###
MAX_NA_PER_SAMPLE=$(grep "^MAX_NA_PER_SAMPLE:" $CONFIG | sed 's/#.*$//' | cut -d ' ' -f2 | sed 's/"//g' || true)

if [[ -z "$MAX_NA_PER_SAMPLE" ]] ; then
  echo -e "\nERROR: You must provide the MAX_NA_PER_SAMPLE in the config_file."
  echo -e "\nExiting.\n"
  exit 1
fi
