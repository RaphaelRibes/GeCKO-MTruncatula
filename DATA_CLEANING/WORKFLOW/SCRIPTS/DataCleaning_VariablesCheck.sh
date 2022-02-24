#!/bin/bash


# Model config file
if [[ "$WORKFLOW_SMK" = "${WORKFLOW}_PairedEnd.smk" ]] ; then
  echo -e "*** Paired end data is expected (PAIRED_END set to TRUE) ***\n"
  model_config="${WORKFLOW_PATH}/SCRIPTS/model_files/config_DataCleaning_PE.yml"
else
  echo -e "*** Single end data is expected (PAIRED_END set to FALSE) ***\n"
  model_config="${WORKFLOW_PATH}/SCRIPTS/model_files/config_DataCleaning_SE.yml"
fi

# Config_file format
sed -i 's/  */ /g' $CONFIG
nb_col_config_file=$(grep -v '#' $CONFIG | grep -v '^[[:space:]]*$' | awk -F ': ' '{print NF}' | sort | uniq)
if [[ $nb_col_config_file != 2 ]] ; then
  echo -e "\nERROR: Your config file (${CONFIG}) does not seem to be properly formated."
  echo "As a reminder:"
  echo "Variables must be specified in the yaml format, with one variable per row, followed by a ':' and a space before the assigned value, for example : VARIABLE_NAME: value."
  echo "Some variables can be left empty, for example: VARIABLE_NAME: \"\" "
  echo -e "\nExiting.\n"
  exit 1
fi


# Config file variables
given_variables=$(grep -v '#' $CONFIG | grep -v '^$' | cut -f1 -d ' ')
expected_variables=$(grep -v '#' $model_config | grep -v '^$' | cut -f1 -d ' ')

for var in $expected_variables ; do
  var_in_config=$(grep "^$var " $CONFIG)
  if [[ -z $var_in_config ]] ; then
    echo -e "\nERROR: The expected variable $var was not found in your config file (${CONFIG})."
    echo "Please make sure to include it, followed by a ':' and a space before the assigned value."
    echo "As a reminder:"
    echo "All expected variables must be specified in the yaml format, with one variable per row, for example : VARIABLE_NAME: value."
    echo "Some variables can be left empty, for example: VARIABLE_NAME: \"\" "
    echo -e "\nExiting.\n"
    exit 1
  fi
done


# Config file variables values

# FASTQ, FASTQ_R1 et FASTQ_R2 pas vides et au format attendu (cohérent avec PAIRED_END)
# OUTPUTS_DIRNAME pas vide
# DEMULT_DIR > si provided le chemin doit exister et contenir des FASTQ (cohérents avec PAIRED_END)
# ADAPT_FILE obligatoire -> doit exister, avoir le bon nombre de colonnes selon PAIRED_END, le bon séparateur, on peut lui enlever les \r
# BARCODE_FILE -> peut être vide si DEMULT_DIR est donné, sinon doit contenir le bon nombre de colonnes, le bon séparateur, on peut lui enlever les \r
# DEMULT_THREADS et DEMULT_SUBSTITUTIONS peuvent être vides si DEMULT_DIR est donné
# TRIMMING_THREADS, TRIMMING_QUAL et TRIMMING_MIN_LENGTH doivent être donnés
