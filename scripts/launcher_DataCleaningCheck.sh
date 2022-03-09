#!/bin/bash

# Model config file
if [[ "$WORKFLOW_SMK" = "${WORKFLOW}_PairedEnd.smk" ]] ; then
  echo -e "\nINFO: Paired end data is expected (PAIRED_END set to TRUE)\n"
  model_config="${workflow_folder}/SCRIPTS/model_files/config_DataCleaning_PE.yml"
else
  echo -e "\nINFO: Single end data is expected (PAIRED_END set to FALSE)\n"
  model_config="${workflow_folder}/SCRIPTS/model_files/config_DataCleaning_SE.yml"
fi


# Config file variables
given_variables=$(grep -v '^#' $CONFIG | grep -v '^$' | cut -f1 -d ' ')
expected_variables=$(grep -v '^#' $model_config | grep -v '^$' | cut -f1 -d ' ')

for var in $expected_variables ; do
  var_in_config=$(grep "^$var " $CONFIG)
  if [[ -z $var_in_config ]] ; then
    echo -e "\nERROR: The expected variable $var was not found in your config file (${CONFIG}). Please make sure to include it."
    echo -e "\nList of expected variables for PAIRED_END set to ${PAIRED_END} (some can be left empty but must appear in the file nonetheless) :\n${expected_variables}"
    echo -e "\nAs a reminder:"
    echo "All expected variables must be specified in YAML format, with one variable per row followed by a ':' and a space before its assigned value. For example: VARIABLE_NAME: value."
    echo "To leave a variable empty, use quotes: VARIABLE_NAME: \"\" "
    echo -e "\nExiting.\n"
    exit 1
  fi
done


### Config file variables values

# OUTPUTS_DIRNAME
OUTPUTS_DIRNAME=$(grep "^OUTPUTS_DIRNAME:" $CONFIG | sed 's/#.*$//' | cut -d ' ' -f2 | sed 's/"//g')
if [[ -z $OUTPUTS_DIRNAME ]] ; then
  echo -e "\nERROR: You must provide the OUTPUTS_DIRNAME in the config_file."
  echo -e "\nExiting.\n"
  exit 1
fi


# FASTQ, FASTQ_R1, FASTQ_R2 and DEMULT_DIR
FASTQ_R1=$(grep "^FASTQ_R1:" $CONFIG | sed 's/#.*$//' | cut -d ' ' -f2 | sed 's/"//g')
FASTQ_R2=$(grep "^FASTQ_R2:" $CONFIG | sed 's/#.*$//' | cut -d ' ' -f2 | sed 's/"//g')
DEMULT_DIR=$(grep "^DEMULT_DIR:" $CONFIG | sed 's/#.*$//' | cut -d ' ' -f2 | sed 's/"//g')
FASTQ=$(grep "^FASTQ:" $CONFIG | sed 's/#.*$//' | cut -d ' ' -f2 | sed 's/"//g')
FASTQ_R1=$(absolutePath $FASTQ_R1)
FASTQ_R2=$(absolutePath $FASTQ_R2)
DEMULT_DIR=$(absolutePath $DEMULT_DIR)
FASTQ=$(absolutePath $FASTQ)


if [[ "$WORKFLOW_SMK" = "${WORKFLOW}_PairedEnd.smk" ]] ; then # if paired end data is expected
  if [[ (! -z "$FASTQ_R1" || ! -z "$FASTQ_R2") && ! -z "$DEMULT_DIR" ]] ; then
    echo -e "\nERROR: You provided both FASTQ_R1 and/or FASTQ_R2 and DEMULT_DIR in the config_file. Please provide either FASTQ_R1 and FASTQ_R2 in case of multiplexed data, or DEMULT_DIR in case of demultiplexed data."
    echo -e "\nExiting.\n"
    exit 1
  fi
  if [[ -z "$FASTQ_R1" || -z "$FASTQ_R2" ]] ; then
    if [[ -z "$DEMULT_DIR" ]] ; then
      echo -e "\nERROR: You must either provide FASTQ_R1 and FASTQ_R2 or DEMULT_DIR in the config_file, but they were all left blank."
      echo -e "\nExiting.\n"
      exit 1
    elif [[ ! -d "$DEMULT_DIR" ]] ; then
      echo -e "\nERROR: the DEMULT_DIR parameter provided in the config file (${DEMULT_DIR}) is not a valid path. Please make sure the directory exists and the path is correctly written."
      echo -e "\nExiting.\n"
      exit 1
    else
      nb_fastq=$(ls ${DEMULT_DIR}/*fastq.gz 2>/dev/null | wc -l)
      if [[ "$nb_fastq" = 0 ]] ; then
        echo -e "\nERROR: The provided DEMULT_DIR (${DEMULT_DIR}) is either empty or the fastq files it contains are not properly named. Input demultiplexed fastq files must end with '.fastq.gz'."
        echo -e "\nExiting.\n"
        exit 1
      else
        fastq_R1_list=$(ls -1 ${DEMULT_DIR}/*R1.fastq.gz 2>/dev/null | xargs -n1 basename)
        nb_fastq_R2_obs=$(ls ${DEMULT_DIR}/*R2.fastq.gz 2>/dev/null | wc -l)
        fastq_R2_list_exp=$(echo $fastq_R1_list | sed 's/.R1./.R2./g')
        nb_fastq_R2_exp=$(echo $fastq_R2_list_exp | wc -w)
        for fastq_r2 in $fastq_R2_list_exp ; do
          if [[ ! -f ${DEMULT_DIR}/$fastq_r2 || $nb_fastq_R2_obs != $nb_fastq_R2_exp ]] ; then
            echo -e "\nERROR: R1 and R2 fastq files provided in the config file do not seem to match. A set of matching *.R1.fastq.gz and *.R2.fastq.gz are expected in the DEMULT folder you provided (${DEMULT_DIR})."
            echo -e "\nExiting.\n"
            exit 1
          fi
        done
      fi
    fi
  elif [[ ! -f "$FASTQ_R1" || ! -f "$FASTQ_R2" ]] ; then
    echo -e "\nERROR: At least one the two input fastq files ${FASTQ_R1} and ${FASTQ_R2} provided in the config file does not exist. Please make sure the files exist and their path is correct."
    echo -e "\nExiting.\n"
    exit 1
  elif [[ "$FASTQ_R1" != *_R1.fastq.gz || "$FASTQ_R2" != *_R2.fastq.gz ]] ; then
      echo -e "\nERROR: Provided FASTQ_R1 and FASTQ_R2 (${FASTQ_R1} and ${FASTQ_R2}) are not properly named. Multiplexed fastq files must end with _R1.fastq.gz and _R2.fastq.gz."
      echo -e "\nExiting.\n"
      exit 1
  else
    R1_basename=$(basename $FASTQ_R1 _R1.fastq.gz)
    R2_basename=$(basename $FASTQ_R2 _R2.fastq.gz)
    if [[ "$R1_basename" != "$R2_basename" ]] ; then
      echo -e "\nERROR: Provided FASTQ_R1 and FASTQ_R2 names do not match: ${R1_basename} != ${R2_basename}."
      echo -e "\nExiting.\n"
      exit 1
    fi
  fi
else # if single end data is expected
  if [[ ! -z "$FASTQ" && ! -z "$DEMULT_DIR" ]] ; then
    echo -e "\nERROR: You provided both FASTQ and DEMULT_DIR in the config_file. Please provide either FASTQ in case of multiplexed data, or DEMULT_DIR in case of demultiplexed data."
    echo -e "\nExiting.\n"
    exit 1
  fi
  if [[ -z "$FASTQ" ]] ; then
    DEMULT_DIR=$(grep "^DEMULT_DIR:" $CONFIG | sed 's/#.*$//' | cut -d ' ' -f2 | sed 's/"//g')
    if [[ -z "$DEMULT_DIR" ]] ; then
      echo -e "\nERROR: You must provide either FASTQ or DEMULT_DIR in the config_file, but they were both left blank."
      echo -e "\nExiting.\n"
      exit 1
    elif [[ ! -d "$DEMULT_DIR" ]] ; then
      echo -e "\nERROR: the DEMULT_DIR parameter provided in the config file (${DEMULT_DIR}) is not a valid path. Please make sure the directory exists and the path is correctly written."
      echo -e "\nExiting.\n"
      exit 1
    else
      nb_fastq=$(ls ${DEMULT_DIR}/*fastq.gz 2>/dev/null | wc -l)
      if [[ "$nb_fastq" = 0 ]] ; then
        echo -e "\nERROR: The provided DEMULT_DIR (${DEMULT_DIR}) is either empty or the fastq files it contains are not properly named. Input demultiplexed fastq files must end with '.fastq.gz'."
        echo -e "\nExiting.\n"
        exit 1
      fi
    fi
  elif [[ ! -f "$FASTQ" ]] ; then
    echo -e "\nERROR: The input fastq file ${FASTQ} provided in the config file does not exist. Please make sure the file exists and its path is correct."
    echo -e "\nExiting.\n"
    exit 1
  elif [[ "$FASTQ" != *.fastq.gz ]] ; then
      echo -e "\nERROR: Provided FASTQ (${FASTQ}) is not properly named. Multiplexed fastq files must end with .fastq.gz"
      echo -e "\nExiting.\n"
      exit 1
  fi
fi

# DEMULT_THREADS and DEMULT_SUBSTITUTIONS
DEMULT_THREADS=$(grep "^DEMULT_THREADS:" $CONFIG | sed 's/#.*$//' | cut -d ' ' -f2 | sed 's/"//g')
DEMULT_SUBSTITUTIONS=$(grep "^DEMULT_SUBSTITUTIONS:" $CONFIG | sed 's/#.*$//' | cut -d ' ' -f2 | sed 's/"//g')
if [[ -z "$DEMULT_DIR" && -z "$DEMULT_THREADS" ]] ; then
  echo -e "\nERROR: A DEMULT_THREADS value must be provided in the config file."
  echo -e "\nExiting.\n"
  exit 1
fi


if [[ -z "$DEMULT_DIR" && -z "$DEMULT_SUBSTITUTIONS" ]] ; then
  echo -e "\nERROR: A DEMULT_SUBSTITUTIONS value must be provided in the config file."
  echo -e "\nExiting.\n"
  exit 1
fi


# BARCODE_FILE and ADAPT_FILE
BARCODE_FILE=$(grep "^BARCODE_FILE:" $CONFIG | sed 's/#.*$//' | cut -d ' ' -f2 | sed 's/"//g')
ADAPT_FILE=$(grep "^ADAPT_FILE:" $CONFIG | sed 's/#.*$//' | cut -d ' ' -f2 | sed 's/"//g')
BARCODE_FILE=$(absolutePath $BARCODE_FILE)
ADAPT_FILE=$(absolutePath $ADAPT_FILE)

if [[ -z "$DEMULT_DIR" && -z "$BARCODE_FILE" ]] ; then
  echo -e "\nERROR: A BARCODE_FILE must be provided in the config file."
  echo -e "\nExiting.\n"
  exit 1
fi
if [[ ! -z "$BARCODE_FILE" && ! -f "$BARCODE_FILE" ]] ; then
  echo -e "\nERROR: The BARCODE_FILE provided in the config file (${BARCODE_FILE}) does not exist. Please make sure the file exists and its path is correct."
  echo -e "\nExiting.\n"
  exit 1
fi

if [[ ! -z "$BARCODE_FILE" ]] ; then
  nb_carriage_returns=$(grep -c $'\r' $BARCODE_FILE)
  if [[ "$nb_carriage_returns" -gt 0 ]] ; then
    echo "Removing windows carriage returns in ${BARCODE_FILE}..."
    sed -i 's/\r$//g' $BARCODE_FILE
    sed -i 's/\r/\n/g' $BARCODE_FILE
  fi
  nb_col_barcode_file=$(awk -F '\t' '{print NF}' $BARCODE_FILE | sort | uniq)
  if [[ "$WORKFLOW_SMK" = "${WORKFLOW}_PairedEnd.smk" && $nb_col_barcode_file != 3 ]] ; then
    echo -e "\nERROR: The barcode file (${BARCODE_FILE}) provided in the config file does not have a proper format. For paired end data, 3 columns are expected, separated by tabs."
    echo -e "\nExiting.\n"
    exit 1
  elif [[ "$WORKFLOW_SMK" = "${WORKFLOW}_SingleEnd.smk" && $nb_col_barcode_file != 2 ]] ; then
    echo -e "\nERROR: The barcode file (${BARCODE_FILE}) provided in the config file does not have a proper format. For single end data, 2 columns are expected, separated by tabs."
    echo -e "\nExiting.\n"
    exit 1
  fi
fi


if [[ -z "$ADAPT_FILE" ]] ; then
  echo -e "\nERROR: An ADAPT_FILE must be provided in the config file."
  echo -e "\nExiting.\n"
  exit 1
fi
if [[ ! -z "$ADAPT_FILE" && ! -f "$ADAPT_FILE" ]] ; then
  echo -e "\nERROR: The ADAPT_FILE provided in the config file (${ADAPT_FILE}) does not exist. Please make sure the file exists and its path is correct."
  echo -e "\nExiting.\n"
  exit 1
fi
nb_carriage_returns=$(grep -c $'\r' $ADAPT_FILE)
if [[ "$nb_carriage_returns" -gt 0 ]] ; then
  echo "Removing windows carriage returns in ${ADAPT_FILE}..."
  sed -i 's/\r$//g' $ADAPT_FILE
  sed -i 's/\r/\n/g' $ADAPT_FILE
fi

nb_col_adapt_file=$(awk -F '\t' '{print NF}' $ADAPT_FILE | sort | uniq)
if [[ "$WORKFLOW_SMK" = "${WORKFLOW}_PairedEnd.smk" && $nb_col_adapt_file != 3 ]] ; then
  echo -e "\nERROR: The adapter file (${ADAPT_FILE}) provided in the config file does not have a proper format. For paired end data, 3 columns are expected, separated by tabs."
  echo -e "\nExiting.\n"
  exit 1
elif [[ "$WORKFLOW_SMK" = "${WORKFLOW}_SingleEnd.smk" && $nb_col_adapt_file != 2 ]] ; then
  echo -e "\nERROR: The adapter file (${ADAPT_FILE}) provided in the config file does not have a proper format. For single end data, 2 columns are expected, separated by tabs."
  echo -e "\nExiting.\n"
  exit 1
fi

# TRIMMING_THREADS, TRIMMING_QUAL and TRIMMING_MIN_LENGTH
TRIMMING_THREADS=$(grep "^TRIMMING_THREADS:" $CONFIG | sed 's/#.*$//' | cut -d ' ' -f2 | sed 's/"//g')
if [[ -z $TRIMMING_THREADS ]] ; then
  echo -e "\nERROR: You must provide the TRIMMING_THREADS in the config_file."
  echo -e "\nExiting.\n"
  exit 1
fi

TRIMMING_QUAL=$(grep "^TRIMMING_QUAL:" $CONFIG | sed 's/#.*$//' | cut -d ' ' -f2 | sed 's/"//g')
if [[ -z $TRIMMING_QUAL ]] ; then
  echo -e "\nERROR: You must provide the TRIMMING_QUAL in the config_file."
  echo -e "\nExiting.\n"
  exit 1
fi

TRIMMING_MIN_LENGTH=$(grep "^TRIMMING_MIN_LENGTH:" $CONFIG | sed 's/#.*$//' | cut -d ' ' -f2 | sed 's/"//g')
if [[ -z $TRIMMING_MIN_LENGTH ]] ; then
  echo -e "\nERROR: You must provide the TRIMMING_MIN_LENGTH in the config_file."
  echo -e "\nExiting.\n"
  exit 1
fi
