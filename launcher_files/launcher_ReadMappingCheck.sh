#!/bin/bash


# Model config file
model_config="${workflow_folder}/SCRIPTS/model_files/config_ReadMapping.yml"


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

## General variables (PAIRED_END) ##
PAIRED_END=$(grep "^PAIRED_END:" $CONFIG | sed 's/#.*$//' | cut -d ' ' -f2 | sed 's/"//g')


# **PAIRED_END**
if [[ -z "$PAIRED_END" ]] ; then
  echo -e "\nERROR: The PAIRED_END variable is missing in your config file (${CONFIG}). Please set it to TRUE or FALSE."
  echo "As a reminder:"
  echo "PAIRED_END: set to TRUE in case of paired end data (R1 + R2), to FALSE in case of single end data."
  echo -e "\nExiting.\n"
  exit 1
elif [[ "$PAIRED_END" == "TRUE" || "$PAIRED_END" == "True" || "$PAIRED_END" == "true" || "$PAIRED_END" == "T" ]] ; then
  echo -e "\nINFO: Paired end data is expected (PAIRED_END set to TRUE)\n"
elif [[ "$PAIRED_END" == "FALSE" || "$PAIRED_END" == "False" || "$PAIRED_END" == "false" || "$PAIRED_END" == "F" ]] ; then
  echo -e "\nINFO: Single end data is expected (PAIRED_END set to FALSE)\n"
else
  echo -e "\nERROR: The PAIRED_END variable is incorrect in your config file (${CONFIG}). Please set it to TRUE or FALSE."
  echo "As a reminder:"
  echo "PAIRED_END: set to TRUE in case of paired end data (R1 + R2), to FALSE in case of single end data."
  echo -e "\nExiting.\n"
  exit 1
fi



## Input data (TRIM_DIR, fastq files, REFERENCE) ##
TRIM_DIRS=$(grep "^TRIM_DIRS:" $CONFIG | sed 's/#.*$//' | cut -d '"' -f2)
REFERENCE=$(grep "^REFERENCE:" $CONFIG | sed 's/#.*$//' | cut -d ' ' -f2 | sed 's/"//g')


# **TRIM_DIR**
for TRIM_DIR in $TRIM_DIRS ; do
  if [[ -z "$TRIM_DIR" ]] ; then
    TRIM_DIR="${HERE}/WORKFLOWS_OUTPUTS/DATA_CLEANING/DEMULT_TRIM"
    if [[ ! -d "$TRIM_DIR" ]] ; then
      echo -e "\nERROR: You did not specify any TRIM_DIR folder in your config file (${CONFIG}), but the default TRIM_DIR folder (${TRIM_DIR}) does not exist. Please specify your TRIM_DIR path in the config_file."
      echo -e "\nExiting.\n"
      exit 1
    fi
  elif [[ ! -d "$TRIM_DIR" ]] ; then
    echo -e "\nERROR: The TRIM_DIR parameter (${TRIM_DIR}/) provided in your config file (${CONFIG}) is not valid. Please make sure the directory exists and the path is correctly written."
    echo -e "\nExiting.\n"
    exit 1
  fi

  if [[ "$PAIRED_END" == "TRUE" || "$PAIRED_END" == "True" || "$PAIRED_END" == "true" || "$PAIRED_END" == "T" ]] ; then
    nb_fastq_R1_obs=$(ls ${TRIM_DIR}/*.R1.fastq.gz 2>/dev/null | wc -l)
    nb_fastq_R2_obs=$(ls ${TRIM_DIR}/*.R2.fastq.gz 2>/dev/null | wc -l)
    if [[ $nb_fastq_R1_obs -eq 0 || $nb_fastq_R2_obs -eq 0 ]] ; then
      echo -e "\nERROR: The TRIM_DIR folder (${TRIM_DIR}/) is either empty or the fastq files it contains are not properly named. Input fastq files must end with '.R1.fastq.gz' and '.R2.fastq.gz'."
      echo -e "\nExiting.\n"
      exit 1
    else
      fastq_R1_list=$(ls -1 ${TRIM_DIR}/*.R1.fastq.gz 2>/dev/null | xargs -n1 basename 2>/dev/null)
      fastq_R2_list_exp=$(echo $fastq_R1_list | sed 's/.R1./.R2./g')
      nb_fastq_R2_exp=$(echo $fastq_R2_list_exp | wc -w)
      if [[ $nb_fastq_R2_obs != $nb_fastq_R2_exp ]] ; then
        echo -e "\nERROR: Input R1 and R2 fastq files do not seem to match. A set of matching *.R1.fastq.gz and *.R2.fastq.gz are expected in the TRIM_DIR folder (${TRIM_DIR})."
        echo -e "\nExiting.\n"
        exit 1
      fi
      for fastq_r2 in $fastq_R2_list_exp ; do
        if [[ ! -f ${TRIM_DIR}/$fastq_r2 ]] ; then
          echo -e "\nERROR: Input R1 and R2 fastq files do not seem to match. A set of matching *.R1.fastq.gz and *.R2.fastq.gz are expected in the TRIM_DIR folder (${TRIM_DIR})."
          echo -e "\nExiting.\n"
          exit 1
        fi
      done
    fi
  fi

  if [[ "$PAIRED_END" == "FALSE" || "$PAIRED_END" == "False" || "$PAIRED_END" == "false" || "$PAIRED_END" == "F" ]] ; then
    nb_fastq=$(ls ${TRIM_DIR}/*.fastq.gz 2>/dev/null | wc -l)
    nb_fastq_R1_obs=$(ls ${TRIM_DIR}/*R1*fastq.gz 2>/dev/null | wc -l)
    nb_fastq_R2_obs=$(ls ${TRIM_DIR}/*R2*fastq.gz 2>/dev/null | wc -l)
    if [[ "$nb_fastq" = 0 ]] ; then
      echo -e "\nERROR: The TRIM_DIR folder (${TRIM_DIR}) is either empty or the fastq files it contains are not properly named. Input fastq files must end with '.fastq.gz'."
      echo -e "\nExiting.\n"
      exit 1
    fi
    if [[ "$nb_fastq_R1_obs" > 0 && "$nb_fastq_R2_obs" > 0 ]] ; then
      echo -e "\nWARNING: You set PAIRED_END to FALSE but it appears that your data may be paired: 'R1' and 'R2' patterns were found in their names."
      echo -e "Your fastq files will be mapped as single end data. If your data is paired, please set PAIRED_END to TRUE.\n"
    fi
  fi
done

# **REFERENCE**
if [[ -z "$REFERENCE" ]] ; then
  echo -e "\nERROR: The REFERENCE variable is missing in your config file (${CONFIG}). Please provide a reference file to map your reads onto."
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

## Zones and extraction (BED, CREATE_SUB_BAMS)
BED=$(grep "^BED:" $CONFIG | sed 's/#.*$//' | cut -d ' ' -f2 | sed 's/"//g')
CREATE_SUB_BAMS=$(grep "^CREATE_SUB_BAMS:" $CONFIG | sed 's/#.*$//' | cut -d ' ' -f2 | sed 's/"//g')
BED_MIN_MEAN_COV=$(grep "^BED_MIN_MEAN_COV:" $CONFIG | sed 's/#.*$//' | cut -d ' ' -f2 | sed 's/"//g')
BED_MIN_DIST=$(grep "^BED_MIN_DIST:" $CONFIG | sed 's/#.*$//' | cut -d ' ' -f2 | sed 's/"//g')
BED_MIN_LENGTH=$(grep "^BED_MIN_LENGTH:" $CONFIG | sed 's/#.*$//' | cut -d ' ' -f2 | sed 's/"//g')

# **BED**
# Targeted zones bed file (optionnal) -> si donné existe-t-il, est-il au bon format
if [[ ! -z "$BED" && ! -f "$BED" ]] ; then
  echo -e "\nERROR: The BED file (${BED}) provided in your config file (${CONFIG}) does not exist. Please make sure the file exists and the path is correctly written, or leave the parameter empty (BED: \"\")."
  echo -e "\nExiting.\n"
  exit 1
fi
if [[ ! -z "$BED" && -f "$BED" ]] ; then
  nb_col_bed=$(awk '{print NF}' $BED | sort | uniq)
  if [[ $nb_col_bed != 3 ]] ; then
    echo -e "\nERROR: Your bed file (${BED}) does not seem to be properly formated."
    echo "As a reminder:"
    echo "The provided bed file must have 3 columns separated by tabs. Each row represents a genomic region of interest, with the first column indicating the reference's sequence name, and the second and third columns the starting and ending positions of the region in the sequence."
    echo -e "\nExiting.\n"
    exit 1
  fi
fi

# **CREATE_SUB_BAMS**
if [[ -z "$CREATE_SUB_BAMS" ]] ; then
  echo -e "\nERROR: The CREATE_SUB_BAMS variable is missing in your config file (${CONFIG}). Please set it to TRUE or FALSE."
  echo "As a reminder:"
  echo "CREATE_SUB_BAMS: set to TRUE in case you provided a bed file AND want to extract the reads mapping onto the specified regions, to FALSE otherwise. If set to TRUE, the extracted reads will be stored in new bams, and a new reference (matching the bams) containing only the bed regions will be created."
  echo -e "\nExiting.\n"
  exit 1
elif [[ "$CREATE_SUB_BAMS" != "TRUE" && "$CREATE_SUB_BAMS" != "True" && "$CREATE_SUB_BAMS" != "true" && "$CREATE_SUB_BAMS" != "T" && "$CREATE_SUB_BAMS" != "FALSE" && "$CREATE_SUB_BAMS" != "False" && "$CREATE_SUB_BAMS" != "false" && "$CREATE_SUB_BAMS" != "F" ]] ; then
  echo -e "\nERROR: The CREATE_SUB_BAMS variable is incorrect in your config file (${CONFIG}). Please set it to TRUE or FALSE."
  echo "As a reminder:"
  echo "CREATE_SUB_BAMS: set to TRUE in case you provided a bed file AND want to extract the reads mapping onto the specified regions, to FALSE otherwise. If set to TRUE, the extracted reads will be stored in new bams, and a new reference (matching the bams) containing only the bed regions will be created."
  echo -e "\nExiting.\n"
  exit 1
fi
if [[ "$CREATE_SUB_BAMS" == "TRUE" || "$CREATE_SUB_BAMS" == "True" || "$CREATE_SUB_BAMS" == "true" || "$CREATE_SUB_BAMS" == "T" ]] ; then
  if [[ -z "$BED" && (-z "$BED_MIN_MEAN_COV" || -z "$BED_MIN_DIST" || -z "$BED_MIN_LENGTH") ]] ; then
    echo -e "\nERROR: CREATE_SUB_BAM was set to TRUE but neither the bed file nor the parameters to automatically create it were provided in your config file (${CONFIG}). Please either provide a bed file, or values for BED_MIN_MEAN_COV, BED_MIN_DIST and BED_MIN_LENGTH, or set CREATE_SUB_BAM to FALSE."
    echo -e "\nExiting.\n"
    exit 1
  fi
  if [[ ! -z "$BED" && (! -z "$BED_MIN_MEAN_COV" || ! -z "$BED_MIN_DIST" || ! -z "$BED_MIN_LENGTH") ]] ; then
    echo -e "\nERROR: You provided both the bed file and the parameters to automatically create it in your config file (${CONFIG}). Please EITHER provide a bed file, OR values for BED_MIN_MEAN_COV, BED_MIN_DIST and BED_MIN_LENGTH."
    echo -e "\nExiting.\n"
    exit 1
  fi
fi


## Mapping parameters ##
MAPPER=$(grep "^MAPPER:" $CONFIG | sed 's/#.*$//' | cut -d ' ' -f2 | sed 's/"//g')
REMOVE_DUP=$(grep "^REMOVE_DUP:" $CONFIG | sed 's/#.*$//' | cut -d ' ' -f2 | sed 's/"//g')

# **MAPPER**
if [[ -z "$MAPPER" ]] ; then
  echo -e "\nERROR: The MAPPER variable is missing in your config file (${CONFIG}). Please provide it."
  echo "As a reminder:"
  echo "MAPPER: the mapper that will be used to map your reads. Currently supported options are 'bwa-mem2_mem', 'bwa_mem', 'bowtie2' and 'minimap2'."
  echo -e "\nExiting.\n"
  exit 1
fi
if [[ "$MAPPER" != "bwa-mem2_mem" && "$MAPPER" != "bwa_mem" && "$MAPPER" != "bowtie2" && "$MAPPER" != "minimap2" ]] ; then
  echo -e "\nERROR: The MAPPER (${MAPPER}) provided in your config file (${CONFIG}) is unknown."
  echo "As a reminder:"
  echo "MAPPER: the mapper that will be used to map your reads. Currently supported options are 'bwa-mem2_mem', 'bwa_mem', 'bowtie2' and 'minimap2'."
  echo -e "\nExiting.\n"
  exit 1
fi

# **REMOVE_DUP**
if [[ "$REMOVE_DUP" != "TRUE" && "$REMOVE_DUP" != "True" && "$REMOVE_DUP" != "true" && "$REMOVE_DUP" != "T" && "$REMOVE_DUP" != "FALSE" && "$REMOVE_DUP" != "False" && "$REMOVE_DUP" != "false" && "$REMOVE_DUP" != "F" ]] ; then
  echo -e "\nERROR: The REMOVE_DUP variable is incorrect in your config file (${CONFIG}). Please set it to TRUE or FALSE."
  echo "As a reminder:"
  echo "REMOVE_DUP: set to TRUE to remove duplicates after mapping (picard MarkDuplicates -REMOVE_DUPLICATES TRUE), to FALSE otherwise."
  echo -e "\nExiting.\n"
  exit 1
fi
