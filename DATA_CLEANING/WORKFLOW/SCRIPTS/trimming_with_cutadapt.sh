#!/usr/bin/env bash


#-----------------------------------------------------------------------------------------------------------
#
#    				TRIM A PAIR OF FASTQ FILES (R1 + R2) TO REMOVE ADAPTERS SEQUENCES AND LOW QUALITY READS
#
#-----------------------------------------------------------------------------------------------------------


# LAUNCHING EXAMPLE:
#./trimming_with_cutadapt.sh --ind Anvergur_A --trimdir {working_directory}/DEMULT_TRIM --R1 Anvergur_A.R1.fastq.gz --R2 Anvergur_A.R2.fastq.gz --adapt_file adapter_file_VIR_Cap001.txt --nodes 1 --qual 30 --min_length 36

  # Input file:
  #> adapter_file example (no header, tab-separated) >> Ind_name; R1_3'_adapter_seq; R2_3'_adapter_seq:
  #Anvergur_M ACTGCTTAGATCGGAAGAGCACACGTCTGAACTCCAGTCACATTACTCGATCTCGTATGCCGTCTTCTGCTTG	ACTCGTAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGGCTATAGTGTAGATCTCGGTGGTCGCCGTATCATT
  #Pescadou_M ATCGTGACAGATCGGAAGAGCACACGTCTGAACTCCAGTCACATTACTCGATCTCGTATGCCGTCTTCTGCTTG	GGACATCAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGGCTATAGTGTAGATCTCGGTGGTCGCCGTATCATT
  #Claudio_B2	AACTCGTAGATCGGAAGAGCACACGTCTGAACTCCAGTCACATTACTCGATCTCGTATGCCGTCTTCTGCTTG	TGCGCTAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGGCTATAGTGTAGATCTCGGTGGTCGCCGTATCATT

  # Output files:
  # ...

## Arguments description:
  # ...

#### ARGUMENTS:

POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
  --adapt_file)
  ADAPTERS_SEQ_FILE="$2"
  shift # past argument
  shift # past value
  ;;
  --ind)
  IND_NAME="$2"
  shift # past argument
  shift # past value
  ;;
  --trimdir)
  TRIM_DIR="$2"
  shift # past argument
  shift # past value
  ;;
  --R1)
  R1_DEMULT="$2"
  shift # past argument
  shift # past value
  ;;
  --R2)
  R2_DEMULT="$2"
  shift # past argument
  shift # past value
  ;;
  --nodes)
  NODES="$2"
  shift # past argument
  shift # past value
  ;;
  --qual)
  QC="$2"
  shift # past argument
  shift # past value
  ;;
  --min_length)
  MIN_LENGTH="$2"
  shift # past argument
  shift # past value
  ;;
  *)    # unknown option
  POSITIONAL+=("$1") # save it in an array for later
  shift # past argument
  ;;
esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters


### Manage file and folder paths (if relative path change it to absolute path)

if [[ ! "$ADAPTERS_SEQ_FILE" = /* ]] ; then
  ADAPTERS_SEQ_FILE=$(readlink -f $ADAPTERS_SEQ_FILE) ;
fi

if [[ ! "$TRIM_DIR" = /* ]] ; then
  TRIM_DIR=$(readlink -f $TRIM_DIR) ;
fi

if [[ ! "$R1_DEMULT" = /* ]] ; then
  R1_DEMULT=$(readlink -f $R1_DEMULT) ;
fi

if [[ ! "$R2_DEMULT" = /* ]] ; then
  R2_DEMULT=$(readlink -f $R2_DEMULT) ;
fi



#### SCRIPT
## Retrieve files and sequences
R1_readthrough_seq=$(grep -w $IND_NAME $ADAPTERS_SEQ_FILE | cut -f2)
R2_readthrough_seq=$(grep -w $IND_NAME $ADAPTERS_SEQ_FILE | cut -f3)

## Run cutadapt
cutadapt --action=trim --quality-cutoff $QC --minimum-length $MIN_LENGTH --no-indels -j $NODES -a $R1_readthrough_seq -A $R2_readthrough_seq -o ${TRIM_DIR}/${IND_NAME}_trimmed.R1.fastq.gz -p ${TRIM_DIR}/${IND_NAME}_trimmed.R2.fastq.gz $R1_DEMULT $R2_DEMULT > ${TRIM_DIR}/trimming_cutadapt_${IND_NAME}.info
