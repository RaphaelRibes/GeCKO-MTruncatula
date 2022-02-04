#!/usr/bin/env bash


#-----------------------------------------------------------------------------------------------------------
#
#    					DEMULTIPLEX A PAIR OF FASTQ FILES (R1 + R2) INTO INDIVIDUAL PAIRED FASTQ FILES
#
#-----------------------------------------------------------------------------------------------------------


# LAUNCHING EXAMPLE:
#./demultiplex_with_cutadapt.sh --demultdir {working_directory}/DEMULT --R1 DEV_Cap005_sub100000_R1.fastq.gz --R2 DEV_Cap005_sub100000_R2.fastq.gz --tag_file tag_file_DEV_Cap005_DT.txt --nodes 1 --auth_subst 0.15
  # Input files:
  #> tag_file example (no header, tab-separated):
  #Tc2208	AGCGCA	AGCGCA
  #Tc2235	CTCAGC	CTCAGC
  #Tc2249	TAGATC	TAGATC

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
  --demultdir)
  DEMULT_DIR="$2"
  shift # past argument
  shift # past value
  ;;
  --R1)
  R1="$2"
  shift # past argument
  shift # past value
  ;;
  --R2)
  R2="$2"
  shift # past argument
  shift # past value
  ;;
  --tag_file)
  TAG_FILE="$2"
  shift # past argument
  shift # past value
  ;;
  --nodes)
  NODES="$2"
  shift # past argument
  shift # past value
  ;;
  --substitutions)
  SUBSTITUTIONS="$2"
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

if [[ ! "$DEMULT_DIR" = /* ]] ; then  			# if the path is relative
  DEMULT_DIR=$(readlink -f $DEMULT_DIR) ;   # make it absolute
fi
if [[ ! "$R1" = /* ]] ; then  			# if the path is relative
  R1=$(readlink -f $R1) ;   # make it absolute
fi
if [[ ! "$R2" = /* ]] ; then
  R2=$(readlink -f $R2) ;
fi
if [[ ! "$TAG_FILE" = /* ]] ; then
  TAG_FILE=$(readlink -f $TAG_FILE) ;
fi


#### SCRIPT:

## 1/ CREATE TAGS FASTA FILES
mkdir ${DEMULT_DIR}/Cutadapt_tmp
awk '{print ">"$1"\n^"$2}' $TAG_FILE > ${DEMULT_DIR}/Cutadapt_tmp/R1_tags_cutadapt.fasta
awk '{print ">"$1"\n^"$3}' $TAG_FILE > ${DEMULT_DIR}/Cutadapt_tmp/R2_tags_cutadapt.fasta


## 2/ DEMULTIPLEX
cutadapt -e $SUBSTITUTIONS --no-indels --pair-adapters -j $NODES -g file:${DEMULT_DIR}/Cutadapt_tmp/R1_tags_cutadapt.fasta -G file:${DEMULT_DIR}/Cutadapt_tmp/R2_tags_cutadapt.fasta -o ${DEMULT_DIR}/{name}.R1.fastq.gz -p ${DEMULT_DIR}/{name}.R2.fastq.gz $R1 $R2 > ${DEMULT_DIR}/demultiplexing_cutadapt.info

## 3/ REMOVE TMP FILES
rm -r ${DEMULT_DIR}/Cutadapt_tmp
