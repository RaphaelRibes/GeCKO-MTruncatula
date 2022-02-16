#!/usr/bin/env bash

#---------------------------------------------------#
#													#
#    		demultiplex_with_cutadapt_PE.sh			#
#													#
#---------------------------------------------------#

# >>> USAGE CONTEXT:
# github: https://github.com/BioInfo-GE2POP-BLE/CAPTURE_PIPELINES_SNAKEMAKE
# This script is intended to be used by the Snakemake workflow "DATA_CLEANING" / Sequencing type: Paired end
# tools: Cutadapt (DOI:10.14806/ej.17.1.200) https://cutadapt.readthedocs.io/en/stable/guide.html#demultiplexing

# >>> OBJECTIVE(S):
# demultiplex a pair of fastq files (R1 + R2) into individual paired fastq files
# Assignment of sequences to each of the genotypes according to the barcode or tag, located at the P5 and P7 ends, which were assigned to them during the construction of the libraries.

# >>> LAUNCHING EXAMPLE AND SETTINGS:
#./demultiplex_with_cutadapt.sh --demultdir {demult_dir} --R1 DEV.R1.fastq.gz --R2 DEV.R2.fastq.gz --tag_file sample_file_DEV.txt --nodes 1 --substitutions 0.15

#--demultdir
	# folders that will contain the demultiplexing output files
#--R1
	# path to the fastq.gz file corresponding to the R1 sequences 
#--R2
	# path to the fastq.gz file corresponding to the R2 sequences
#--tag_file
	# path to the sample_file.txt containing the list of samples names (column 1) , sequences of barcode/tags for read 1 (P5) (column 2) and sequences of barcode/tags for read 2 (P7) (column 3)
	# example (no header, tab-separated):
	#Tc2208a	AGCGCA	AGCGCA
	#Tc2235a	CTCAGC	CTCAGC
	#Tc2249a	TAGATC	TAGATC		
#--nodes
	# number of nodes to be allocated on cluster 
#--substitutions
	# percentage of substitution by barcode (tag). Example: 1 substitution on a barcode of 8pb, note 0.15
# this script launches cutapdat with two options by default: --no-indels --pair-adapters

# >>> OUTPUTS
# two files by sample, listed in the sample file, contain sequences : sample.R1.fastq.gz and sample.R1.fastq.gz >>> storage in: {OUTPUTS_DIRNAME}/DATA_CLEANING/DEMULT
# a report created automatically by Cutadap: demultiplexing_cutadapt.info >>> storage in: {OUTPUTS_DIRNAME}/DATA_CLEANING/DEMULT/REPORTS/CUTADAPT_INFOS


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
  *)    			 # unknown option
  POSITIONAL+=("$1") # save it in an array for later
  shift 			 # past argument
  ;;
esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters


# Manage file and folder paths (if relative path change it to absolute path)

if [[ ! "$DEMULT_DIR" = /* ]] ; then  			
  DEMULT_DIR=$(readlink -f $DEMULT_DIR) ;   
fi
if [[ ! "$R1" = /* ]] ; then  			
  R1=$(readlink -f $R1) ;   
fi
if [[ ! "$R2" = /* ]] ; then
  R2=$(readlink -f $R2) ;
fi
if [[ ! "$TAG_FILE" = /* ]] ; then
  TAG_FILE=$(readlink -f $TAG_FILE) ;
fi


#### SCRIPT:

## 1/ CREATE BARCODE/TAGS FASTA FILES
mkdir ${DEMULT_DIR}/Cutadapt_tmp
awk '{print ">"$1"\n^"$2}' ${TAG_FILE} > ${DEMULT_DIR}/Cutadapt_tmp/R1_tags_cutadapt.fasta
awk '{print ">"$1"\n^"$3}' ${TAG_FILE} > ${DEMULT_DIR}/Cutadapt_tmp/R2_tags_cutadapt.fasta

## 2/ RUN CUTADAPT - DEMULTIPLEXING
cutadapt -e ${SUBSTITUTIONS} --no-indels --pair-adapters -j ${NODES} -g file:${DEMULT_DIR}/Cutadapt_tmp/R1_tags_cutadapt.fasta -G file:${DEMULT_DIR}/Cutadapt_tmp/R2_tags_cutadapt.fasta -o ${DEMULT_DIR}/{name}.R1.fastq.gz -p ${DEMULT_DIR}/{name}.R2.fastq.gz ${R1} ${R2} > ${DEMULT_DIR}/demultiplexing_cutadapt.info

## 3/ REMOVE TMP FILES
rm -r ${DEMULT_DIR}/Cutadapt_tmp
