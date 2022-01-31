#!/usr/bin/env bash


#-----------------------------------------------------------------------------------------------------------
#
#    					COUNT READS IN ALL FASTQ.GZ FILES OF FOLDER
#
#-----------------------------------------------------------------------------------------------------------


# LAUNCHING EXAMPLE:
#./script_fastq_read_count.sh --fastq_dir /home/ardissonm/scratch/VIR/DEMULT --output /home/ardissonm/scratch/VIR/RESULTS/DEMULT/fastq_read_count_demult.txt

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
    --fastq_dir)
    FASTQ_DIR="$2"
    shift # past argument
    shift # past value
    ;;
	--output)
    OUTPUT="$2"
    shift # past argument
    shift # past value
    ;;
    --default)
    DEFAULT=YES
    shift # past argument
    ;;
    *)    # unknown option
    POSITIONAL+=("$1") # save it in an array for later
    shift # past argument
    ;;
esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters


### Manage file and folder paths (if relative path change it to absolute path)

if [[ ! "$FASTQ_DIR" = /* ]] ; then  
  FASTQ_DIR=$(readlink -f $FASTQ_DIR) ;
fi
if [[ ! "$OUTPUT" = /* ]] ; then 	
  OUTPUT=$(readlink -f $OUTPUT) ; 
fi



#### SCRIPT:

for gz in ${FASTQ_DIR}/*fastq.gz ; do rows=$(gunzip -c $gz | wc -l) ; ((reads=${rows}/4)) ; name=$(basename $gz) ; echo -e ${name}"\t"${reads} ; done | sort > ${OUTPUT}
