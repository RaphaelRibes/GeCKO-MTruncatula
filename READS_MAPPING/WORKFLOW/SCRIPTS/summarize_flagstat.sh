#!/bin/bash

while [[ $# -gt 0 ]]
do
  key="$1"
  case $key in
    --pairedEnd)
    PAIRED="TRUE"
    shift # past argument
    ;;
    --singleEnd)
    PAIRED="FALSE"
    shift # past argument
    ;;
    --flagstat_folder)
    FLAGSTAT_FOLDER="$2"
    shift
    shift
    ;;
esac
done



if [[ "$PAIRED" = "TRUE" ]] ; then
  echo -e "Sample\tR1_mapped\tReads_mapped\t%_reads_mapped" > $FLAGSTAT_FOLDER/summary_flagstat.tsv
  for flagstat in $(ls $FLAGSTAT_FOLDER/flagstat_*) ; do
    sample=$(basename $flagstat | sed 's/flagstat_//')
    paired_mapped=$(cat $flagstat | grep 'with itself and mate mapped' | cut -f1 -d' ')
    R1_mapped=$(echo "$paired_mapped/2" | bc)
    reads_mapped=$(cat $flagstat | grep 'primary mapped (' | cut -f1 -d' ')
    prc_mapped=$(awk 'FNR==7' $flagstat | sed 's/.*(//' | sed 's/%.*//')
    echo -e $sample"\t"$R1_mapped"\t"$reads_mapped"\t"$prc_mapped >> $FLAGSTAT_FOLDER/summary_flagstat.tsv
  done
fi

if [[ "$PAIRED" = "FALSE" ]] ; then
  echo -e "Sample\tReads_mapped\t%_reads_mapped" > $FLAGSTAT_FOLDER/summary_flagstat.tsv
  for flagstat in $(ls $FLAGSTAT_FOLDER/flagstat_*) ; do
    sample=$(basename $flagstat | sed 's/flagstat_//')
    reads_mapped=$(cat $flagstat | grep 'primary mapped (' | cut -f1 -d' ')
    prc_mapped=$(awk 'FNR==7' $flagstat | sed 's/.*(//' | sed 's/%.*//')
    echo -e $sample"\t"$reads_mapped"\t"$prc_mapped >> $FLAGSTAT_FOLDER/summary_flagstat.tsv
  done
fi
