#!/bin/bash

set -e -o pipefail

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
    --stats_folder)
    STATS_FOLDER="$2"
    shift
    shift
    ;;
    --output)
    OUTPUT="$2"
    shift
    shift
    ;;
esac
done


if [[ "$PAIRED" = "TRUE" ]] ; then
  echo -e "Sample\tR1_mapped\tReads_mapped\t%_reads_mapped" > $OUTPUT
  for stats in $(ls $STATS_FOLDER/stats_*) ; do
    sample=$(basename $stats | sed 's/stats_//')
    paired_mapped=$(cat $stats | grep 'mapped and paired:' | cut -f3)
    R1_mapped=$(echo "$paired_mapped/2" | bc)
    reads_mapped=$(cat $stats | grep 'reads mapped:' | cut -f3)
    total_reads=$(cat $stats | grep 'raw total sequences:' | cut -f3)
    prc_mapped=$(echo "scale=2; 100*$reads_mapped/$total_reads" | bc -l)
    echo -e $sample"\t"$R1_mapped"\t"$reads_mapped"\t"$prc_mapped >> $OUTPUT
  done
fi

if [[ "$PAIRED" = "FALSE" ]] ; then
  echo -e "Sample\tReads_mapped\t%_reads_mapped" > $OUTPUT
  for stats in $(ls $STATS_FOLDER/stats_*) ; do
    sample=$(basename $stats | sed 's/stats_//')
    reads_mapped=$(cat $stats | grep 'reads mapped:' | cut -f3)
    total_reads=$(cat $stats | grep 'raw total sequences:' | cut -f3)
    prc_mapped=$(echo "scale=2; 100*$reads_mapped/$total_reads" | bc -l)
    echo -e $sample"\t"$reads_mapped"\t"$prc_mapped >> $STATS_FOLDER/summary_stats.tsv
  done
fi
