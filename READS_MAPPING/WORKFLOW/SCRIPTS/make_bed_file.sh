#!/bin/bash

#{scripts_dir}/make_bed_file.sh --input_bams_dir $input_bams_dir --output_dir $output_dir --min_cov 190 --max_dist 100 --min_length 100



#### ARGUMENTS:

while [[ $# -gt 0 ]]
do
  key="$1"
  case $key in
    --input_bams_dir)
    INPUT_BAMS_DIR="$2"
    shift
    shift
    ;;
    --output_dir)
    OUTPUT_DIR="$2"
    shift
    shift
    ;;
    --min_cov)
    MIN_COV="$2"
    shift
    shift
    ;;
    --max_dist)
    MAX_DIST="$2"
    shift
    shift
    ;;
    --min_length)
    MIN_LENGTH="$2"
    shift
    shift
    ;;
esac
done



### Manage file and folder paths (if relative path change it to absolute path)

if [[ ! "$INPUT_BAMS_DIR" = /* ]] ; then
  INPUT_BAMS_DIR=$(readlink -f $INPUT_BAMS_DIR) ;
fi
if [[ ! "$OUTPUT_DIR" = /* ]] ; then
  OUTPUT_DIR=$(readlink -f $OUTPUT_DIR) ;
fi



# ------------------------------------------------------------------------------------------------------ #


# Merge all bams into one
samtools merge ${OUTPUT_DIR}/all_merged.bam ${INPUT_BAMS_DIR}/*bam || { (>&2 echo 'The bed creation failed (bams merging step)') ; exit 1; }


# Compute coverage and only keep zones with enough reads
bedtools genomecov -ibam ${OUTPUT_DIR}/all_merged.bam -bg | awk -v m=$MIN_COV '{if ($4>=m) print $0}' > ${OUTPUT_DIR}/all_covered_zones_mincov.bed || { (>&2 echo 'The bed creation failed (genomecov step)') ; exit 1; }


# Merge overlapping and close enough zones
awk -v M=$MAX_DIST 'BEGIN{OFS="\t"; chr=0; start=0; end=0}{
  if ($1==chr && $2<=end+M) {
    end=$3
  }
  else {
    if(end !=0) {print chr, start, end} ;
    chr=$1 ; start=$2 ; end=$3
  }
} END{if(end !=0) {print chr, start, end}}' ${OUTPUT_DIR}/all_covered_zones_mincov.bed > ${OUTPUT_DIR}/all_covered_zones_mincov_collapsed.bed || { (>&2 echo 'The bed creation failed (zones merging step)') ; exit 1; }


# Remove zones that are too short
awk -v m=$MIN_LENGTH '{if ($3-$2+1 >= m){print $0}}' ${OUTPUT_DIR}/all_covered_zones_mincov_collapsed.bed > ${OUTPUT_DIR}/zones.bed || { (>&2 echo 'The bed creation failed (removing short zones step)') ; exit 1; }


rm ${OUTPUT_DIR}/all_merged.bam ${OUTPUT_DIR}/all_covered_zones_mincov.bed ${OUTPUT_DIR}/all_covered_zones_mincov_collapsed.bed
