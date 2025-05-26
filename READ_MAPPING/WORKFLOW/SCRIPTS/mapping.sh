#!/bin/bash

# Example usage (as in your comments):
# {scripts_dir}/mapping.sh --fastq_R1 {input.fastq_paired_R1} --fastq_R2 {input.fastq_paired_R2} --ref {input.ref} --technology {technology} --output_dir {bams_dir} --sample {wildcards.base} --MD_options {params.picard_markduplicates_options} --MD_java_options {params.picard_markduplicates_java_options} --reports_dir {bams_reports_dir} --umi {dedupUMI} --gpu {num_gpus_for_parabricks}

set -e -o pipefail

#### ARGUMENTS:

# Initialize variables to prevent unbound variable errors later
FASTQ_R1=""
FASTQ_R2=""
FASTQ_U="" # Keeping for parsing if passed, but not used in Parabricks block
FASTQ=""   # Keeping for parsing if passed, but not used in Parabricks block
REF=""
TECHNOLOGY=""
OUTPUT_DIR=""
SAMPLE=""
MD_OPTIONS=""
MD_JAVA_OPTIONS=""
REPORTS_DIR=""
UMI=""
GPU="1" # Default to 1 GPU if not specified

while [[ $# -gt 0 ]]
do
  key="$1"

  case $key in
    --fastq_R1)
    FASTQ_R1="$2"
    shift # past argument
    shift # past value
    ;;
    --fastq_R2)
    FASTQ_R2="$2"
    shift # past argument
    shift # past value
    ;;
    --fastq_U) # This argument is parsed but its content won't be used by the Parabricks paired-end path
    FASTQ_U="$2"
    shift # past argument
    shift # past value
    ;;
    --fastq) # This argument is parsed but its content won't be used by the Parabricks paired-end path
    FASTQ="$2"
    shift # past argument
    shift # past value
    ;;
    --ref)
    REF="$2"
    shift # past argument
    shift # past value
    ;;
    --technology)
    TECHNOLOGY="$2"
    shift # past argument
    shift # past value
    ;;
    --output_dir)
    OUTPUT_DIR="$2"
    shift # past argument
    shift # past value
    ;;
    --sample)
    SAMPLE="$2"
    shift # past argument
    shift # past value
    ;;
    --MD_options)
    MD_OPTIONS="$2"
    shift # past argument
    shift # past value
    ;;
    --MD_java_options)
    MD_JAVA_OPTIONS="$2"
    shift # past argument
    shift # past value
    ;;
    --reports_dir)
    REPORTS_DIR="$2"
    shift # past argument
    shift # past value
    ;;
    --umi)
    UMI="$2"
    shift # past argument
    shift # past value
    ;;
    --gpu)
    GPU="$2"
    shift # past argument
    shift # past value
    ;;
    *) # unknown option
    echo "Unknown option: $1"
    exit 1
    ;;
esac
done


### Manage file and folder paths (if relative path change it to absolute path)
# All path conversion should happen before variables are used.
if [[ ! "$OUTPUT_DIR" = /* ]] ; then
  OUTPUT_DIR=$(readlink -f $OUTPUT_DIR) ;
fi
if [[ ! -z "$FASTQ_R1" && ! "$FASTQ_R1" = /* ]] ; then
  FASTQ_R1=$(readlink -f "$FASTQ_R1") ;
fi
if [[ ! -z "$FASTQ_R2" && ! "$FASTQ_R2" = /* ]] ; then
  FASTQ_R2=$(readlink -f "$FASTQ_R2") ;
fi

if [[ ! -z "$REPORTS_DIR" && ! "$REPORTS_DIR" = /* ]] ; then
  REPORTS_DIR=$(readlink -f "$REPORTS_DIR") ;
fi
if [[ ! "$REF" = /* ]] ; then # Ensure reference path is absolute
  REF=$(readlink -f "$REF") ;
fi


# --- Parabricks specific variables ---
# Construct Read Group (RG) arguments for Parabricks.
PARABRICKS_RG_ID="${SAMPLE}"
PARABRICKS_RG_SM="${SAMPLE}"
PARABRICKS_RG_PL="${TECHNOLOGY}" # e.g., ILLUMINA, ONT (ensure Parabricks supports this value)

# Parabricks options based on parsed arguments.
# Remove any BWA-specific mapper options that are no longer relevant.
PARABRICKS_OPTIONS="--num-gpus ${GPU} --num-threads 8" # Example: adjust --num-threads as needed

echo "Starting paired-end alignment for sample ${SAMPLE} using Parabricks fq2bam..."

# Parabricks fq2bam for paired-end reads.
# Output is directly a sorted BAM file.
parabricks fq2bam \
  --in-fq "${FASTQ_R1}" "${FASTQ_R2}" \
  --ref "${REF}" \
  --out-bam "${OUTPUT_DIR}/${SAMPLE}_raw.bam" \
  --read-group-id "${PARABRICKS_RG_ID}" \
  --read-group-sm "${PARABRICKS_RG_SM}" \
  --read-group-pl "${PARABRICKS_RG_PL}" \
  ${PARABRICKS_OPTIONS} # Pass Parabricks specific options

if [ $? -eq 0 ]; then
  echo "Successfully completed paired-end alignment for ${SAMPLE}. Output: ${OUTPUT_DIR}/${SAMPLE}_raw.bam"
  # IMPORTANT: Parabricks fq2bam outputs a sorted BAM. You likely need to index it.
  echo "Indexing the BAM file..."
  samtools index "${OUTPUT_DIR}/${SAMPLE}_raw.bam"
  if [ $? -ne 0 ]; then
    echo "Error indexing BAM file for ${SAMPLE}."
    exit 1
  fi
else
  echo "Error during Parabricks fq2bam execution for ${SAMPLE}. Please check logs."
  exit 1
fi

# --- Downstream Processing (Picard MarkDuplicates) ---
# This part now correctly receives a sorted BAM from Parabricks.

# Remove the redundant samtools sort command:
# samtools sort -o ${OUTPUT_DIR}/${SAMPLE}_sortcoord.bam ${OUTPUT_DIR}/${SAMPLE}_raw.sam

# Renamed output from _raw.sam to _raw.bam for clarity and consistency.
# Use the _raw.bam directly for Picard as it's already sorted.
# If you want to rename it to _sortcoord.bam for consistency with previous pipelines,
# you can add an 'mv' command here, but it's not strictly necessary.
# Example if you want to rename: mv "${OUTPUT_DIR}/${SAMPLE}_raw.bam" "${OUTPUT_DIR}/${SAMPLE}_sortcoord.bam"
# And then use "${OUTPUT_DIR}/${SAMPLE}_sortcoord.bam" in the Picard command.

echo "Running Picard MarkDuplicates..."
if [ "${UMI}" == "FALSE" ] ; then
  picard ${MD_JAVA_OPTIONS} MarkDuplicates \
    -I "${OUTPUT_DIR}/${SAMPLE}_raw.bam" \
    -O "${OUTPUT_DIR}/${SAMPLE}_markedDup.bam" \
    -VALIDATION_STRINGENCY SILENT \
    ${MD_OPTIONS} \
    -REMOVE_DUPLICATES FALSE \
    -M "${REPORTS_DIR}/DUPLICATES/${SAMPLE}.bam.metrics"
  if [ $? -ne 0 ]; then
    echo "Error during Picard MarkDuplicates for ${SAMPLE}."
    exit 1
  else
    echo "Successfully completed Picard MarkDuplicates for ${SAMPLE}."
  fi
fi

# Further processing like indexing the marked duplicates BAM, etc.
# samtools index "${OUTPUT_DIR}/${SAMPLE}_markedDup.bam"