#!/bin/bash

#{scripts_dir}/mapping.sh --paired_end --fastq_R1 {input.fastq_paired_R1} --fastq_R2 {input.fastq_paired_R2} --ref {input.ref} --mapper {mapper} --mapper_params {mapper_params} --technology {technology} --output_dir {bams_dir} --reports_dir {bams_reports_dir} --sample {wildcards.base} --rm_dup {rm_dup} --markdup_params {} --index_params {}
#ou
#{scripts_dir}/mapping.sh --single_end --fastq {input.fastq_single} --ref {input.ref} --mapper {mapper} --mapper_params {mapper_params} --technology {technology} --output_dir {bams_dir} -reports_dir {bams_reports_dir} --sample {wildcards.base} --rm_dup {rm_dup} --markdup_params {markdup_params} --index_params {index_params}


#### ARGUMENTS:

while [[ $# -gt 0 ]]
do
  key="$1"

  case $key in
    --paired_end)
    PAIRED="TRUE"
    shift # past argument
    ;;
    --single_end)
    PAIRED="FALSE"
    shift # past argument
    ;;
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
    --fastq)
    FASTQ="$2"
    shift # past argument
    shift # past value
    ;;
    --ref)
    REF="$2"
    shift # past argument
    shift # past value
    ;;
    --mapper)
    MAPPER="$2"
    shift # past argument
    shift # past value
    ;;
    --mapper_params)
    MAPPER_PARAMS="$2"
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
    --reports_dir)
    REPORTS_DIR="$2"
    shift # past argument
    shift # past value
    ;;
    --sample)
    SAMPLE="$2"
    shift # past argument
    shift # past value
    ;;
    --rm_dup)
    RM_DUP="$2"
    shift # past argument
    shift # past value
    ;;
    --markdup_params)
    MARKDUP_PARAMS="$2"
    shift # past argument
    shift # past value
    ;;
    --index_params)
    INDEX_PARAMS="$2"
    shift # past argument
    shift # past value
    ;;
esac
done


### Manage file and folder paths (if relative path change it to absolute path)

if [[ ! "$OUTPUT_DIR" = /* ]] ; then
  OUTPUT_DIR=$(readlink -f $OUTPUT_DIR) ;
fi
if [[ ! "$REPORTS_DIR" = /* ]] ; then
  REPORTS_DIR=$(readlink -f $REPORTS_DIR) ;
fi
if [[ ! -z "$FASTQ_R1" && ! "$FASTQ_R1" = /* ]] ; then
  FASTQ_R1=$(readlink -f $FASTQ_R1) ;
fi
if [[ ! -z "$FASTQ_R2" && ! "$FASTQ_R2" = /* ]] ; then
  FASTQ_R2=$(readlink -f $FASTQ_R2) ;
fi
if [[ ! -z "$FASTQ" && ! "$FASTQ" = /* ]] ; then
  FASTQ=$(readlink -f $FASTQ) ;
fi


# Mapping

if [ "${MAPPER}" = "bwa-mem2_mem" ] ; then
  if [ "${PAIRED}" = "TRUE" ] ; then
    bwa-mem2 mem ${MAPPER_PARAMS} -R $(echo "@RG\tID:${SAMPLE}\tPL:${TECHNOLOGY}\tSM:${SAMPLE}") ${REF} ${FASTQ_R1} ${FASTQ_R2} > ${OUTPUT_DIR}/${SAMPLE}.sam
  else
    bwa-mem2 mem ${MAPPER_PARAMS} -R $(echo "@RG\tID:${SAMPLE}\tPL:${TECHNOLOGY}\tSM:${SAMPLE}") ${REF} ${FASTQ} > ${OUTPUT_DIR}/${SAMPLE}.sam
  fi
fi

if [ "${MAPPER}" = "bwa_mem" ] ; then
  if [ "${PAIRED}" = "TRUE" ] ; then
    bwa mem ${MAPPER_PARAMS} -R $(echo "@RG\tID:${SAMPLE}\tPL:${TECHNOLOGY}\tSM:${SAMPLE}") ${REF} ${FASTQ_R1} ${FASTQ_R2} > ${OUTPUT_DIR}/${SAMPLE}.sam
  else
    bwa mem ${MAPPER_PARAMS} -R $(echo "@RG\tID:${SAMPLE}\tPL:${TECHNOLOGY}\tSM:${SAMPLE}") ${REF} ${FASTQ} > ${OUTPUT_DIR}/${SAMPLE}.sam
  fi
fi

if [ "${MAPPER}" = "bowtie2" ] ; then
  REF_INDEX=$(echo $REF | sed 's/.fasta//' | sed 's/.fa//')
  if [ "${PAIRED}" = "TRUE" ] ; then
    #bowtie2 --no-unal -p n -x index_name -1 reads_1.fastq -2 reads_2.fastq -S output.sam
    bowtie2 ${MAPPER_PARAMS} --rg-id ${SAMPLE} --rg "PL:${TECHNOLOGY}" --rg "SM:${SAMPLE}" -x ${REF_INDEX} -1 ${FASTQ_R1} -2 ${FASTQ_R2} -S ${OUTPUT_DIR}/${SAMPLE}.sam
  else
    bowtie2 ${MAPPER_PARAMS} --rg-id ${SAMPLE} --rg "PL:${TECHNOLOGY}" --rg "SM:${SAMPLE}" -x ${REF_INDEX} -U ${FASTQ} -S ${OUTPUT_DIR}/${SAMPLE}.sam
  fi
fi

if [ "${MAPPER}" = "minimap2" ] ; then
  REF_INDEX=$(echo $REF | sed 's/.fasta/.mmi/' | sed 's/.fa/.mmi/')
  if [ "${PAIRED}" = "TRUE" ] ; then
    minimap2 ${MAPPER_PARAMS} -R $(echo "@RG\tID:${SAMPLE}\tPL:${TECHNOLOGY}\tSM:${SAMPLE}") -a ${REF_INDEX} ${FASTQ_R1} -2 ${FASTQ_R2} > ${OUTPUT_DIR}/${SAMPLE}.sam
  else
    #minimap2 -x splice:hq -t8 -a $REF $file > ${OUTPUT_DIR}/${SAMPLE}.sam
    minimap2 ${MAPPER_PARAMS} -R $(echo "@RG\tID:${SAMPLE}\tPL:${TECHNOLOGY}\tSM:${SAMPLE}") -a ${REF_INDEX} ${FASTQ} > ${OUTPUT_DIR}/${SAMPLE}.sam
  fi
fi

#else
#  echo "The provided mapper is unknown. Implemented mappers are 'bwa-mem2_mem', 'bowtie' and 'minimap2'."
#fi


samtools view -Sb -o ${OUTPUT_DIR}/${SAMPLE}.bam ${OUTPUT_DIR}/${SAMPLE}.sam
rm ${OUTPUT_DIR}/${SAMPLE}.sam

# Fill in mate information
if [ "${MAPPER}" = "bowtie2" ] ; then
  mv ${OUTPUT_DIR}/${SAMPLE}.bam ${OUTPUT_DIR}/${SAMPLE}.fix.bam
else
  samtools fixmate ${OUTPUT_DIR}/${SAMPLE}.bam ${OUTPUT_DIR}/${SAMPLE}.fix.bam
  rm ${OUTPUT_DIR}/${SAMPLE}.bam
fi

# Sort
picard SortSam -I ${OUTPUT_DIR}/${SAMPLE}.fix.bam -O ${OUTPUT_DIR}/${SAMPLE}.sort.bam -SO coordinate -VALIDATION_STRINGENCY SILENT
rm ${OUTPUT_DIR}/${SAMPLE}.fix.bam

# Remove duplicates
if [ "${RM_DUP}" = "True" ] ; then
	mkdir -p ${REPORTS_DIR}
  mkdir -p ${REPORTS_DIR}/DUPLICATES
  picard MarkDuplicates -I ${OUTPUT_DIR}/${SAMPLE}.sort.bam -O ${OUTPUT_DIR}/${SAMPLE}.bam -VALIDATION_STRINGENCY SILENT ${MARKDUP_PARAMS} -REMOVE_DUPLICATES TRUE -M ${REPORTS_DIR}/DUPLICATES/${SAMPLE}.bam.metrics
  rm ${OUTPUT_DIR}/${SAMPLE}.sort.bam
else
	mv ${OUTPUT_DIR}/${SAMPLE}.sort.bam ${OUTPUT_DIR}/${SAMPLE}.bam
fi


# Create index
samtools index ${INDEX_PARAMS} ${OUTPUT_DIR}/${SAMPLE}.bam
