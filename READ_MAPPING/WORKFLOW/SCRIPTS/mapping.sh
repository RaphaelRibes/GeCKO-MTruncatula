#!/bin/bash

#{scripts_dir}/mapping.sh --paired_end --fastq_R1 {input.fastq_paired_R1} --fastq_R2 {input.fastq_paired_R2} --ref {input.ref} --mapper {mapper} --mapper_options {mapper_options} --technology {technology} --output_dir {bams_dir} --sample {wildcards.base}
#ou
#{scripts_dir}/mapping.sh --single_end --fastq {input.fastq_single} --ref {input.ref} --mapper {mapper} --mapper_options {mapper_options} --technology {technology} --output_dir {bams_dir} --sample {wildcards.base}


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
    --fastq_U)
    FASTQ_U="$2"
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
    --mapper_options)
    MAPPER_OPTIONS="$2"
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
esac
done


### Manage file and folder paths (if relative path change it to absolute path)

if [[ ! "$OUTPUT_DIR" = /* ]] ; then
  OUTPUT_DIR=$(readlink -f $OUTPUT_DIR) ;
fi
if [[ ! -z "$FASTQ_R1" && ! "$FASTQ_R1" = /* ]] ; then
  FASTQ_R1=$(readlink -f $FASTQ_R1) ;
fi
if [[ ! -z "$FASTQ_R2" && ! "$FASTQ_R2" = /* ]] ; then
  FASTQ_R2=$(readlink -f $FASTQ_R2) ;
fi
if [[ ! -z "$FASTQ_U" && ! "$FASTQ_U" = /* ]] ; then
  FASTQ_U=$(readlink -f $FASTQ_U) ;
fi
if [[ ! -z "$FASTQ" && ! "$FASTQ" = /* ]] ; then
  FASTQ=$(readlink -f $FASTQ) ;
fi

# Mapping

if [ "${MAPPER}" = "bwa-mem2_mem" ] ; then
  if [ "${PAIRED}" = "TRUE" ] ; then
    bwa-mem2 mem ${MAPPER_OPTIONS} -R $(echo "@RG\tID:${SAMPLE}\tPL:${TECHNOLOGY}\tSM:${SAMPLE}") ${REF} ${FASTQ_R1} ${FASTQ_R2} > ${OUTPUT_DIR}/${SAMPLE}_paired.sam || { (>&2 echo 'The mapping of ${SAMPLE} failed') ; exit 1; }
    if [ ! -z "$FASTQ_U" ] ; then
      bwa-mem2 mem ${MAPPER_OPTIONS} -R $(echo "@RG\tID:${SAMPLE}\tPL:${TECHNOLOGY}\tSM:${SAMPLE}") ${REF} ${FASTQ_U} > ${OUTPUT_DIR}/${SAMPLE}_unpaired.sam || { (>&2 echo 'The mapping of ${SAMPLE} failed') ; exit 1; }
      cat <(samtools view -h ${OUTPUT_DIR}/${SAMPLE}_paired.sam) <(samtools view ${OUTPUT_DIR}/${SAMPLE}_unpaired.sam) > ${OUTPUT_DIR}/${SAMPLE}.sam
    else
      mv ${OUTPUT_DIR}/${SAMPLE}_paired.sam ${OUTPUT_DIR}/${SAMPLE}.sam
    fi
  else
    bwa-mem2 mem ${MAPPER_OPTIONS} -R $(echo "@RG\tID:${SAMPLE}\tPL:${TECHNOLOGY}\tSM:${SAMPLE}") ${REF} ${FASTQ} > ${OUTPUT_DIR}/${SAMPLE}.sam || { (>&2 echo 'The mapping of ${SAMPLE} failed') ; exit 1; }
  fi
fi

if [ "${MAPPER}" = "bwa_mem" ] ; then
  if [ "${PAIRED}" = "TRUE" ] ; then
    bwa mem ${MAPPER_OPTIONS} -R $(echo "@RG\tID:${SAMPLE}\tPL:${TECHNOLOGY}\tSM:${SAMPLE}") ${REF} ${FASTQ_R1} ${FASTQ_R2} > ${OUTPUT_DIR}/${SAMPLE}_paired.sam || { (>&2 echo 'The mapping of ${SAMPLE} failed') ; exit 1; }
    if [ ! -z "$FASTQ_U" ] ; then
      bwa mem ${MAPPER_OPTIONS} -R $(echo "@RG\tID:${SAMPLE}\tPL:${TECHNOLOGY}\tSM:${SAMPLE}") ${REF} ${FASTQ_U} > ${OUTPUT_DIR}/${SAMPLE}_unpaired.sam
      cat <(samtools view -h ${OUTPUT_DIR}/${SAMPLE}_paired.sam) <(samtools view ${OUTPUT_DIR}/${SAMPLE}_unpaired.sam) > ${OUTPUT_DIR}/${SAMPLE}.sam
    else
      mv ${OUTPUT_DIR}/${SAMPLE}_paired.sam ${OUTPUT_DIR}/${SAMPLE}.sam
    fi
  else
    bwa mem ${MAPPER_OPTIONS} -R $(echo "@RG\tID:${SAMPLE}\tPL:${TECHNOLOGY}\tSM:${SAMPLE}") ${REF} ${FASTQ} > ${OUTPUT_DIR}/${SAMPLE}.sam || { (>&2 echo 'The mapping of ${SAMPLE} failed') ; exit 1; }
  fi
fi

if [ "${MAPPER}" = "bowtie2" ] ; then
  REF_INDEX=$(echo $REF | sed 's/.fasta//' | sed 's/.fas//' | sed 's/.fa//')
  if [ "${PAIRED}" = "TRUE" ] ; then
    bowtie2 ${MAPPER_OPTIONS} --rg-id ${SAMPLE} --rg "PL:${TECHNOLOGY}" --rg "SM:${SAMPLE}" -x ${REF_INDEX} -1 ${FASTQ_R1} -2 ${FASTQ_R2} -S ${OUTPUT_DIR}/${SAMPLE}_paired.sam || { (>&2 echo 'The mapping of ${SAMPLE} failed') ; exit 1; }
    if [ ! -z "$FASTQ_U" ] ; then
      bowtie2 ${MAPPER_OPTIONS} --rg-id ${SAMPLE} --rg "PL:${TECHNOLOGY}" --rg "SM:${SAMPLE}" -x ${REF_INDEX} -U ${FASTQ_U} -S ${OUTPUT_DIR}/${SAMPLE}_unpaired.sam
      cat <(samtools view -h ${OUTPUT_DIR}/${SAMPLE}_paired.sam) <(samtools view ${OUTPUT_DIR}/${SAMPLE}_unpaired.sam) > ${OUTPUT_DIR}/${SAMPLE}.sam
    else
      mv ${OUTPUT_DIR}/${SAMPLE}_paired.sam ${OUTPUT_DIR}/${SAMPLE}.sam
    fi
  else
    bowtie2 ${MAPPER_OPTIONS} --rg-id ${SAMPLE} --rg "PL:${TECHNOLOGY}" --rg "SM:${SAMPLE}" -x ${REF_INDEX} -U ${FASTQ} -S ${OUTPUT_DIR}/${SAMPLE}.sam || { (>&2 echo 'The mapping of ${SAMPLE} failed') ; exit 1; }
  fi
fi

if [ "${MAPPER}" = "minimap2" ] ; then
  REF_INDEX=$(echo $REF | sed 's/.fasta/.mmi/' | sed 's/.fas/.mmi/' | sed 's/.fa/.mmi/')
  if [ "${PAIRED}" = "TRUE" ] ; then
    minimap2 ${MAPPER_OPTIONS} -R $(echo "@RG\tID:${SAMPLE}\tPL:${TECHNOLOGY}\tSM:${SAMPLE}") -a ${REF_INDEX} ${FASTQ_R1} -2 ${FASTQ_R2} > ${OUTPUT_DIR}/${SAMPLE}_paired.sam || { (>&2 echo 'The mapping of ${SAMPLE} failed') ; exit 1; }
    if [ ! -z "$FASTQ_U" ] ; then
      minimap2 ${MAPPER_OPTIONS} -R $(echo "@RG\tID:${SAMPLE}\tPL:${TECHNOLOGY}\tSM:${SAMPLE}") -a ${REF_INDEX} ${FASTQ_U} > ${OUTPUT_DIR}/${SAMPLE}_unpaired.sam
      cat <(samtools view -h ${OUTPUT_DIR}/${SAMPLE}_paired.sam) <(samtools view ${OUTPUT_DIR}/${SAMPLE}_unpaired.sam) > ${OUTPUT_DIR}/${SAMPLE}.sam
    else
      mv ${OUTPUT_DIR}/${SAMPLE}_paired.sam ${OUTPUT_DIR}/${SAMPLE}.sam
    fi
  else
    minimap2 ${MAPPER_OPTIONS} -R $(echo "@RG\tID:${SAMPLE}\tPL:${TECHNOLOGY}\tSM:${SAMPLE}") -a ${REF_INDEX} ${FASTQ} > ${OUTPUT_DIR}/${SAMPLE}.sam || { (>&2 echo 'The mapping of ${SAMPLE} failed') ; exit 1; }
  fi
fi



samtools view -Sb -o ${OUTPUT_DIR}/${SAMPLE}.bam ${OUTPUT_DIR}/${SAMPLE}.sam || { (>&2 echo 'The transformation of ${SAMPLE} sam to bam failed') ; exit 1; }
rm -f ${OUTPUT_DIR}/${SAMPLE}_paired.sam ${OUTPUT_DIR}/${SAMPLE}_unpaired.sam ${OUTPUT_DIR}/${SAMPLE}.sam

# Fill in mate information
samtools sort -n ${OUTPUT_DIR}/${SAMPLE}.bam -o ${OUTPUT_DIR}/${SAMPLE}.sortn.bam || { (>&2 echo 'The samtools sort command of ${SAMPLE} failed') ; exit 1; }

if [ "${MAPPER}" = "bowtie2" ] ; then
  mv ${OUTPUT_DIR}/${SAMPLE}.sortn.bam ${OUTPUT_DIR}/${SAMPLE}.fix.bam
else
  #samtools fixmate -m ${OUTPUT_DIR}/${SAMPLE}.sortn.bam ${OUTPUT_DIR}/${SAMPLE}.fix.bam || { (>&2 echo 'The fixmate command of ${SAMPLE} failed') ; exit 1; }
  #rm ${OUTPUT_DIR}/${SAMPLE}.sortn.bam
  mv ${OUTPUT_DIR}/${SAMPLE}.sortn.bam ${OUTPUT_DIR}/${SAMPLE}.fix.bam
fi
rm ${OUTPUT_DIR}/${SAMPLE}.bam
