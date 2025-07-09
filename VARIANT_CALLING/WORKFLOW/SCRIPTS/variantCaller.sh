#!/bin/bash
#SBATCH --gres=gpu:a100:1
#SBATCH --threads=12

module () {
    eval `/usr/bin/modulecmd bash $*`
}

#pbrun haplotypecaller \
#    --ref "$1" \
#    --in-bam "$2" \
#    --out-variants "$3" \
#    #--haplotypecaller-options "$4"

pbrun deepvariant \
    --gvcf \
    --ref "$1" \
    --in-bam "$2" \
    --out-variants "$3" \
    --verbose \
    --mode "$4"
echo ""
echo "DeepVariant completed successfully."

# Remove the unwanted .vcf.gz file if it exists
BASENAME=$(basename "$3" .g.vcf.gz)
DIRNAME=$(dirname "$3")
rm -f "${DIRNAME}/${BASENAME}.vcf.gz"

# Index the output GVCF file
echo "Indexing GVCF file: $3"
pbrun indexgvcf \
    --input "$3"
echo ""
echo "Indexing of GVCF completed successfully."
