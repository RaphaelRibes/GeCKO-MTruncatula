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

# Index the output GVCF file
echo "Indexing GVCF file: $3"
pbrun indexgvcf \
    --input "$3"
echo ""
echo "Indexing of GVCF completed successfully."

# singularity exec utils/singularity_image/parabricks.sif ./VARIANT_CALLING/WORKFLOW/SCRIPTS/variantCaller.sh /lustre/ribesr/genom-asm-4-pg-haploid/results/ESP099_results/02_final_assembly/hap/ragtag_scafold/sorted/M.truncatula.fa /lustre/ribesr/GeCKO-MTruncatula/WORKFLOWS_OUTPUTS/READ_MAPPING/BAMS /lustre/ribesr/GeCKO-MTruncatula/test shortread