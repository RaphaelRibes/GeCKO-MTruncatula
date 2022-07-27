# 20220719

module load egglib/3.1  

cd /home/ardissonm/scratch/DEV/ANALYSES_2022_WORKFLOW/DEV_Cap009_and_Cap010/WORKFLOWS_OUTPUTS/VCF_FILTERING/EGGLIB

03_PopGenStatsSampleLocus_Filtered.vcf
egglib_vcf_stats.py
labels_groups.txt

python egglib_vcf_stats.py



# avec egglib_vcf_stats2.py


python egglib_vcf_stats.py --bed_file /home/ardissonm/scratch/DEV/ANALYSES_2022_WORKFLOW/DEV_Cap009_and_Cap010/WORKFLOWS_OUTPUTS/READS_MAPPING/MAPPING_ZAVITAN_WITH_DUPLICATES/BAMS/all_covered_zones_119_extra100pb_collapsed_min100.bed --vcf_file 03_PopGenStatsSampleLocus_Filtered.vcf --labels_file labels_groups.txt

