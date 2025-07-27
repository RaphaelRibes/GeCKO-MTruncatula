#!/bin/bash
# -i: input directory
# -c: data cleaning
# -m: read mapping
# -v: variant calling
# -f: variant filtering
# -j: number of jobs

SOURCE_DIR="/storage/replicated/cirad/projects/GE2POP/2024_AGRODIV/01_raw_data/01_preliminary_data/RAW_FASTQS"
dirname=""
reference=false
data_cleaning=false
read_mapping=false
variant_calling=false
variant_filtering=false
jobs=100

while getopts "i:rcmvfaj:" opt; do
  case $opt in
    i) dirname="$OPTARG" ;;
    r) reference=true ;;
    c) data_cleaning=true ;;
    m) read_mapping=true ;;
    v) variant_calling=true ;;
    f) variant_filtering=true ;;
    a) data_cleaning=true
       read_mapping=true
       variant_calling=true
       variant_filtering=true ;;
    j) jobs="$OPTARG" ;;
    *) echo "Usage: $0 [-i input_dir] [-c] [-m] [-v] [-f] [-j jobs]"
       exit 1 ;;
  esac
done

input_dir="$SOURCE_DIR/$dirname"

if [ "$data_cleaning" = "true" ]; then
  echo "Running data cleaning..."
  chmod +x .config/DATA_CLEANING/makeadapters.sh
  chmod +x runGeCKO.sh
  cd .config/DATA_CLEANING/ && ./makeadapters.sh "$input_dir" && cd ../..
  ./runGeCKO.sh --workflow DataCleaning --config-file .config/DATA_CLEANING/config.yml --cluster-profile .config/DATA_CLEANING/SLURM/ --jobs $jobs
fi

if [ "$reference" = true ]; then
  reference=/storage/replicated/cirad/projects/GE2POP/REFERENCES/MEDICAGO/truncat/genomeA17v5/MtrunA17r5.0-20161119-ANR.genome.fasta
else
  scaffold_dir=/storage/replicated/cirad/projects/GE2POP/2024_AGRODIV/02_results/03_results_on_assembly/03_scaffolds/"$dirname"/sorted
  reference_file="$dirname"-M.truncatula.fa
  reference=$scaffold_dir/$reference_file
fi
echo "Using reference: $reference"

if [ "$read_mapping" = "true" ]; then
  echo "Running read mapping..."
  # Change the path of the reference in .config/READ_MAPPING/config.yml to the assembly directory
  sed -i "s|REFERENCE: .*|REFERENCE: $reference|" .config/READ_MAPPING/config.yml
  ./runGeCKO.sh --workflow ReadMapping --config-file .config/READ_MAPPING/config.yml --cluster-profile .config/READ_MAPPING/SLURM/ --jobs $jobs
fi

if [ "$variant_calling" = "true" ]; then
  echo "Running variant calling..."
  # Change the path of the reference in .config/VARIANT_CALLING/config.yml to the assembly directory
  sed -i "s|REFERENCE: .*|REFERENCE: $reference|" .config/VARIANT_CALLING/config.yml
  ./runGeCKO.sh --workflow VariantCalling --config-file .config/VARIANT_CALLING/config.yml --cluster-profile .config/VARIANT_CALLING/SLURM/ --jobs $jobs
fi

if [ "$variant_filtering" = "true" ]; then
  echo "Running variant filtering..."
  ./runGeCKO.sh --workflow VCF_FILTERING --config-file .config/VCF_FILTERING/config.yml --cluster-profile .config/VCF_FILTERING/SLURM/ --jobs $jobs
fi