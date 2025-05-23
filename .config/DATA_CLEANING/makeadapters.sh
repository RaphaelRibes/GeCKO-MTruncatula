rm -f adapters.txt
touch adapters.txt

if [ -z "$1" ]; then
    echo "Usage: $0 <directory>"
    exit 1
fi

for file in "$1"/*; do
    # Check if the file is a .fastq.gz file
    if [[ "$file" == *.fastq.gz ]]; then
        # If there is R2, skip
        if [[ "$file" == *".R2.fastq.gz" ]]; then
            continue
        fi
        # Remove .R1.fastq.gz from the filename
        base_file_name=$(basename "$file" | sed 's/.R1.fastq.gz//')
        # write to adapters.txt
        echo -e "${base_file_name}\tAGATCGGAAGAGCACACGTCTGAACTCCAGTCA\tAGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" >> adapters.txt
    fi
done