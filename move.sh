#!/bin/bash
#
# This script copies GeCKO workflow outputs (READ_MAPPING or VARIANT_CALLING)
# to a specified results directory.
# It requires specifying the data type and the target reference folder.
#

#SBATCH --mem=1G
#SBATCH --job-name=move_gecko_data
#SBATCH --move-%j.out
#SBATCH --move-%j.err

# --- Usage Function ---
# Displays help information if the script is called incorrectly.
usage() {
    echo "A script to copy GeCKO mapping or calling data."
    echo ""
    echo "Usage: $0 [-m | -c] [-r | -a] FOLDER_NAME"
    echo "------------------------------------------------------------------"
    echo "  Data Type (select one):"
    echo "    -m : Copy READ_MAPPING data."
    echo "    -c : Copy VARIANT_CALLING data."
    echo ""
    echo "  Target Reference (select one):"
    echo "    -r : Target the 'A17-vs-long_reads' results folder."
    echo "    -a : Target the 'assembly-vs-long_reads' results folder."
    echo ""
    echo "  Required Argument:"
    echo "    FOLDER_NAME : The name of the source sample folder (e.g., F11023)."
    echo "------------------------------------------------------------------"
    exit 1
}

# --- Argument Parsing ---
# Initialize variables to hold the selected options.
MODE=""
REF=""

# Use getopts to parse command-line flags.
while getopts "mcra" opt; do
    case ${opt} in
        m)
            # Ensure -m and -c are not used together.
            [ -n "$MODE" ] && { echo "Error: Cannot specify both -m and -c."; usage; }
            MODE="mapping"
            ;;
        c)
            # Ensure -c and -m are not used together.
            [ -n "$MODE" ] && { echo "Error: Cannot specify both -c and -m."; usage; }
            MODE="calling"
            ;;
        r)
            # Ensure -r and -a are not used together.
            [ -n "$REF" ] && { echo "Error: Cannot specify both -r and -a."; usage; }
            REF="A17"
            ;;
        a)
            # Ensure -a and -r are not used together.
            [ -n "$REF" ] && { echo "Error: Cannot specify both -a and -r."; usage; }
            REF="assembly"
            ;;
        \?)
            # Handle invalid options.
            usage
            ;;
    esac
done
# Shift the parsed options out of the argument list.
shift "$((OPTIND-1))"

# --- Argument Validation ---
# Verify that one mode option was chosen.
if [ -z "$MODE" ]; then
    echo "Error: You must specify a data type (-m for mapping or -c for calling)."
    usage
fi

# Verify that one reference option was chosen.
if [ -z "$REF" ]; then
    echo "Error: You must specify a target reference (-r for A17 or -a for assembly)."
    usage
fi

# Verify that the FOLDER_NAME argument was provided.
if [ -z "$1" ]; then
    echo "Error: You must provide a FOLDER_NAME as the final argument."
    usage
fi

FOLDER_NAME="$1"

# --- Path Construction ---
# Determine the source and destination directory names based on the selected options.

# Set source directory name based on the mode.
if [ "$MODE" == "mapping" ]; then
    SOURCE_TYPE="READ_MAPPING"
elif [ "$MODE" == "calling" ]; then
    SOURCE_TYPE="VARIANT_CALLING"
fi

# Set destination directory names based on mode and reference.
if [ "$MODE" == "mapping" ]; then
    DEST_RESULTS_DIR="04_results_on_mapping"
elif [ "$MODE" == "calling" ]; then
    DEST_RESULTS_DIR="05_results_on_calling"
fi

if [ "$REF" == "A17" ]; then
    DEST_REF_DIR="A17-vs-long_reads"
elif [ "$REF" == "assembly" ]; then
    DEST_REF_DIR="assembly-vs-long_reads"
fi

# Assemble the full source and destination paths.
SOURCE_PATH="/storage/simple/users/ribesr/nv_scratch/GeCKO/GeCKO-MTruncatula-${FOLDER_NAME}-lr/WORKFLOWS_OUTPUTS/${SOURCE_TYPE}"
DEST_PATH="/storage/replicated/cirad/projects/GE2POP/2024_AGRODIV/02_results/${DEST_RESULTS_DIR}/${DEST_REF_DIR}/${FOLDER_NAME}"

# --- Execution ---
echo "================================================="
echo "Starting data copy process..."
echo "  - Data Type:   ${MODE} (${SOURCE_TYPE})"
echo "  - Reference:   ${REF} (${DEST_REF_DIR})"
echo "  - Sample Folder: ${FOLDER_NAME}"
echo ""
echo "  - Source:      ${SOURCE_PATH}"
echo "  - Destination: ${DEST_PATH}"
echo "================================================="

# Check if the source directory actually exists before proceeding.
if [ ! -d "$SOURCE_PATH" ]; then
    echo "Error: Source directory not found: $SOURCE_PATH"
    exit 1
fi

# If ${DEST_PATH} does not exist, send an error message.
if [ ! -d "$DEST_PATH" ]; then
    echo "Error: Destination directory does not exist: $DEST_PATH"
    exit 1
fi


# Execute the copy command. The source directory itself will be copied into the destination.
echo "Copying data..."
cp -r "${SOURCE_PATH}" "${DEST_PATH}"

echo ""
echo "Copy complete."
echo "================================================="
