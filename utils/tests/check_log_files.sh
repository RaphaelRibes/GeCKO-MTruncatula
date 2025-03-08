#!/usr/bin/env bash

# Usage: check_log_files.sh <path/to/GeCKO>
# Depends on a test_subfolders.txt file expected in a CONFIG directory next to the script

set -euo pipefail

# ---- ARGUMENTS ---- #

GeCKOpath=$(realpath $1)


# --------- DEPENDENCIES ---------- #

source ${GeCKOpath}/utils/utils.sh



# -------- FUNCTIONS -------- #

printSmkLogErrors(){
    echo -e "\n"
    for testSubfolder in "${testSubfolders[@]}" ; do
        for log in ${testSubfolder}/slurm*out ; do
            echo $log":"
            grep error $log || true
            grep '%' $log | tail -1 || true
            echo -e "\n"
        done
    done
}


# ---------- MAIN --------- #

setErrorExitMsg

mapfile -t testSubfolders < $(dirname $0)"/CONFIG/test_subfolders.txt"

printSmkLogErrors



