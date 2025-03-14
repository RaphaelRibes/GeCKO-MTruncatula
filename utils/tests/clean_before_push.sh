#!/usr/bin/env bash

# Usage: clean_before_push.sh --gecko-path .
# Depends on a test_subfolders.txt file expected in a CONFIG directory next to the script

set -euo pipefail

# to-do:
# also rm reference index files for READ_MAPPING datasets and .dict, .fai and Reference_zones_intervals_for_GATK.list for the VARIANT_CALLING ones

# ------- CONFIG ------- #

configPath=$(dirname $0)"/CONFIG"
dirsFile=${configPath}/test_subfolders.txt


# -------- ARGUMENTS -------- #

while [[ $# -gt 0 ]]
    do
    key="$1"
    case $key in
        --gecko-path)
        geckoPath=$(realpath "$2")
        shift
        shift
        ;;
        -*|*)
        echo -e "\nWARNING: $1 option is unknown and will be ignored.\n"
        shift
        ;;
    esac
done

# -------- FUNCTIONS -------- #

sourceDependencies(){
    local geckoPath=$1
    source ${geckoPath}/utils/utils.sh
}

cleanDir(){
    local testDir=$1
    cd ${testDir}
    echo "Cleaning ${testDir}..."
    rm -rf Logs_* .snakemake .java .config .conda WORKFLOWS_OUTPUTS slurm*
    cdSilent -
}

cleanDirs(){
    echo -e "\n"
    cleanDir ${geckoPath}/utils/getToolVersions
    for testDir in "${testDirs[@]}" ; do
        cleanDir ${testDir}
    done
}

# ---------- MAIN --------- #

sourceDependencies $geckoPath

setErrorExitMsg

importTestDirs $dirsFile

cleanDirs
