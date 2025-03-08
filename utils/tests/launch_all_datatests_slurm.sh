#!/usr/bin/env bash

# Usage: launch_all_datatests_slurm.sh <path/to/GeCKO> agap_normal 1.2.0
# Depends on two files (test_subfolders.txt and test_subfolder_names.txt) expected in a CONFIG directory next to the script

set -euo pipefail


# ---- ARGUMENTS ---- #
GeCKOpath=$(realpath $1)
partition=$2
SingularityImageVersion=$3

# ------- CONFIG ------- #
configFilesDirPath=$(dirname $0)"/CONFIG"
dirsToCleanFile=${configFilesDirPath}/test_subfolders.txt
expectedTestSubfolderNamesFile=${configFilesDirPath}/test_subfolder_names.txt

debugTestsDir=${GeCKOpath}/utils/tests/data_tests

# --------- DEPENDENCIES ---------- #
source ${GeCKOpath}/utils/launching/utils.sh
source ${GeCKOpath}/utils/utils.sh


# -------- FUNCTIONS -------- #

importConfFiles(){
    mapfile -t dirsToClean < ${dirsToCleanFile}
    mapfile -t expectedTestSubfolderNames < ${expectedTestSubfolderNamesFile}
}

cleanOldFiles() {
    local GeCKOpath=$1
    cd $GeCKOpath
    for dir in "${dirsToClean[@]}" ; do
        rm -fr ${dir}/slurm*.out ${dir}/Logs_*Workflow
    done
    cdSilent -
}

sbatchRunGeCKO(){
    local runGeCKO=$1
    local WFname=$2
    local WFprefix=$3
    local partition=$4
    shift 4
    local extraOptions=("$@")
    
    runGeCKOcmd="sbatch --partition=${partition} --wrap=\"${runGeCKO} --workflow ${WFname} --config-file CONFIG/config_${WFname}.yml --cluster-profile CONFIG/${WFprefix}_CLUSTER_PROFILE_SLURM ${extraOptions}\""
    echo $runGeCKOcmd
    eval $runGeCKOcmd
}

sbatchRunGeCKOtest(){
    local runGeCKO=$1
    local WFname=$2
    local WFprefix=$3
    local partition=$4
    sbatchRunGeCKO $runGeCKO $WFname $WFprefix $partition "--jobs 20"
}

listTestSubfoldersInDir(){
    dir=$1
    find $dir -maxdepth 1 -type d -printf "%f\n" | grep -Fxf <(printf "%s\n" "${expectedTestSubfolderNames[@]}")
}


launchTestFolderTests(){
    local GeCKOpath=$1
    local testPath=$2
    local WFname=$3
    local WFprefix=$4
    local hasSubfolders=$5
    local partition=$6

    runGeCKO=${GeCKOpath}/runGeCKO.sh

    if [ $hasSubfolders == "TRUE" ] ; then
        subfolders=$(listTestSubfoldersInDir ${testPath})
    elif [ $hasSubfolders == "FALSE" ] ; then
        subfolders="."
    else
        echo "ERROR: Unknown hasSubfolders argument (expected TRUE or FALSE)"
        exit 1
    fi
    
    for subfolder in ${subfolders} ; do
        cd ${testPath}/${subfolder}
        echo -e "\n"${testPath}/${subfolder}":"
        sbatchRunGeCKOtest ${runGeCKO} ${WFname} ${WFprefix} ${partition}
    done
}




# ---------- MAIN --------- #

setErrorExitMsg

importConfFiles

cleanOldFiles $GeCKOpath

dlImageSylabs "library://ge2pop_gecko/gecko/gecko:${SingularityImageVersion}" "utils/singularity_image/GeCKO.sif"


declare -A testFolders
testFolders["DataCleaning_main"]="${GeCKOpath}/DATA_CLEANING/EXAMPLE DataCleaning DC TRUE" # fourth element = whether the workflow's test folder has subfolders or not
testFolders["ReadMapping_main"]="${GeCKOpath}/READ_MAPPING/EXAMPLE ReadMapping RM TRUE"
testFolders["VariantCalling_main"]="${GeCKOpath}/VARIANT_CALLING/EXAMPLE VariantCalling VC FALSE"
testFolders["VcfFiltering_main"]="${GeCKOpath}/VCF_FILTERING/EXAMPLE VcfFiltering VF FALSE"
testFolders["VcfFiltering_debug"]="${debugTestsDir}/VCF_FILTERING VcfFiltering VF TRUE"


for testFolder in "${!testFolders[@]}"; do
    read -r testPath WFname WFprefix hasSubfolders <<< "${testFolders[$testFolder]}"
    launchTestFolderTests $GeCKOpath $testPath $WFname $WFprefix $hasSubfolders $partition
done






