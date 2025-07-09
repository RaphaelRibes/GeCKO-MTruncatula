#!/bin/bash

ulimit -Sn 65536

mkdir -p "$1";
# glnexus_cli is the executable within the GLnexus container
glnexus_cli --list "$2" --config "$3" --thread "$4" > "$1"/joint_calls.bcf

bcftools convert -O z -o "$1"/joint_calls.vcf.gz ."$1"/joint_calls.bcf
bcftools index -t "$1"/joint_calls.vcf.gz