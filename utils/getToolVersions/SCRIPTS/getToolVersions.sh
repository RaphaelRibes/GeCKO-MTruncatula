input=$1
output=$2


isPicardInList=$(grep "^picard$" $input | wc -l)
if [ $isPicardInList -eq 1 ] ; then
    echo -e -n "picard\t" >> ${output}
    find $(conda info --base) -name "picard.jar" | grep "GeCKO_env" | sed -e 's,/picard.jar,,' -e 's,.*/,,'  >> ${output}
fi

isBwaInList=$(grep "^bwa$" $input | wc -l)
if [ $isBwaInList -eq 1 ] ; then
    echo -e -n "bwa\t" >> ${output}
    bwa 2>&1 | grep "Version:"  >> ${output}
fi

isBwaMem2InList=$(grep "^bwa-mem2$" $input | wc -l)
if [ $isBwaMem2InList -eq 1 ] ; then
    echo -e -n "bwa-mem2\t" >> ${output}
    bwa-mem2 version >> ${output}
fi


for tool in $(grep -E -v "\b(picard|bwa|bwa-mem2)\b" ${input}) ; do
    echo -e -n ${tool}"\t" >> ${output}
    ${tool} --version 2>&1 > /dev/null
    if [ $? -eq 0 ] ; then
        ${tool} --version | head -1 >> ${output}
    fi
done
