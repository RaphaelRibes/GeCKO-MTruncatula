#!/bin/bash

vcf_gz_input=$1
tsv_output=$2
outdir=$3

# if the file is very big, sample 100000 rows
nb_rows=$(zcat $vcf_gz_input | grep -c -v '#')
if [[ $nb_rows -gt 100000 ]] ; then
  zcat $vcf_gz_input | grep -v '#' | awk -v rows=$nb_rows 'BEGIN {srand()} {if (rand() <= 100000/rows) print $0}' > ${vcf_gz_input}_sample
  gzip ${vcf_gz_input}_sample
  vcf_gz_input=${vcf_gz_input}_sample.gz
fi

# list all variables from column 8
zcat $vcf_gz_input | grep -v '#' | cut -f8 | tr ';' '\n' | cut -d '=' -f1 | sort | uniq > ${outdir}/tmp_variables
cat ${outdir}/tmp_variables | tr '\n' '\t' | sed 's/\t$//' | awk '{print "Contig\tPos\tQual\t"$0}' > $tsv_output


# for every site retrieve the values for all variables (and mark it NA if the info is missing)
nvar=$(cat ${outdir}/tmp_variables | wc -l)
awk -v nvar=$nvar -F"\t|;" '{
  if (NR==FNR)
    {vars[NR]=$1}
  else {
    print $1"\t"$2"\t"$3 ;
    for(j=1;j<=nvar;j++){
      for(i=4;i<=NF;i++){
        if($i ~ "^"vars[j]"=") {sub(vars[j]"=","",$i); print $i; break} else if (i == NF){print "NA"}
      }
    }
  }
}' ${outdir}/tmp_variables <(zcat $vcf_gz_input | grep -v '#' | cut -f1,2,6,8) | awk -v nvar=$nvar 'NR % (nvar+1) {printf("%s\t", $0); next} {print $0}' >> $tsv_output

rm ${outdir}/tmp_variables
rm -f $1"_sample.gz"
