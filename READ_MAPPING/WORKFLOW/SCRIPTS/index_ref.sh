#!/bin/bash

#./index_ref.sh {input.ref} {mapper}

REF=$1
MAPPER=$2

# Index reference
if [[ "$MAPPER" == "bwa-mem2_mem" ]] ; then
  bwa-mem2 index $REF  || { (>&2 echo 'The indexing of $REF failed') ; exit 1; }
fi

if [[ "$MAPPER" == "bwa_mem" ]] ; then
  bwa index $REF || { (>&2 echo 'The indexing of $REF failed') ; exit 1; }
fi

if [[ "$MAPPER" == "bowtie2" ]] ; then
  REF_INDEX=$(echo $REF | sed 's/.fasta//' | sed 's/.fas//' | sed 's/.fa//')
  bowtie2-build $REF $REF_INDEX || { (>&2 echo 'The indexing of $REF failed') ; exit 1; }
fi

if [[ "$MAPPER" == "minimap2" ]] ; then
  REF_INDEX=$(echo $REF | sed 's/.fasta/.mmi/' | sed 's/.fas//' | sed 's/.fa/.mmi/')
  minimap2 -d $REF_INDEX $REF || { (>&2 echo 'The indexing of $REF failed') ; exit 1; }
fi
