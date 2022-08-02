infile_vcf=$1
infile_extraHeader=$2
outfile=$3

#tmpVCFfile=$(mktemp __addPopGenStat.XXXXXXXXX)

#focus on SNP with exactly two alleles and add some extra information
  awk -v OFS="\t" '{ nbSamples=0; A1=0;A2=0;A1A1=0; A2A2=0; A1A2=0; biAllelic=1;\
    if (! ($1 ~ /^#/) )
    {\
      for (i=10; i<=NF; i++){ \
        nbSamples +=1;\
        split($i,tmp,":");geno=tmp[1]; \
        if (geno == "./.") miss += 1; \
        else if (geno == "0/0" ||geno == "0|0") {A1A1 += 1; A1+=2;} \
        else if (geno == "1/1" || geno == "1|1") {A2A2 += 1; A2+=2;} \
        else if (geno == "0/1" || geno == "1/0" || geno == "0|1" || geno == "1|0") {A1A2 += 1; A1+=1;A2+=1;}\
        else biAllelic=0;
      }\
      if( (biAllelic==1) && (A1>0) && (A2>0))
      {\
        nbG=A1A1+A1A2+A2A2; \
        p=A1/(A1+A2);\
        q=A2/(A1+A2);\
        He=(2*p*q);\
        F=(He-(A1A2/nbG))/He;\
        pcMiss=(nbSamples-nbG)/(nbSamples);\
        res="A1A1="A1A1";A2A2="A2A2";A1A2="A1A2";nbG="nbG";pcNA="pcMiss";p="p";q="q";He="He";F="F;\
        $8=$8";"res;\
        print; 
      }
    }
  else
    {
      print $0
    }
  }' $infile_vcf | awk -v toAdd=$infile_extraHeader 'BEGIN{first=1}/^##FORMAT/ {if(first==1){system (" cat "toAdd)} first=0} { print; }' > $outfile

   
    
    
    
#filter position poorly genotyped
#vcftools --vcf __indiv_filtered.recode.vcf --min-alleles 2 --max-missing 0.95 --maf 0.05 --recode --recode-INFO-all --out __indiv_site_filter --min-meanDP 20
#https://gist.github.com/darencard/1f7270dba015a1afa89eb7dd762f161e
#vcftools --recode --recode-INFO-all --stdout --vcf __indiv_site_filter.recode.vcf | \
#bcftools query -f '%CHROM\t%POS[\t%GT]\n' - | \

#snakemake  --use-singularity --use-conda --snakefile $snakePipe --configfile config.yml --jobs 200 --printshellcmds --use-envmodules 


