import egglib
import argparse

bases = set('ACGT')

# interface
parser = argparse.ArgumentParser(description="")
parser.add_argument("-i", "--input", dest="input", help="input VCF file", required=True)
parser.add_argument("-o", "--output", dest="output", help="output VCF file", required=True)
args = parser.parse_args()

# open VCF
vcf = egglib.io.VcfParser(args.input)
input_vcf = open(args.input)

# open output VCF and write header
header = """##INFO=<ID=SNP,Number=0,Type=Integer,Description="Variant with at least 2 mononucleotide alleles">
##INFO=<ID=INDEL,Number=0,Type=Integer,Description="Variant with at least 2 alleles of different lengths">
##INFO=<ID=nbGS,Number=1,Type=Integer,Description="Number of genotyped samples">
##INFO=<ID=pNA,Number=1,Type=Float,Description="Frequency of missing genotypes for this site">
##INFO=<ID=nbG,Number=1,Type=Integer,Description="Number of genotypes">
##INFO=<ID=nbAll,Number=1,Type=Integer,Description="Number of alleles">
##INFO=<ID=All,Number=.,Type=String,Description="list of alleles">
##INFO=<ID=AllCount,Number=.,Type=Int,Description="Counts of all alleles">
##INFO=<ID=AllFreq,Number=.,Type=Float,Description="Frequency of all alleles">
##INFO=<ID=MinHomo,Number=1,Type=Integer,Description="Number of the least frequent homozygous genotype">
##INFO=<ID=MAF,Number=1,Type=Float,Description="Minor allele frequency">
##INFO=<ID=MBF,Number=1,Type=Float,Description="Minor base frequency">
##INFO=<ID=He,Number=1,Type=Float,Description="Nei expected Heterozygosity">
##INFO=<ID=Fis,Number=1,Type=Float,Description="Inbreeding coefficient">
"""
output_vcf = open(args.output, 'w')
found_info = False
wrote_info = False
for line in input_vcf:
    if line[:6] == '#CHROM':
        if not wrote_info:
            output_vcf.write(header)
        output_vcf.write(line)
        break
    if line[:6] == '##INFO':
        output_vcf.write(line)
        found_info |= True
    else:
        assert line[:2] == '##'
        if found_info and not wrote_info:
            output_vcf.write(header)
            wrote_info = True
        output_vcf.write(line)

# prepare computestats
struct = egglib.struct_from_samplesizes([vcf.num_samples], ploidy=2)
cs = egglib.stats.ComputeStats(multi_hits=True, struct=struct)
cs.add_stats('He', 'Fis')

# static objects reused for each site
site = egglib.Site()
frq = egglib.Freq()
alph = egglib.Alphabet(cat='char', name='DNAgap', expl='ACTG-', miss='?', case_insensitive=False)

# process all lines of VCF
cur = None
for ch, pos, na in vcf:
    if ch != cur:
        cur = ch
        print(ch)
    vcf.get_genotypes(dest=site) # get variant as Site object

    # if DNA, replace with custom alphabet treating gap as an allele
    if site.alphabet.name == 'DNA':
        site.from_list(site.as_list(), alphabet=alph)

    # compute stats using ComputeStats
    stats = cs.process_site(site)
    He = stats['He']
    Fis = stats['Fis']

    # set missing data if no valid data / no polymorphism
    if He is None:
        He = '.'
        Fis = '.'
    elif He == 0.0:
        Fis = '.'
    else:
        assert Fis is not None

    # compute allele/genotype frequencies using Freq
    frq.from_site(site, struct=struct)
    nbAll = frq.num_alleles
    nbG = frq.num_genotypes
    All = [frq.allele(i) for i in range(nbAll)]
    genos = [tuple(sorted(frq.genotype(i))) for i in range(nbG)]
    AllCount = [frq.freq_allele(i) for i in range(nbAll)]
    t = sum(AllCount)
    if t > 0: AllFreq = [i/t for i in AllCount]
    else: AllFreq = None
    pG = [frq.freq_genotype(i) for i in range(nbG)]

    # minimum frequency of homozygote genotypes
    if nbAll > 0:
        dG = dict(zip(genos, pG))
        MinHomo = min([dG.get((a, a), 0) for a in All])
    else:
        MinHomo = '.'

    nbGS = sum(pG) # number of genotypes with exploitable data
    pNA = site.num_missing / site.ns if site.ns > 0 else '.' # proportion of missing data

    # site type (SNP and/or INDEL)
    alls_bases = bases.intersection(All)
    alls_diff = set(All).difference(bases)
    SNP = len(alls_bases) > 1
    INDEL = len(alls_diff) > 1 or (len(alls_diff) > 0 and len(alls_bases) > 0)
    if not SNP and not INDEL: assert nbAll == 1

    # MAF
    if nbAll > 1:
        MAF = sorted(AllFreq)[-2]
    else:
        MAF = '.'
    if SNP:
        dA = dict(zip(All, AllCount))
        AllCount_b = [dA[i] for i in alls_bases]
        MBF = sorted(AllCount_b)[-2] / sum(AllCount_b)
    else:
        MBF = '.'

    # insert new INFO in VCF
    line = input_vcf.readline()
    cols = line.split('\t')
    cols[7] += ';'
    cols[7] += 'SNP=1;' if SNP else 'SNP=0;'
    cols[7] += 'INDEL=1;' if INDEL else 'INDEL=0;'
    #if SNP: cols[7] += 'SNP;'
    #if INDEL: cols[7] += 'INDEL;'
    cols[7] += f'nbGS={nbGS};'
    cols[7] += 'pNA=.;' if pNA == '.' else f'pNA={pNA:.6f};'
    cols[7] += f'nbG={nbG};'
    cols[7] += f'nbAll={nbAll};'
    cols[7] += f'All={",".join(All)};'
    cols[7] += 'AllCount=.' if na==0 else f'AllCount={",".join(map(str, AllCount))};'
    cols[7] += 'AllFreq=.' if AllFreq is None else f'AllFreq={",".join(f"{i:.4f}" for i in AllFreq)};'
    cols[7] += 'MinHomo=.;' if MinHomo == '.' else f'MinHomo={MinHomo};'
    cols[7] += 'MAF=.;' if MAF == '.' else f'MAF={MAF:.6f};'
    cols[7] += 'MBF=.;' if MBF == '.' else f'MBF={MBF:.6f};'
    cols[7] += 'He=.;' if He == '.' else f'He={He:.6f};'
    cols[7] += 'Fis=.' if Fis == '.' else f'Fis={Fis:.6f}'

    # export line
    output_vcf.write('\t'.join(cols))
