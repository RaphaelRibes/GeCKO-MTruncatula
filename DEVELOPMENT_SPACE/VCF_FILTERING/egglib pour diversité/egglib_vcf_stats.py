### egglib_vcf_stats.py ###

import egglib
import argparse

# variables
parser = argparse.ArgumentParser(description="")
parser.add_argument("-b", "--bed_file", dest="bed_file", help="calculate population genetics parameters with Egglib from vcf file")
parser.add_argument("-v", "--vcf_file", dest="vcf_file", help="path to vcf file Filtered")
parser.add_argument("-l", "--labels_file", dest="labels_file", help="path to file containing labels groups")
args = parser.parse_args()
bed_file = args.bed_file
vcf_file = args.vcf_file
labels_file = args.labels_file

stats_tot = 'nseff', 'S', 'thetaW', 'Pi', 'D', 'He', 'FstWC', 'Aing'
stats_pop = 'nseff', 'S', 'thetaW', 'Pi', 'D', 'He', 'Fis'
defaults = {'nseff': 'NA', 'S': 0, 'thetaW': 0, 'Pi': 0, 'D': 'NA', 'He': 0, 'FstWC': 'NA', 'Fis': 'NA', 'Aing': 0} # default values if nseff=0
undefined = {'D', 'Fis', 'FstWC'} # undefined if S=0
outfile1 = 'egglib_stats.txt'
outfile2 = 'egglib_stats_pairwise.txt'
outfile3 = 'egglib_outliers.txt'



vcf = egglib.io.VcfParser(vcf_file)
samples = [vcf.get_sample(i) for i in range(vcf.num_samples)]

# import labels and create master structure
labels = []
indivs = []
mapping1 = {}
mapping2 = {}
with open(labels_file, "r") as f:
    for row in f:
        pop, name, label = row.split()
        labels.append(pop)
        indivs.append(name)
        if pop not in mapping1: mapping1[pop] = []
        mapping1[pop].append(name)
        mapping2[name] = label

# import BED positions
bed = {}
with open(bed_file) as f:
    for line in f:
        ch, start, end = line.split()
        if ch not in bed: bed[ch] = []
        bed[ch].append((int(start), int(end)))

# create structure objects
pops = list(set(labels))
idx_pops = {}
for p in pops:
    idx_pops[p] = {}
    for n in mapping1[p]:
        rk = samples.index(mapping2[n])
        idx_pops[p][n] = (2*rk, 2*rk+1)
d = {None: {p: idx_pops[p] for p in pops}}
struct = egglib.struct_from_dict(d, None)

# create master CS
cs = egglib.stats.ComputeStats(struct=struct, multi=True)
cs.add_stats(*stats_tot)

# create one CS per pair of populations
cs_pops = {}
cs_pairs = {}
for i in range(len(pops)):
    d = {None: {pops[i]: idx_pops[pops[i]]}}
    cs_pops[pops[i]] = egglib.stats.ComputeStats(
                    struct=egglib.struct_from_dict(d, None), multi=True)
    cs_pops[pops[i]].add_stats(*stats_pop)

    for j in range(i+1, len(pops)):
        d = {None:
                {pops[i]: idx_pops[pops[i]],
                 pops[j]: idx_pops[pops[j]]}}
        cs_pairs[(pops[i], pops[j])] = egglib.stats.ComputeStats(
                    struct=egglib.struct_from_dict(d, None), multi=True)
        cs_pairs[(pops[i], pops[j])].add_stats('FstWC')

# iterate over VCF positions
cur_c = None # cur_c and cur_i identify the current window / contig
cur_i = None
ignored = 0
results = {} # per window (total and per-pop) stats index by window identifies

print('Running', end='', flush=True)
c = 0
for ch, pos, nall in vcf:

    # if position out of current window (contig)
    if ch != cur_c or pos < bed[ch][cur_i][0] or pos >= bed[ch][cur_i][1]:

        # next window on same chromosome
        if ch == cur_c and cur_i < len(bed[ch])-1 and pos >= bed[ch][cur_i+1][0] and pos < bed[ch][cur_i+1][1]:
            next_c = cur_c
            next_i += 1

        # any other window
        else:
            for i, (start, end) in enumerate(bed[ch]):
                if pos >= start and pos < end:
                    next_c = ch
                    next_i = i
                    break
            else:
                ignored += 1
                continue

        # store results
        if cur_c != None:
            results[(cur_c, cur_i)] = (cs.results(),
                            {i: cs_pops[i].results() for i in cs_pops})
        if next_c != cur_c:
            print('\n'+next_c, end=' ', flush=True)
        else:
            if cur_i//1000 > c:
                print('.', end='', flush=True)
                c = cur_i//1000
        cur_c = next_c
        cur_i = next_i

    # analyze site
    cs.process_site(vcf.get_genotypes())
    for i in cs_pops.values(): i.process_site(vcf.get_genotypes())
    for i in cs_pairs.values(): i.process_site(vcf.get_genotypes())

# store results of last window
results[(cur_c, cur_i)] = cs.results(), {i: cs_pops[i].results() for i in cs_pops}

# report results per window
with open(outfile1, 'w') as f:
    header = ['chrom', 'start', 'end'] + list(stats_tot)
    for p in pops: header.extend([p + '_' + s for s in stats_pop])
    f.write('\t'.join(header) + '\n')
    for c in bed:
        for i, (start, end) in enumerate(bed[c]):
            row = [c, start, end]
            if (c, i) in results:
                stats, stats_pops = results[(c, i)]
                if stats['nseff'] is None:
                    for k in stats: stats[k] = defaults[k]
                    stats['nseff'] = 0
                elif stats['S'] == 0:
                    for k in stats:
                        if k in undefined:
                            stats[k] = 'NA'
                for s in stats_tot: row.append(stats[s])
                for pop in pops:
                    if stats_pops[pop]['S'] == 0:
                        for k in stats_pops[pop]:
                            if k in undefined: stats_pops[pop][k] = 'NA'
                    for s in stats_pop:
                        row.append(stats_pops[pop][s])
            else:
                for s in stats_tot: row.append(defaults[s])
                for p in pops:
                    for s in stats_pop: row.append(defaults[s])

            f.write('\t'.join(map(str, row)) + '\n')

# report pairwise FstWC
with open(outfile2, 'w') as f:
    f.write('Pop1\tpop2\tFstWC\n')
    for (pop1, pop2), cs in cs_pairs.items():
        f.write(f'{pop1}\t{pop2}\t{cs.results()["FstWC"]:.5f}\n')

# report files outside BED regions (contigs)
with open(outfile3, 'w') as f:
    f.write(f'VCF sites outside BED regions: {ignored}\n')

print('\nDone')
