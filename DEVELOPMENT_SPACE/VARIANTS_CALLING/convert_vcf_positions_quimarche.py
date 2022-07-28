import os
import re
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description="convert positions of the polymorphic sites on the complete genomic reference")
parser.add_argument("-v", "--vcf", dest="vcf", help="path to input vcf raw file ")
parser.add_argument("-c", "--vcf_converted", dest="vcf_converted", help="path to vcf file converted")
args = parser.parse_args()
vcf = args.vcf
vcf_converted = args.vcf_converted



flag = True


chr_size_df = pd.read_csv('reference_chr_size.txt', sep="\t",header=None)
chr_size = list(chr_size_df.to_records(index=False))

with open(vcf, 'r', errors="replace") as f1:
    with open(vcf_converted, 'w') as f2:
        for line in f1:
            if flag and re.match('##source', line):
                flag = False
                for k, L in chr_size:
                    f2.write('##contig=<ID='+str(k)+',length='+str(L)+'>\n')
            if not re.match('##contig', line): f2.write(line)
            if re.match('#CHROM', line):
                header = list(line.strip().split("\t"))
                break
        conv_pos = []
        for line in f1:
            bits = line.split()
            ch, start, end = bits[0].rsplit('_', 2)
            start = int(start)
            end = int(end)
            pos = int(start) + int(bits[1]) - 1
            bits[0] = ch
            bits[1] = str(pos)
            conv_pos.append(tuple(bits))
        conv_pos_df = pd.DataFrame(conv_pos, columns = header)
        conv_pos_order_df = conv_pos_df.sort_values(by=['#CHROM','POS'], ascending = (True, True))
        conv_pos_order_df.reset_index(inplace=True, drop=True)
        for row in range(len(conv_pos_order_df)):
            f2.write('\t'.join(conv_pos_order_df.loc[row])+"\n")
