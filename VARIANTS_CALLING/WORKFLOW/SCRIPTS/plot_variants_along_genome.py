import pandas as pd
import matplotlib.pyplot as plt
import argparse
import sys

## Arguments
parser = argparse.ArgumentParser(description="Plot the detected SNPs along the reference genome")
parser.add_argument("-i", "--contigs-lengths", dest="lengths", help="path to DP input tsv file")
parser.add_argument("-j", "--snp-pos", dest="pos", help="path to GT input tsv file")
parser.add_argument("-o", "--output", dest="pdf", help="name for output pdf file")
args = parser.parse_args()
lengths_input = args.lengths
pos_input = args.pos
pdf_output = args.pdf


contigs_lengths = pd.read_csv(lengths_input, sep = "\t")
snp_pos = pd.read_csv(pos_input, sep = "\t")

nb_contigs = len(contigs_lengths)

# Don't do the plot if there are too many contigs
if nb_contigs > 50:
    fig, ax = plt.subplots()
    ax.text(0.5, 0.5, '[Too many contigs to show]', horizontalalignment='center', verticalalignment='center')
    plt.tick_params(left = False, labelleft = False , labelbottom = False, bottom = False)
    plt.savefig(pdf_output)
    sys.exit(0)

# Plot contigs
fig, ax = plt.subplots()
for i in range(0, nb_contigs):
    ax.plot([1, contigs_lengths["length"][i]], [i+1, i+1], color = "black")

plt.yticks(range(1, nb_contigs+1), contigs_lengths["contig"])

# Add SNPs
for i, contig in enumerate(contigs_lengths["contig"]):
    snp_contig = snp_pos.loc[snp_pos["contig"] == contig]
    for pos in snp_contig["pos"]:
        ax.plot([pos, pos], [i+1-0.3, i+1+0.3], color = "#990033", linewidth = 0.8, alpha = 0.6)

plt.title("Representation of the detected variants along the reference chromosomes", fontsize = 10)
plt.xlabel("Position (pb)")

plt.savefig(pdf_output)
