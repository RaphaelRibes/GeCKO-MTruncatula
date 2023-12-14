import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import argparse

## Arguments
parser = argparse.ArgumentParser(description="Plot boxplot of DP values from vcf genotypes fields")
parser.add_argument("-i", "--input-DP", dest="DP", help="path to DP input tsv file")
parser.add_argument("-j", "--input-GT", dest="GT", help="path to GT input tsv file")
parser.add_argument("-o", "--output", dest="pdf", help="name for output pdf file")
args = parser.parse_args()
DP_input = args.DP
GT_input = args.GT
pdf_output = args.pdf

## Read input file
df=pd.read_csv(DP_input, header=None, sep="\t")
GT=pd.read_csv(GT_input, header=None, sep="\t")
df = df.mask(GT=="./.", GT).mask(GT==".|.", GT).mask(GT==".", GT)
df.drop(df.columns[[0, 1]], axis = 1, inplace = True)

## Make it numeric
df = df.apply(pd.to_numeric, errors='coerce')
df['mean'] = df.mean(axis=1)
m = df['mean'].mean()

# How many SNPs in the file
nb_SNPs = df.shape[0]
footnote_text = "The boxplot is based on data from "+str(nb_SNPs)+" positions. If the input vcf file contained more\nthan 100000 positions, approximately 100000 loci were randomly sampled to plot the boxplot from."

# NA % in the dataframe
na_p=round(100*df.drop('mean', axis=1).isnull().sum().sum()/df.drop('mean', axis=1).size,2)
na_text="NA percent: "+str(na_p)+" %"

# Boxplot
fig, ax = plt.subplots()
boxplot = ax.boxplot(df['mean'])
ax.set_ylim(bottom=0)

for line in boxplot['medians']:
    x, y = line.get_xydata()[1]
    text = ' μ={:.2f}'.format(m)
    ax.annotate(text, xy=(x, y))

plt.title("Mean DP per genotype (SNP X sample)")
plt.tick_params(labelbottom = False, bottom = False)
plt.figtext(0.5, 0.01, footnote_text, horizontalalignment='center', fontsize=8)
plt.subplots_adjust(left=None, bottom=0.14, right=None, top=None, wspace=None, hspace=None)
plt.text(0.6, 0.8*max(df['mean']), na_text, bbox=dict(boxstyle='round', facecolor='wheat'))
plt.savefig(pdf_output)
