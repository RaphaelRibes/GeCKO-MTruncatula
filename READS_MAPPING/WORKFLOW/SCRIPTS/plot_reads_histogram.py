#!/usr/bin/env python

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys


dir=sys.argv[1]

sns.set(style="darkgrid")


df=pd.read_csv(dir+"/summary_flagstat.tsv", sep='\t')

fig, axs = plt.subplots(2, 1)
sns.histplot(data=df, x="Reads_mapped", color="skyblue", ax=axs[0])
sns.histplot(data=df, x="%_reads_mapped", color="teal", ax=axs[1])
plt.show()
fig.tight_layout(pad=1.5)
plt.savefig(dir+"/reads_histograms.pdf")
