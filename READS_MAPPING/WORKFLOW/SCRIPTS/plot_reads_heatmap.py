#!/usr/bin/env python

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys

dir=sys.argv[1]
#samples_groups=sys.argv[2]
#zones_groups=sys.argv[3]

sns.set(style="darkgrid")

# Input file
df=pd.read_csv(dir+"/mean_depth_per_zone_per_sample.tsv", sep='\t')

## Heatmap 1: sorted by depth
# Compute columns and rows mean values and reorder the dataframe
df_reindexedRows=df.reindex(df.mean(axis=1).sort_values(ascending=False).index, axis=0)
df_reindexed=df_reindexedRows.reindex(df_reindexedRows.mean(axis=0).sort_values(ascending=False).index, axis=1)


# Plot heatmap
sns.heatmap(df_reindexed, yticklabels=False, xticklabels=False, cmap="BuPu").set(xlabel='Samples', ylabel='Genomic zones')
plt.title("Mean depth per base", fontsize =20)
plt.savefig(dir+"/mean_depth_per_zone_per_sample_heatmap.pdf")

## Heatmap 2: showing groups
#if (len(samples_groups)):




# ajouter éventuel fichier avec des groupes pour grouper sur heatmap
# -> un fichier de groupes de zones (par ex: chr) et un fichier de groupes de génotypes
