#python plot_variants_infos_histograms.py --input variants_infos.tsv --output variants_infos_histograms.pdf --bins 1000

import argparse
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy
import seaborn as sns

## Arguments
parser = argparse.ArgumentParser(description="Plot histograms of numeric variables from vcf infos")
parser.add_argument("-i", "--input", dest="tsv", help="path to input tsv file")
parser.add_argument("-o", "--output", dest="pdf", help="name for output pdf file")
args = parser.parse_args()
tsv_input = args.tsv
pdf_output = args.pdf

## Read input file
infos = pd.read_csv(tsv_input, sep="\t")

# Remove contig and pos columns and columns with the same value for all variants
infos.drop(columns = infos.columns[[0,1]], axis = 1, inplace= True)
nunique = infos.nunique()
cols_to_drop = nunique[nunique == 1].index
infos.drop(cols_to_drop, axis=1, inplace= True)

# Only keep numeric columns
numerics = ['float64', 'int64']
others = infos.select_dtypes(exclude=numerics)
infos = infos.select_dtypes(include=numerics)

# Reformat non-numeric columns to try and make it numeric
others_converted_to_float = {}
for column in others:
    list_of_list_of_values = [str(values).split(',') for values in others[column].tolist()]
    list_of_values = ['NaN' if value=="." else value for elem in list_of_list_of_values for value in elem]
    isFloatColumn = True
    list_of_float_values = []
    for value in list_of_values:
        try:
            void = float(value)
            list_of_float_values.append(float(value))
        except ValueError:
            isFloatColumn=False
            break
    if isFloatColumn:
        others_converted_to_float[column] = list_of_float_values


# Are there Qual values less than 1
transform_with_log = True
if (hasattr(infos, 'Qual') and min(infos["Qual"]) < 1):
    transform_with_log = False

# How many SNPs in the file
nb_SNPs = infos.shape[0]
footnote_text = "The histograms are based on data from "+str(nb_SNPs)+" positions. If the input vcf file contained more\nthan 100000 positions, approximately 100000 loci were randomly sampled to plot the histograms from."

# Plot one histogram per column and save it to the output pdf
with PdfPages(pdf_output) as pdf:
    for column in infos:
        plt.figure()
        if (column == "Qual" and transform_with_log):
            fig, ax = plt.subplots()
            histplot = sns.histplot(numpy.log10(infos[column]), kde=True)
            y_min, y_max = plt.gca().get_ylim()
            plt.text(0, 1.1, "Qual thresholds", color="#2F4F4F", transform = ax.transAxes)
            for qual in 10, 20, 30, 40:
                plt.axvline(numpy.log10(qual), color="#2F4F4F", linestyle='dashed', linewidth=1)
                plt.text(1.01*numpy.log10(qual), qual/100*y_max, qual, color="#2F4F4F")
            plt.title("log10(Qual)")
        else:
            histplot = sns.histplot(infos[column], kde=True)
            plt.title(column)
        histplot.set(xlabel=None)
        if (len(histplot.lines) > 0):
            histplot.lines[0].set_color('#e67300')
        plt.figtext(0.5, 0.01, footnote_text, horizontalalignment='center', fontsize=8)
        plt.subplots_adjust(left=None, bottom=0.14, right=None, top=None, wspace=None, hspace=None)
        pdf.savefig(histplot.get_figure())
    for key in others_converted_to_float.keys():
        plt.figure()
        np_arr = numpy.array(others_converted_to_float[key])
        histplot = sns.histplot(np_arr, kde=True)
        plt.title(key)
        if (len(histplot.lines) > 0):
            histplot.lines[0].set_color('#e67300')
        plt.figtext(0.5, 0.01, footnote_text, horizontalalignment='center', fontsize=8)
        plt.subplots_adjust(left=None, bottom=0.14, right=None, top=None, wspace=None, hspace=None)
        pdf.savefig(histplot.get_figure())
