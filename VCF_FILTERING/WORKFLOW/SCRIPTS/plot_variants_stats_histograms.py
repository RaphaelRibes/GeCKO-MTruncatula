#python plot_variants_infos_histograms.py --input variants_infos.tsv --output variants_infos_histograms.pdf --bins 1000

import argparse
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy

## Arguments
parser = argparse.ArgumentParser(description="Plot histograms of numeric variables from vcf infos")
parser.add_argument("-i", "--input", dest="tsv", help="path to input tsv file")
parser.add_argument("-o", "--output", dest="pdf", help="name for output pdf file")
args = parser.parse_args()
tsv_input = args.tsv
pdf_output = args.pdf


## Read input file
infos = pd.read_csv(tsv_input, sep="\t")

# Remove contig and pos columns
infos.drop(columns = infos.columns[[0,1]], axis = 1, inplace= True)

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


# Choose the bins number depending on the number of variants
nb_rows = len(infos.index)
if (nb_rows < 200):
    bins = 20
elif (nb_rows < 400):
    bins = 50
else:
    bins = 100


# Are there Qual values less than 1
transform_with_log = True
if (hasattr(infos, 'Qual') and min(infos["Qual"]) < 1):
    transform_with_log = False


# Plot one histogram per column and save it to the output pdf
with PdfPages(pdf_output) as pdf:
    for column in infos:
        if (column == "Qual" and transform_with_log):
            fig, ax = plt.subplots()
            plt.hist(numpy.log10(infos["Qual"]), bins=bins)
            y_min, y_max = plt.gca().get_ylim()
            plt.text(0, 1.1, "Qual thresholds", color="#2F4F4F", transform = ax.transAxes)
            for i in 10, 20, 30, 40:
              plt.axvline(numpy.log10(i), color="#2F4F4F", linestyle='dashed', linewidth=1)
              plt.text(1.01*numpy.log10(i), i/100*y_max, i, color="#2F4F4F")
            plt.title("log10(Qual)")
            pdf.savefig(fig)
        else:
            fig, ax = plt.subplots()
            infos.hist(column, ax=ax, bins=bins)
            pdf.savefig(fig)
    for key in others_converted_to_float.keys():
        fig, ax = plt.subplots()
        np_arr = numpy.array(others_converted_to_float[key])
        plt.hist(np_arr, bins=bins)
        plt.title(key)
        pdf.savefig(fig)
