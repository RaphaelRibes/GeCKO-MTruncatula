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
    list_of_list_of_values = [values.split(',') for values in others[column].tolist()]
    list_of_values = [value for elem in list_of_list_of_values for value in elem]
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

# Plot one histogram per column and save it to the output pdf
with PdfPages(pdf_output) as pdf:
    for column in infos:
        fig, ax = plt.subplots()
        infos.hist(column, ax=ax, bins=bins)
        pdf.savefig(fig)
    for key in others_converted_to_float.keys():
        fig, ax = plt.subplots()
        np_arr = numpy.array(others_converted_to_float[key])
        plt.hist(np_arr, bins=bins)
        plt.title(key)
        pdf.savefig(fig)
