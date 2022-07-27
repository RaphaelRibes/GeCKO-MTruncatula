import pandas as pd
import re

input = "test70.vcf"
output = "test70_order.vcf"


with open(input, 'r') as f2:
    for line in f2:
        if re.match('#CHROM', line):
            header = list(line.split("\t"))
            break

with open(input, 'r') as f2:
    with open(output, 'w') as f3:
        df = pd.DataFrame(columns=header)
        i = 0
        for line in f2:
            if line.startswith("#"):
                f3.write(line)
            else:
                df.loc[i] = list(line.split("\t"))
                i = i+1
        print(df)
        df2 = df.sort_values(by=['#CHROM','POS'], ascending = (True, True))
        print(df2)
        df2.reset_index(inplace=True, drop=True)
        print(df2)
        for row in range(len(df2)):
            f3.write('\t'.join(df2.loc[row]))
