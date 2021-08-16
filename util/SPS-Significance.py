#!/usr/bin/python3
from scipy.stats import wilcoxon
from scipy.stats import norm
from scipy.stats import kruskal
from statistics import mean
#import scikit_posthocs as sp
import sys
from io import StringIO
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# Read in replicate CSV with multiple tables seperated by empty lines and
# convert to a single table suitable for pandas.  The single table joins the
# multiple tables by prefixing a column titled "replicate_group" which holds
# a 1-based index of the individual tables.
replicate_idx = 0
file_data = ""
with open(sys.argv[1]) as infile:
    for line in infile:
        if line.startswith("AMA"):
            replicate_idx += 1
            line = "replicate_group," + line
            if ( replicate_idx == 1 ):
                file_data += line
        elif not line.strip():
            continue
        else:
            line = str(replicate_idx) + "," + line
            file_data += line

CSVDATA = StringIO(file_data)
df = pd.read_csv(CSVDATA)

# Select data from the Qscore tool, specifically the SPS data (QScore_Q).
sel_cols = [col for col in df.columns if 'QScore_Q' in col and 'gput' not in col]
data = []
pdat = []

dft = df[df['replicate_group'] == 1]
for gput in dft['QScore_Q:gput'].values.tolist():
    for replicate in range(1,11):
        dfl=df[(df['replicate_group'] == replicate) & (df['QScore_Q:gput'] == gput)]
        row = []
        for col1_idx in range(0,len(sel_cols)):
            vals = dfl.loc[:,[sel_cols[col1_idx]]].values.tolist()
            row.append(vals[0][0])
        #print("gput=" + str(gput) + " rep="+str(replicate)+ ": ", end='' )
        #print(row)
        data.append(row)

#                   tool1  tool2  tool3 ...
# gput=100 rep=1:  [0.998, 0.993, 0.998, 0.997, 0.998, 0.998, 0.912]
# gput=100 rep=2:  [0.999, 0.980, 0.998, 0.998, 0.997, 0.998, 0.912]
# gput=100 rep=3:  [0.998, 0.961, 0.998, 0.998, 0.996, 0.998, 0.921]
# gput=100 rep=4:  [0.997, 0.995, 0.998, 0.997, 0.996, 0.998, 0.927]
# gput=100 rep=5:  [0.997, 0.975, 0.997, 0.996, 0.998, 0.998, 0.899]
# gput=100 rep=6:  [0.998, 0.991, 0.998, 0.997, 0.995, 0.998, 0.888]
# gput=100 rep=7:  [0.998, 0.990, 0.998, 0.996, 0.996, 0.998, 0.923]
# gput=100 rep=8:  [0.996, 0.989, 0.997, 0.997, 0.997, 0.998, 0.919]
# gput=100 rep=9:  [0.998, 0.985, 0.998, 0.998, 0.998, 0.998, 0.896]
# gput=100 rep=10: [0.996, 0.993, 0.998, 0.997, 0.998, 0.998, 0.893]
# gput=250 rep=1:  [0.994, 0.970, 0.992, 0.987, 0.992, 0.993, 0.832]
# ...


#for replicate in range(0,10):
#    row = []
#    for col1_idx in range(0,len(sel_cols)):
#        dfl = df[df['replicate_group'] == replicate+1] \
#               .loc[:,[ 'QScore_Q:gput',sel_cols[col1_idx]]]
#        # gputs
#        gputs = dfl['QScore_Q:gput'].values.tolist()
#        # sps
#        algs = dfl[sel_cols[col1_idx]].values.tolist()
#        #alg_area = np.trapz(algs,gputs)
#        #row.append(alg_area)
#        #pdat.append([alg_area,sel_cols[col1_idx][9:]])
#        for alg in algs:
#           row.append(alg)
#    data.append(row)
## Data: rows=replicates-gput cols=methods value=sps

#print(str(data))
#exit(0)
nDF = pd.DataFrame(data, columns=sel_cols)

#pDF = pd.DataFrame(pdat, columns=['area','alg'])
#plt.figure(figsize=(15,8))
#sns.set(style="whitegrid")
#ax = sns.boxplot(x="alg",y="area", data=pDF, showfliers = False)
#ax = sns.swarmplot(x="alg",y="area", data=pDF, color=".25")
#plt.savefig("stats.svg")
#plt.savefig("stats.png")
#plt.savefig("stats.pdf")

print()
print ("SPS Statistics for " + sys.argv[1])
print ("  - Calculated on all replicate-parameter values")
print()

# Calculate Kruskal-Wallis H-test
kDat =  nDF.T.values.tolist()
h,p = kruskal(*kDat)
print("Kruskal-Wallis H-test")
print("    H="+str(h)+", p="+str(p))
print()


print ("Wilcoxon signed rank test: mean_diff [p-val]")
print ("    " + "{:19}".format(""),end='')
for col2 in range(0,len(nDF.columns)):
    print ("{:17}".format(nDF.columns[col2].replace("QScore_Q:","").replace("-padded","")+" "), end='')
print()
for col1 in range(0,len(nDF.columns)):
    print ("    " + "{:17}".format(nDF.columns[col1].replace("QScore_Q:","").replace("-padded","")), end='')
    for col2 in range(0,len(nDF.columns)):
        if ( col2 != col1 ):
            data1 = nDF[nDF.columns[col1]].tolist()
            data2 = nDF[nDF.columns[col2]].tolist()
            minW, p = wilcoxon(data1, data2)
            print("{:6.2f}".format(mean(data1)-mean(data2)) + " " + "[{:5.1e}]".format(p) + " ",end="")
        else:
            print("{:17}".format(""), end='')
    print('')
print()

# CSV
#print ("Wilcoxon signed rank test:")
#for col2 in range(1,len(nDF.columns)):
#    print ("," + nDF.columns[col2], end='')
#print ('')
#for col1 in range(0,len(nDF.columns)):
#    print (nDF.columns[col1], end='')
#    for col2 in range(0,len(nDF.columns)):
#        if ( col2 > col1 ):
#            data1 = nDF[nDF.columns[col1]].tolist()
#            data2 = nDF[nDF.columns[col2]].tolist()
#            minW, p = wilcoxon(data1, data2)
#            print(",%.3f"%p,end="")
#        else:
#            print(",", end='')
#    print('')
#

