#!/usr/bin/python3
from scipy.stats import wilcoxon
from scipy.stats import norm
from scipy.stats import kruskal
from statistics import mean
#import scikit_posthocs as sp
from io import StringIO
import sys
import re
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

# Control if we generate CSV or a text format
generateCSV = 1

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
#vs_refmsacons:refiner
sel_cols = [col for col in df.columns if 'vs_refmsacons' in col and 'gput' not in col]
data = []
pdat = []

dft = df[df['replicate_group'] == 1]
for gput in dft['vs_refmsacons:gput'].values.tolist():
    for replicate in range(1,11):
        dfl=df[(df['replicate_group'] == replicate) & (df['vs_refmsacons:gput'] == gput)]
        row = []
        for col1_idx in range(0,len(sel_cols)):
            vals = dfl.loc[:,[sel_cols[col1_idx]]].values.tolist()
            row.append(vals[0][0])
        #print("gput=" + str(gput) + " rep="+str(replicate)+ ": ", end='' )
        #print(row)
        data.append(row)

#for replicate in range(0,10):
#    row = []
#    for col1_idx in range(0,len(sel_cols)):
#        dfl = df[df['replicate_group'] == replicate+1] \
#               .loc[:,[ 'vs_refmsacons:gput',sel_cols[col1_idx]]]
#        gput = dfl['vs_refmsacons:gput'].values.tolist()
#        alg = dfl[sel_cols[col1_idx]].values.tolist()
#        alg_area = np.trapz(alg,gput)
#        row.append(alg_area)
#        pdat.append([alg_area,sel_cols[col1_idx][9:]])
#    data.append(row)
## Data: rows=replicates cols=methods value=area

nDF = pd.DataFrame(data, columns=sel_cols)


#paper-data/DNATransTree-1-Tigger1-R3S-eval/replicates.csv
#paper-data/DNATransTree-1-Tigger1-R3S-gput100-mfl2-eval/replicates.csv
mats = re.match("^.*\/(DNATransTree|LINETree)-(\d+)-(\S+)-R3S-eval\/.*",sys.argv[1])
if ( not mats ):
    mats = re.match("^.*\/(DNATransTree|LINETree)-(\d+)-(\S+)-R3S-gput(\d+)-mfl2-eval\/.*",sys.argv[1])
if ( not mats ):
    print("Could not find " + sys.argv[1])
    exit(1)

print()
if ( mats.lastindex > 3 ):
    print("Derived Consensus Statistics for " + mats.group(3) + " Sequence Fragmentation Analysis" )
else:
    print("Derived Consensus Statistics for " + mats.group(3) + " Sequence Divergence Analysis" )
print()
if ( mats.group(1) == "DNATransTree"):
    print("Tree Simulation: DNA Transposon Tree #" + mats.group(2))
else:
    print("Tree Simulation: LINE Tree #" + mats.group(2))
if ( mats.lastindex > 3 ):
    if ( mats.group(4) == "100" ):
        print("Fragmentation Simulation, Low Divergence Sequences [gput100]")
    if ( mats.group(4) == "1500" ):
        print("Fragmentation Simulation, Medium Divergence Sequences [gput1500]")
    if ( mats.group(4) == "3000" ):
        print("Fragmentation Simulation, High Divergence Sequences [gput3000]")
print("Data File: " + sys.argv[1])
print()


#pDF = pd.DataFrame(pdat, columns=['area','alg'])
#plt.figure(figsize=(15,8))
#sns.set(style="whitegrid")
#ax = sns.boxplot(x="alg",y="area", data=pDF, showfliers = False)
#ax = sns.swarmplot(x="alg",y="area", data=pDF, color=".25")
#plt.savefig("stats.svg")
#plt.savefig("stats.png")
#plt.savefig("stats.pdf")

# Calculate Kruskal-Wallis H-test
kDat = nDF.T.values.tolist()
h,p = kruskal(*kDat)
print("Kruskal-Wallis H-test")
if ( generateCSV ):
    print(",H =,"+str(h)+",p =,"+str(p))
else:
    print("    H="+str(h)+", p="+str(p))
print()

# Wilcoxon Signed Rank Test
print ("Wilcoxon signed rank test: mean_diff [p-val]")
if ( not generateCSV ):
    print ("    " + "{:19s}".format(""),end='')
    for col2 in range(0,len(nDF.columns)):
        print ("{:17s}".format(nDF.columns[col2].replace("vs_refmsacons:","")), end='')
    print()
    for col1 in range(0,len(nDF.columns)):
        print ("    " + "{:17s}".format(nDF.columns[col1].replace("vs_refmsacons:","")), end='')
        for col2 in range(0,len(nDF.columns)):
            if ( col2 != col1 ):
                data1 = nDF[nDF.columns[col1]].tolist()
                data2 = nDF[nDF.columns[col2]].tolist()
                minW, p = wilcoxon(data1, data2)
                # The area being evaluated here is being minimized ( the values
                # are the distance from the ideal alignment score ).
                print("{:6.2f}".format(mean(data1)-mean(data2)) + " " + "[{:5.1e}]".format(p) + " ",end="")
            else:
                print("{:17s}".format(""), end='')
        print('')
    print()
else:
    # CSV format
    for col2 in range(0,len(nDF.columns)):
        print ("," + nDF.columns[col2].replace("vs_refmsacons:",""), end='')
    print ('')
    for col1 in range(0,len(nDF.columns)):
        print (nDF.columns[col1].replace("vs_refmsacons:",""), end='')
        for col2 in range(0,len(nDF.columns)):
            if ( col2 != col1 ):
                data1 = nDF[nDF.columns[col1]].tolist()
                data2 = nDF[nDF.columns[col2]].tolist()
                minW, p = wilcoxon(data1, data2)
                print(",{:6.2f}".format(mean(data1)-mean(data2)) + " " + "[{:5.1e}]".format(p),end="")
            else:
                print(",", end='')
        print('')

