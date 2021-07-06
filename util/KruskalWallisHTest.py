#!/usr/bin/python3
from scipy.stats import wilcoxon
from scipy.stats import norm
from scipy.stats import kruskal
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
with open("/u1/home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/paper-data/DNATransTree-1-Tigger1-R3S-eval-new2/replicates.csv") as infile:
#with open("/u1/home/rhubley/projects/DNAMultipleAlignment/TEForwardEvolve/paper-data/LINETree-1-L2-R3S-eval-new/replicates.csv") as infile:
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
#sel_cols = [col for col in df.columns if 'QScore_Q' in col]
sel_cols = [col for col in df.columns if 'QScore_Q' in col and 'gput' not in col]
data = []
pdat = []
for replicate in range(0,10):
    row = []
    for col1_idx in range(0,len(sel_cols)):
        dfl = df[df['replicate_group'] == replicate+1] \
               .loc[:,[ 'QScore_Q:gput',sel_cols[col1_idx]]]
        gput = dfl['QScore_Q:gput'].values.tolist()
        alg = dfl[sel_cols[col1_idx]].values.tolist()
        alg_area = np.trapz(alg,gput)
        row.append(alg_area)
        pdat.append([alg_area,sel_cols[col1_idx][9:]])
        #print("" + sel_cols[col1_idx] + " = " + str(alg_area))
    data.append(row)



#print(pdat)
nDF = pd.DataFrame(data, columns=sel_cols)
print(nDF)
pDF = pd.DataFrame(pdat, columns=['area','alg'])
plt.figure(figsize=(15,8))
sns.set(style="whitegrid")
ax = sns.boxplot(x="alg",y="area", data=pDF, showfliers = False)
ax = sns.swarmplot(x="alg",y="area", data=pDF, color=".25")
plt.savefig("test.svg")
plt.savefig("test.png")
plt.savefig("test.pdf")

# Calculate Kruskal-Wallis
#for col1_idx in range(0,len(sel_cols)):
foo=  nDF.T.values.tolist()
#print(foo)
h,p = kruskal(*foo)
print("H="+str(h)+", p="+str(p))

exit(0)


# Calculate area using trapezoid method (since we have non-uniform spacing)
#  result = np.trapz(y,x)
#    - One result per algorithm per replicate
#               alg1 alg2
#        rep1
#        rep2
#
# Create box plot with points
#    sns.set(style="whitegrid")
#    tips = sns.load_dataset("tips")
#
#    ax = sns.boxplot(x="day", y="total_bill", data=tips, showfliers = False)
#    ax = sns.swarmplot(x="day", y="total_bill", data=tips, color=".25")
#
#    plt.show()
#    plt.savefig("test.svg")
#    plt.savefig("test.png")
#    plt.savefig("test.pdf")
#
# Calculate Kruskal-Wallis
#   h,p = kruskal(alg1_array,alg2_array,...)
#

# Select data from the vs_refmsacons analysis.
# Generate all pairs of treatments ( algorithms ) and perform Wilcoxon
# signed-rank test on each pair.
sel_cols = [col for col in df.columns if 'vs_refmsacons' in col and 'gput' not in col]
for gput in [100,250,500,750,1000,1250,1500,1750,2000,2250,2500,2750,3000,3250,3500,4000,5000,6000]:
    print("Wilcoxon Signed-Rank Test (significance=0.05, 2-tailed hypothesis): VS Reference MSA Consensus (NW score fraction) gput=" + str(gput))
    for col1_idx in range(0,len(sel_cols)-1):
        print(","+sel_cols[col1_idx],end="")
    print("")
    for col1_idx in range(0,len(sel_cols)):
        print (sel_cols[col1_idx], end="")
        for col2_idx in range(0,col1_idx):
            dfl = df[df['AMA_predictive_value:gput'] == gput] \
                   .loc[:,[ sel_cols[col1_idx],sel_cols[col2_idx]]]
            data1 = dfl[sel_cols[col1_idx]].values.tolist()
            data2 = dfl[sel_cols[col2_idx]].values.tolist()
            minW, p = wilcoxon(data1, data2)
            print(", %.3f"%p,end="")
        print("")
    print("\n")


# Z-score calculation for reference
#minW, p = wilcoxon(data1, data2)
#n = len(data1)
#z=norm.isf(p/2)
#if minW-((n*(n+1))/4) < 0:
#    z *= -1
#print('n=%.3f, min(W-,W+)=%.3f, z=%.3f, p=%.3f' % (n, minW, z, p))
#
## interpret
#alpha = 0.05
#if p > alpha:
#    print('Same distribution (fail to reject H0)')
#else:
#    print('Different distribution (reject H0)')
#
