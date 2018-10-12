from Bio import SeqIO # https://biopython.org/wiki/SeqIO
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy import stats
import seaborn as sns
sns.set(color_codes=True)

from ggplot import * # http://ggplot.yhathq.com/



# Importing the amino acids
prot_list = []
for record in SeqIO.parse("PP_step1_trn.fasta", "fasta"):
    seq = str(record.seq)
    prot_list.append(seq)
prot_series = pd.Series(prot_list)

# Creating a length series for the amino acids
seq_len = []
for element in prot_list:
    seq_len.append(len(element))
len_series = pd.Series(seq_len)

# Importing the labels
label_list = []
for record in SeqIO.parse("PP_step1_trn copy.fasta", "fasta"):
    seq = str(record.seq)
    label_list.append(seq)
label_series = pd.Series(label_list)
label_series.describe() # Non_Bind: 651; Bind: 536


# Combining length series and labels into one data frame
data = {"Labels" : label_list, "Length" : seq_len}
prot_df = pd.DataFrame(data, columns = ["Labels", "Length"])


# Hypothesis: protein length is in relation to binding properties
ax = sns.boxplot(x="Labels", y="Length", data=prot_df)
plt.show()
plt.savefig("boxplot_len-and-labels")