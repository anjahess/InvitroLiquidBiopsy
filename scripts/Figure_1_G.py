"""

Scripts and source data to reproduce figures of the manuscript:
Hess, Anja et al. 2025: Non-disruptive monitoring of cellular dynamics with cell-free DNA methylation

Max Planck Institute for Molecular Genetics, Berlin, Germany

Author: Anja Hess
Date: 2025 APRIL 07

"""
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

####################################################################################
# FIG 1G - Delta plot for WGBS
####################################################################################
source = "./sourcedata/FIG_1/1G/1d_plus/1kb_tiles_averages_shared_mean.tsv"

df = pd.read_table(source, sep='\t', index_col=0)
df2 = df.pivot(index="day", columns="condition", values="mean_meth")
df2["delta"] = (df2["cfDNA-1"]-df2["gDNA"]).abs()
mean_delta = df2["delta"].mean()
sns.lineplot(data=df, x='day', y='mean_meth', hue="condition",
             palette=["#85ada3", "#2d435b"])
plt.title(f"Av Methylation Delta g/cf = {mean_delta}")
plt.ylim(0, 1)
plt.ylabel("Average CpG Methylation")
plt.savefig(source.replace(".tsv", f".pdf"))
plt.close()

# END OF SCRIPT