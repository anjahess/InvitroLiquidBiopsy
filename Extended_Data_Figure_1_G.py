"""

Scripts and source data to reproduce figures of the manuscript:
Hess, Anja et al. 2025: Non-disruptive monitoring of cellular dynamics with cell-free DNA methylation

Max Planck Institute for Molecular Genetics, Berlin, Germany

Author: Anja Hess
Date: 2025 APRIL 07

"""
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

####################################################################################
# CpG coverage plotting
####################################################################################
df = pd.read_csv("./sourcedata/FIG_1EXT/E1G/FigE1G_WGBS-filter-stats.csv",
                 index_col=0)
print(df)
plt.figure(figsize=(2,3))
df["specimen"] = df["sample"].str.split('_', expand=True)[0]
sns.barplot(data=df, x="specimen", y="n_cpg", palette=["#5f9ea0","#2c4e60", "salmon"],  edgecolor="black")
sns.swarmplot(data=df, x="specimen", y="n_cpg", palette=["#5f9ea0","#2c4e60","salmon"], linewidth=0.5,edgecolor="black", s=4)
sns.barplot(data=df, x="specimen", y="n_cpg_filtered", palette=["#5f9ea0","#2c4e60", "salmon"], edgecolor="black")
sns.swarmplot(data=df, x="specimen", y="n_cpg_filtered", palette=["#5f9ea0","#2c4e60","salmon"], linewidth=0.5,edgecolor="black", s=4)
plt.ylim(0, 2.5e7)
plt.tight_layout()
plt.savefig("./sourcedata/FIG_1EXT/E1G/filter_stats.pdf", bbox_inches="tight")
plt.close()

# END OF SCRIPT