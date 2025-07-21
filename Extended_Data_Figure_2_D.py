"""

Scripts and source data to reproduce figures of the manuscript:
Hess, Anja et al. 2025: Non-disruptive monitoring of cellular dynamics with cell-free DNA methylation

Max Planck Institute for Molecular Genetics, Berlin, Germany

Author: Anja Hess
Date: 2025 MAY 22

"""
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from plots import *

results_file = "./sourcedata/FIG_2EXT/E2D/Bioreactor-hepatocytes_ONT-filter-stats.csv"
plot_dir = "./sourcedata/FIG_2EXT/E2D/"


####################################################################################
# FILTER STATS PER REPLICATE
####################################################################################
df = pd.read_csv(results_file)
df = df.replace(pd.NA, np.nan)
df["input_ng_normalized_bef-adapter-lig"] = (
    df["input_ng_normalized_bef-adapter-lig"].astype(float))
col_val = "biol_rep_group" #"experiment"
####################################################################################
# FILTER STATS PER EXPERIMENR
####################################################################################
fig, ax = plt.subplots(2,2)
palette = ["#5f9ea0",  "#2c4e60", "salmon", "blue"]
order = ["gDNA","cell-free DNA","reference"]
sns.relplot(
    data=df, x="filtering", y="n_cpg_filtered",
    hue="specimen", hue_order=order, style="specimen",
    style_order=order,
    col=col_val, kind="line", palette=palette,
    linewidth=0.5, height=5, aspect=.75, facet_kws=dict(sharex=False),
    ax=ax[0]
)
plt.ylim(0, None)

plt.savefig(results_file.replace(".csv", f"_by-{col_val}_line.pdf"),
            bbox_inches="tight")
plt.close()

####################################################################################
# BAR PLOTS FOR NANOGRAMS INPUT
####################################################################################
df = df[df["specimen"] != "reference"]
sns.catplot(
    data=df, palette=["#85ada3", "#2d435b"],
    x="input_ng_normalized_bef-adapter-lig", y="specimen", col=col_val,
    kind="bar", orient="h", sharex=False, margin_titles=True,
    errwidth=0.1, capsize=0.3, errorbar="sd", edgecolor="black",
    aspect=1, ax=ax[1])
plt.savefig(results_file.replace(".csv", f"_by-{col_val}_bar.pdf"),
            bbox_inches="tight")
plt.close()
# END OF SCRIPT