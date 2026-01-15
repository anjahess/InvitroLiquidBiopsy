"""

Scripts and source data to reproduce figures of the manuscript:
Hess, Anja et al. 2025: Non-disruptive monitoring of cellular dynamics with cell-free DNA methylation

Max Planck Institute for Molecular Genetics, Berlin, Germany

Author: Anja Hess
Date: 2025 APRIL 07

"""
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from utils.plots import R_violins
from utils.constants import SOURCE_DATA_DIR, FIGURE_DIR, TABLE_DIR

results_file = f"{SOURCE_DATA_DIR}/FIG_2/2DE/Bioreactor-hepatocytes_hep-enhancers-ES-filtered.tsv"
anno_file = f"{SOURCE_DATA_DIR}/FIG_2/2DE/Bioreactor-hepatocytes_hep-enhancers-ES-filtered_anno.tsv"
plot_dir = FIGURE_DIR

####################################################################################
# FIG 2C - Violin plot enhancer methylation
####################################################################################
# For demonstration, lower the subsampling to e.g. 1K, will be faster
R_violins(results_file, plot_dir=plot_dir)

####################################################################################
# FIG 2D - Heatmap enhancer methylation
####################################################################################
cmap=sns.diverging_palette(220, 20, as_cmap=True)
df = pd.read_table(anno_file, sep="\t", index_col=0)
df["index_no"] = df.index
df.set_index("gene", inplace=True)
df.dropna(inplace=True)
title="Hepatocyte Enhancers"
cols_of_interest = [e for e in df.columns if "DNA" in e]
df = df[cols_of_interest]
y_labels_ticks = df.index.tolist()
y_labels_new = []
for gene in y_labels_ticks:
    if type(gene) == float:
        y_labels_new.append(" ")
    else:
        y_labels_new.append(gene)
print(len(df), " total regions in clustermap.")
cm = sns.clustermap(df, rasterized=True, cmap=cmap, vmin=0, vmax=1,
                    yticklabels=y_labels_new)
plt.savefig(f"{FIGURE_DIR}Bioreactor-hepatocytes_cluster_n{len(df)}.pdf",
            bbox_inches="tight")
plt.close()

