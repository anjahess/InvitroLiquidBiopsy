"""

Scripts and source data to reproduce figures of the manuscript:
Hess, Anja et al. 2025: Non-disruptive monitoring of cellular dynamics with cell-free DNA methylation

Max Planck Institute for Molecular Genetics, Berlin, Germany

Author: Anja Hess
Date: 2025 APRIL 07

"""
from plots import *
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import glob
results_file = "./sourcedata/FIG_2EXT/E2E-F/Bioreactor-hepatocytes_hep-markers-ES-filtered.tsv"
plot_dir = "./sourcedata/FIG_2EXT/E2E-F/"

####################################################################################
# FIG E2E -cell type marker methylation (based on Loyfer et al, 2023)
####################################################################################
# For demonstration, lower the subsampling to e.g. 1K, will be faster
R_violins(results_file, plot_dir=plot_dir, subsampling=100000)

####################################################################################
# FIG E2E - Delta plot
####################################################################################
for file in glob.glob(plot_dir + "*.tsv"):
    if "mean" in file or "difference" in file:
        continue

    if not os.path.isfile(file.replace('.tsv', '_mean.tsv')):
        df = pd.read_table(file, sep='\t', index_col=0)

        try:
            df = df.drop(columns=['0','1','2'])
        except:
            try:
                df = df.drop(columns=['start', 'end', '3'])
            except:
                try:
                    df = df.drop(columns=['start', 'end'])
                except:
                    df = pd.read_table(file, sep='\t',
                                       )
                    df = df.drop(columns=['0', '1', '2'])
        mean_df = df.mean().to_frame()
        mean_df["sample"] = mean_df.index.astype(str)
        mean_df["sample"] = mean_df["sample"
        ].str.split(".", expand=True)[0].str.split("_merge", expand=True)[0]
        mean_df["condition"] = mean_df["sample"].str.split("_", expand=True)[1]
        mean_df["day"] = mean_df["sample"].str.rsplit("DNA_", expand=True)[1]
        mean_df["day"] = mean_df["day"].str.rsplit("D", expand=True)[1]
        mean_df["day"] = mean_df["day"].str.rsplit("replace",
                                                   expand=True)[0].astype(int)
        print(mean_df["day"])
        mean_df.rename(columns={0: "mean_meth"}, inplace=True)
        mean_df.to_csv(file.replace('.tsv', '_mean.tsv'), sep='\t')

    if os.path.isfile(file.replace('.tsv', '_mean.tsv')):
        df = pd.read_table(file.replace('.tsv', '_mean.tsv'),
                           sep='\t', index_col=0)
        df2 = df.pivot(index="day", columns="condition", values="mean_meth")
        df2["delta"] = (df2["cfDNA"]-df2["gDNA"]).abs()
        mean_delta = df2["delta"].mean()
        print(mean_delta)
        df2.to_csv(file.replace(".tsv", f"_difference.tsv"))
        sns.lineplot(data=df, x='day', y='mean_meth', hue="condition",
                         palette=["#2d435b","#85ada3", ], )
        plt.title(f"Av Methylation Delta g/cf = {mean_delta}")
        plt.ylim(0, 1)
        plt.ylabel("Average CpG Methylation")
        plt.savefig(file.replace(".tsv", f"_delta.pdf"))
        plt.show()
        plt.close()

####################################################################################
# FIG E2E - Heatmap cell type marker methylation (based on Loyfer et al, 2023)
####################################################################################
cmap=sns.diverging_palette(220, 20, as_cmap=True)
df = pd.read_table(results_file, sep="\t", index_col=0)
df.dropna(inplace=True)
title="Hepatocyte Markers (based on Loyfer et al, 2023)"
cols_of_interest = [e for e in df.columns if "DNA" in e]
df = df[cols_of_interest]
y_labels_ticks = df.index.tolist()
y_labels_new = []
for gene in y_labels_ticks:
    if type(gene) == float:
        y_labels_new.append(" ")
    else:
        y_labels_new.append(gene)
cm = sns.clustermap(df, rasterized=True, cmap=cmap, vmin=0, vmax=1,
                    yticklabels=y_labels_new)
plt.savefig(results_file.replace(".tsv", f"_cluster.pdf"),
            bbox_inches="tight")
plt.close()

# END OF SCRIPT