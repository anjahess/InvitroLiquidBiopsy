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
from plots import *

results_file = "./sourcedata/FIG_2/2B/Bioreactor-hepatocytes_1kb_tiles_averages_shared.tsv"
plot_dir = "./sourcedata/FIG_2/2B/"
filter_file = "./sourcedata/FIG_2/2B/filter_stats.csv"

####################################################################################
# FIG 2B - Smooth scatter plots
####################################################################################
df = pd.read_csv(results_file, sep="\t")
samples = [e for e in list(df.columns) if e not in ["chrom", "start", "end"]
           and not e.isdigit()]
df = df[samples]
df.dropna(inplace=True)
print(df.shape)

samples = ['FoPHep_cfDNA_D10_merge.chm13', 'FoPHep_gDNA_D0_merge.chm13',
           'FoPHep_cfDNA_D20_merge.chm13', 'FoPHep_gDNA_D20_merge.chm13',
           'FoPHep_cfDNA_D0_merge.chm13', 'FoPHep_gDNA_D10_merge.chm13']
smooth_scatter(df, plot_dir=plot_dir, samples=samples, region_id="Tiles")

####################################################################################
# Fig 2B - CpG coverage plotting
####################################################################################
df = pd.read_csv(filter_file)
print(df.shape)
plt.figure(figsize=(2,3))
df["specimen"] = df["sample"].str.split('_', expand=True)[0]

x_value = "sample"
palette = ["salmon", "#2c4e60","#2c4e60","#2c4e60","#5f9ea0","#5f9ea0","#5f9ea0"]
sns.barplot(data=df, x=x_value, y="n_cpg", palette=palette, edgecolor="black")
sns.barplot(data=df, x=x_value, y="n_cpg_filtered", palette=palette, edgecolor="black")
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig(filter_file.replace(".csv", ".pdf"),
            bbox_inches="tight")
plt.close()

# END OF SCRIPT