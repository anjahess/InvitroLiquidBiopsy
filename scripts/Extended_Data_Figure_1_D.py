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
# Electropherogram plotting
####################################################################################
source = "./sourcedata/FIG_1EXT/E1D/"
df = pd.read_csv(f"{source}/FigE1D_sourcedata.csv", index_col=None)
print(df)
x="bp_pos"
y="normalized_fluorescent_units"
hue="timepoint"
style="line"
df["category"] = df["treatment"].map({"NONE":"untreated",
                                      "EDTA_50uM":"stabilized",
                                      "EDTA_25uM":"stabilized",
                                      "70°C":"stabilized","gDNA":"gDNA"})
palette=["#85ada3","#2d435b"]

####################################################################################
#  Show degradation under normal conditions (of cf and added gDNA)
####################################################################################
g = sns.FacetGrid(df, col="category", hue="timepoint", palette=palette)
g.map(sns.lineplot, x,y)
g.set(xscale='log')
plt.legend()
plt.savefig(f"{source}overview.pdf", bbox_inches='tight')
plt.close()


####################################################################################
#  Show stabilizing effect per treatment (70°C, EDTA)
####################################################################################
g = sns.FacetGrid(df, hue="treatment", row="timepoint")
g.map(sns.lineplot, x,y).add_legend()
g.set(xscale='log')
plt.savefig(f"{source}stabilizer.pdf", bbox_inches='tight')
plt.close()

####################################################################################
#  Show information per mESC cell line
####################################################################################
g = sns.FacetGrid(df, col="category", hue="timepoint", palette=palette,
                  row="line")
g.map(sns.lineplot, x,y)
g.set(xscale='log')
plt.legend()
plt.savefig(f"{source}treatment_per_cat.pdf", bbox_inches='tight')
plt.close()

# END OF SCRIPT
