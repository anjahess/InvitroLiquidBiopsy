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
from utils.constants import SOURCE_DATA_DIR, FIGURE_DIR

####################################################################################
# Read coverage plotting
####################################################################################
source_dir = f"{SOURCE_DATA_DIR}/FIG_1EXT/E1G/WGBS-read-cov-metrics.txt"
palettediff2 = [ "#7db1b3", "#355269"]
df = pd.read_table(source_dir, header=0)
plt.figure(figsize=(2,3))
sns.violinplot(data=df, x="specimen", y="mean", hue="specimen", palette=palettediff2,
               order=["cell-free DNA", "gDNA"])
sns.swarmplot(data=df, x="specimen", y="mean", hue="specimen", palette=palettediff2,
              order=["cell-free DNA", "gDNA"], linewidth=0.5,edgecolor="black", s=4)
plt.ylim(0, 50)
plt.ylabel("Average coverage (x)")
plt.savefig(f"{FIGURE_DIR}mESC_WGBS_reads.pdf",
            bbox_inches="tight")
# END OF SCRIPT