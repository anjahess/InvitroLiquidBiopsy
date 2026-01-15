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
from utils.constants import SOUCE_DATA_DIR, FIGURE_DIR

####################################################################################
# qPCR expression plotting
####################################################################################
source_dir = f"{SOUCE_DATA_DIR}/FIG_2EXT/E2A/Bioreactor_hepatocytes_qPCR.csv"

palettediff2 = ["#7db1b3", "#355269"]
order = [""]
df = pd.read_table(source_dir, header=0)
plt.figure(figsize=(2,3))
sns.barplot(data=df, x="Gene", y="2^(-Delta Ct)", palette=palettediff2,
            errwidth=0.1, capsize=0.3, errorbar="sd", edgecolor="black")
sns.swarmplot(data=df, x="Gene", y="2^(-Delta Ct)", hue="Experiment",
              palette=palettediff2, linewidth=0.5, edgecolor="black", s=4)
plt.ylabel("mRNA expression 2^(-Delta Ct)")
plt.tight_layout()
plt.savefig(f"{FIGURE_DIR}Bioreactor_hepatocytes_qPCR.pdf",
            bbox_inches="tight")
# END OF SCRIPT