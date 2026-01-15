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
from utils.plots import R_violins, smooth_scatter
from utils.constants import SOURCE_DATA_DIR, FIGURE_DIR


####################################################################################
# FIG 1C - Violin plot
####################################################################################
# ATTENTION - THIS FILE NEEDS TO BE DOWNLOADED FIRST FROM GEO (GSE293866).
results_file = f"{SOURCE_DATA_DIR}/FIG_1/1CD/WGBS_DNMT1i_1kb_tiles_averages_shared.tsv"
plot_dir = FIGURE_DIR

R_violins(results_file, plot_dir=plot_dir)

####################################################################################
# FIG 1D - Smooth scatter plots
####################################################################################
df = pd.read_csv(results_file, sep="\t")
samples = [e for e in list(df.columns) if e not in ["chrom", "start", "end"]
           and not e.isdigit()]
df = df[samples]
df.dropna(inplace=True)
samples = ['gDNA_Day0', 'cfDNA_Day0', 'gDNA_Day4_DNMT1i', 'cfDNA_Day4_DNMT1i']
smooth_scatter(df, plot_dir=plot_dir, samples=samples, region_id="Tiles")

# END OF SCRIPT