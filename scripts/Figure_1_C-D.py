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

####################################################################################
# FIG 1C - Violin plot
####################################################################################
results_file = "./sourcedata/FIG_1/1CD/WGBS_DNMT1i_1kb_tiles_averages_shared.tsv"
plot_dir = "./sourcedata/FIG_1/1CD/"

# For demonstration, do random subsampling to e.g. 1K, to be faster
R_violins(results_file, plot_dir=plot_dir, subsampling=1000)

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