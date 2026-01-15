"""

Scripts and source data to reproduce figures of the manuscript:
Hess, Anja et al. 2025: Non-disruptive monitoring of cellular dynamics with cell-free DNA methylation

Max Planck Institute for Molecular Genetics, Berlin, Germany

Author: Anja Hess
Date: 2026 JAN 15

"""
import os
import pandas as pd
from utils.constants import SOURCE_DATA_DIR, FIGURE_DIR
from utils.plots import smooth_scatter
import shutil

# Define directories and unzip the data
results_file =  f"{SOURCE_DATA_DIR}/FIG_2EXT/E2F/Bioreactor-hepatocytes_5hmC_1kb_tiles.zip"
unzip_dir = f"{SOURCE_DATA_DIR}/FIG_2EXT/E2F/"
unzipped_file = f"{results_file.replace('.zip','.tsv')}"
# Unzip the file
if not os.path.isfile(unzipped_file):
    shutil.unpack_archive(results_file,unzip_dir)
plot_dir = f"{FIGURE_DIR}/"

####################################################################################
# ED 2F - Smooth scatter plots
####################################################################################
df = pd.read_csv(unzipped_file, sep="\t")
samples = [e for e in list(df.columns) if e not in ["chrom", "start", "end"]
           and not e.isdigit()]
df = df[samples]
df.dropna(inplace=True)

samples = [['FoPHep_cfDNA_D0_merge_hydroxy.chm13', 'FoPHep_gDNA_D0_merge_hydroxy.chm13'],
           ['FoPHep_cfDNA_D10_merge_hydroxy.chm13', 'FoPHep_gDNA_D10_merge_hydroxy.chm13'],
           ['FoPHep_cfDNA_D20_merge_hydroxy.chm13', 'FoPHep_gDNA_D20_merge_hydroxy.chm13']]
for day in samples:
    smooth_scatter(df, plot_dir=plot_dir, samples=day, region_id="Tiles")
# END OF SCRIPT