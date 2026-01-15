"""

## Plotting functions

Scripts and source data to reproduce figures of the manuscript:
Hess, Anja et al. 2025: Non-disruptive monitoring of cellular dynamics with cell-free DNA methylation

Max Planck Institute for Molecular Genetics, Berlin, Germany

Author: Anja Hess
Date: 2025 APRIL 07

"""

import os.path
import sys
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from scipy.stats import pearsonr
pandas2ri.activate()
import rpy2.robjects as ro
sys.setrecursionlimit(100000)

from .constants import UTILS_DIR, FIGURE_DIR


def R_violins(filename, plot_dir="", samples=[]):
    """
    To show the correlation for every single CpG/region.
    :return:
    """


    #########################################################################
    # 0. Create results folder
    #########################################################################
    file_id = filename.rsplit("/",1)[1].replace(".tsv","")
    results_dir = f"{FIGURE_DIR}{file_id}_vio_methylation.pdf"
    if not os.path.isfile(filename):
        print("R Violin df: Not existing..", filename)
        exit()

    #########################################################################
    # 1. Sneak sample names
    #########################################################################
    sample_names = []
    sneak = pd.read_table(filename,  header=0, nrows=1)
    for e in sneak.columns:
        try:
            e = int(e)
        except ValueError:
            if "Unnamed" not in e and "chrom" not in e and "start" not in e and "end" not in e:
                sample_names.append(e)
    sample_names = sorted(sample_names)
    print(sample_names)

    if "3" in sneak.columns:
        cols_to_plot = sample_names + ["3"]
    else:
        cols_to_plot = sample_names

    df = pd.read_table(filename, header=0, usecols=cols_to_plot, index_col=None)
    original_number = len(df)
    nrows, ncols = df.shape

    ##########################################################################
    # Remove Positions where *all samples* have a NA
    ##########################################################################
    df = df.dropna(subset=sample_names, how="all")
    df_nona = df.dropna() # - see everything
    len_nona, n_cols = df_nona.shape
    title = (f"{filename.replace('.tsv', '').rsplit('/' ,1)[1]} \n"
             f"n = {nrows} (out of {original_number} possible incl chrX,Y)")
    outdir = filename.replace('.tsv', f'.pdf')
    print(nrows, ",", len_nona, "are shared across all samples.")
    if "3" in sneak.columns:
        try:
            df["3"] = df["3"].str.split(":" ,expand=True)[0]
        except:
            ""
        subvars = list(df["3"].unique())
        df = pd.melt(df, id_vars=["3"], value_vars=sample_names)
    else:
        subvars = False
        df = pd.melt(df) # will always give you variable, value

    #################################################################
    # 2. Call R-version
    #################################################################
    # Loading the function defined in R.
    r = robjects.r
    r['source'](f"{UTILS_DIR}/R_plots.R")
    plotvio= robjects.globalenv['plotvio']

    #################################################################
    # 2.1 Plot everything
    #################################################################
    df_r = ro.conversion.py2rpy(df)
    plotvio(df_r, results_dir, title, sample_names, subvars)

    #################################################################
    # 2.1 Plot only shared
    #################################################################
    if "3" in sneak.columns:
        try:
            df_nona["3"] = df_nona["3"].str.split(":" ,expand=True)[0]
        except:
            ""
        subvars = list(df_nona["3"].unique())
        df_nona = pd.melt(df_nona, id_vars=["3"], value_vars=sample_names)
    else:
        subvars = False
        df_nona = pd.melt(df_nona) # will always give you variable, value
    df_r_nona = ro.conversion.py2rpy(df_nona)
    title = (f"{filename.replace('.tsv', '').rsplit('/' ,1)[1]} \n"
             f"n = {len_nona} (out of {original_number} possible incl chrX,Y)")
    results_dir_nona = results_dir.replace(".pdf" ,"_nona.pdf")

    try:
        plotvio(df_r_nona, results_dir_nona, title, sample_names, subvars)
    except:
        print("Unable to plot only shared (no left)")
    # END OF FUNCTION



def smooth_scatter(df, plot_dir="", samples=[], subsampling=100000, region_id=""):
    """
    To show the concordance of DNA(hydroxy-)methylation for  individual CpGs or regions.
    :param df: str
    :param plot_dir: str
    :param samples: list
    :param subsampling: int, only for visualization, not concordance calculation
    :param region_id: str
    :return: Creates smooth scatter correlation plots with metrics to the plot_dir.
    """

    sns.set_theme(style='ticks')

    #########################################################################
    # Random sampling (speed)
    #########################################################################
    try:
        df_cut = df.sample(n=subsampling)
        # this is ONLY for visualization not concordance calculation!
    except:
        print(f"Not enough regions for subsampling to {subsampling}")
        print(f"Please choose a lower number for smooth_scatter parameter 'subsampling'.")
        exit()

    sample_pairs = []

    for sample in samples:
        for other_sample in [e for e in samples if e not in sample]:

            if f"{sample}-{other_sample}" in sample_pairs:
                print(f"{sample}-{other_sample} exists")
                continue
            sample_pairs.append(f"{sample}{other_sample}")
            sample_pairs.append(f"{other_sample}-{sample}")

            #################################################################
            # Concordance
            # Calculate the percentage located on the +-0.1 diagonal
            # ! ! ! NOTE: result will be valid for then ENTIRE dataset
            # while the plot (for memory reasons) will be on the subset
            #################################################################
            df["delta"] = abs(df[sample] - df[other_sample])
            df["delta=<0.1"] = np.where(df['delta'] <= 0.1, True, False)

            result_df = df["delta=<0.1"].value_counts()
            results_df_rel = round(df["delta=<0.1"].value_counts(normalize=True),2)
            result = pd.concat([result_df, results_df_rel],
                               axis=1).reindex(result_df.index)
            print(f"--- Results for {sample} vs {other_sample} ----")
            print(result)

            # Calculate correlation
            corr = pearsonr(df[sample], df[other_sample])
            corr = [np.round(c, 3) for c in corr]
            print(corr)
            print("Total = ", len(df))
            print("-----------")

            ##################################################################
            # Pearson correlation
            ##################################################################
            title = \
                (f"{region_id}\n{sample} vs. {other_sample}, "
                 f"n = {len(df[sample])} regions ({subsampling} plotted)\n"
                f"{result}\n"
                 f"person corr {corr}")

            #################################################################
            # Call R-version --> NOTE SUBSAMPLING FOR PLOT DOES NOT AFFECT
            # CONCORDANCE RESULT
            #################################################################
            # Loading the function defined in R.
            r = robjects.r
            r['source'](f"{UTILS_DIR}/R_plots.R")
            plot_dir = plot_dir + f"{region_id}_{sample}_vs_{other_sample}_R.pdf"
            plotsmoothscatter = robjects.globalenv['plotsmoothscatter']
            df_r = ro.conversion.py2rpy(df_cut)
            plotsmoothscatter(df_r, plot_dir, title, sample, other_sample)
            # END OF SAMPLE COMP LOOP
    # END OF FUNCTION


def clustermap(file, method="complete", metric="euclidean", outdir=""):
    """
    Clustermap seaborn
    :param file: str
    :param method: str
    :param metric: str
    :param outdir: str
    :return:
    """

    cmap = sns.diverging_palette(220, 20, as_cmap=True)
    df = pd.read_table(file, sep="\t", index_col=0)
    df.dropna(inplace=True)
    cols_of_interest = [e for e in df.columns if "DNA" in e]
    df = df[cols_of_interest]
    title = f"n = {len(df)}"
    sns.clustermap(df, rasterized=True, cmap=cmap, method=method, metric=metric, vmin=0, vmax=1)
    plt.title(title)
    file_id = file.rsplit("/",1)[1].replace(".tsv","")
    plt.savefig(f"{outdir}{file_id}_cluster.pdf", bbox_inches="tight")
    plt.close()
    # END OF FUNCTION

# END OF SCRIPT