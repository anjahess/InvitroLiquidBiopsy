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
import numpy as np
import matplotlib
from utils.constants import SOURCE_DATA_DIR, FIGURE_DIR
####################################################################################
# EXT FIG 1C - Line plots visualizing cf yield dynamics with apoptosis rates
####################################################################################

def plot_csv(path_to_file, delimiter=',', palette=["grey"]):
    """
    :param path_to_file:
    :return:
    """

    df_orig = pd.read_csv(path_to_file, delimiter=delimiter, header=0)
    df_orig.replace("#DIV/0!", np.nan, inplace=True)
    df_orig.replace("na", np.nan, inplace=True)
    #################################################################################
    # Iterate through experiments
    #################################################################################
    for experiment in df_orig["Experiment"].unique():
        print(f"--- Analyzing experiment {experiment} ---")
        df = df_orig[df_orig["Experiment"] == experiment]
        df["ng_per_10K_cells"] = df["ng_per_10K_cells"].astype(float)
        df["pg_per_cell"] = df["pg_per_cell"].astype(float)
        df["total_ng"] = df["total_ng"].astype(float)
        try:
            df["time_point"] =  df["time_point"].astype(int)
            df["tp_treat"] = (df["time_point"].astype(str) + " " +
                              df["treatment"].astype(str))
            df.sort_values(by="tp_treat", inplace=True)
        except:
            ""
        ##############################################################################
        # Generate the line plots
        ##############################################################################
        g = sns.PairGrid(df, y_vars=["FACS_apo_percent",  "total_ng",
                                     "ng_per_10K_cells","pg_per_cell"],
                         x_vars=["time_point"], aspect=3)
        g.map(sns.lineplot, data=df, hue="treatment", style="cell_line",
              style_order=["V6.5","F1G4","KH2"], markersize=5, linewidth=0.5,
              markeredgecolor="black", markeredgewidth=0.1,
              palette=palette, err_style="band", errorbar=('ci',95),
              markers=True, dashes=False).add_legend()
        g.fig.set_size_inches(2, 4)

        ##############################################################################
        # Customize limitis
        ##############################################################################
        lims = [(0, 40), (0, 250), (0,5), (0,0.5)]
        loc = matplotlib.ticker.MultipleLocator(0.5)
        tick_inc = [2, 2, 2, 2]
        tick_inc_y = [10, 100, 2, 0.1]

        for ax, ylims, xticks, yticks in zip(g.axes.flat, lims, tick_inc, tick_inc_y):
            ax.set_ylim(ylims)
            ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(xticks))
            ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(yticks))

        plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
        plt.savefig(f"{FIGURE_DIR}/mESC_apo_cfDNAyield.pdf",
                    bbox_inches="tight")
        plt.close()

        ##############################################################################
        # Calculate the slope (cfDNA yield / 10K cells / hour)
        ##############################################################################
        x = df["time_point"]
        y = df["ng_per_10K_cells"]
        slope_intercept = np.polyfit(x, y, 1)
        print(slope_intercept) # [0.10809674 0.04921136]
        # END OF FUNCTION

root_dir = f"{SOURCE_DATA_DIR}FIG_1EXT/E1C/FigE1C_sourcedata.csv"
plot_csv(root_dir, delimiter="\t", palette=["#477b80","#b85c54"])

# END OF SCRIPT