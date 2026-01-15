"""

Scripts and source data to reproduce figures of the manuscript:
Hess, Anja et al. 2025: Non-disruptive monitoring of cellular dynamics with cell-free DNA methylation

Max Planck Institute for Molecular Genetics, Berlin, Germany

Author: Anja Hess
Date: 2025 OCTOBER 20

"""
import os
import glob
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
import scikit_posthocs as sp
from scipy.stats import f_oneway
from utils.constants import SOURCE_DATA_DIR, FIGURE_DIR, TABLE_DIR
from utils.plots import clustermap, R_violins

file_dir = f"{SOURCE_DATA_DIR}FIG_3EXT/"
palette= ["#355269", "#7db1b3"]

def run_stats(df, variable="", category=""):
    """
    Function to perform scipy.stats' the per group test to infer statistical
    significance
    :param df: pandas.DataFrame
    :param variable: continuous variable
    :param category: categorical variable
    :return: statistics per group in a dataframe

    """
    test_performed = "1W-ANOVA" # with fdr_bh correction"
    sdata = []
    #####################################################################
    # 1. Collect numerical values
    #####################################################################
    groups = []
    sdict = {}
    names = []
    p_value = signi = results = None

    for group in df[category].unique():
        group_data = df[df[category] == group][variable]
        group_data = list(group_data)
        if not group_data:
                print(f"No data found for group {group}.")
                continue
        groups.append(group_data)
        sdict.update({str(group): group_data})
        names.append(str(group))

    ##################################################################
    # 2. Run test
    ##################################################################
    stats, p_value = f_oneway(*groups)
    if p_value > 0.05:
        signi = False
    # 2. If the test says groups are different do a posthoc
    if p_value < 0.05:
        signi = True
        groups_for_posthoc = np.asarray(groups, dtype="object")
        results = sp.posthoc_conover(groups_for_posthoc)
        results.columns = names
        results["condition"] = names
        results.set_index("condition", inplace=True)
    sdata.append([test_performed, p_value, signi, results,
                         sdict])
    #####################################################################
    # 2. Generate df from storage
    #####################################################################
    sdf = pd.DataFrame(sdata, columns=["test_performed", "p_value",
                                       "p<0.05", "posthoc_p_values", "groups"])
    return sdf
    # END OF STATS FUNCTION


def delta_plot_stats(file, outdir_plots="", outdir_tables=""):
    """
    Function for the delta line plots (difference of mean methylation)
    and statistical comparison of replicates
    :param df: path to file
    :param outdir_plots: str
    :param outdir_tables: str
    :return: creates line + violin plots of average methylation and statistics
    tables
    """
    ################################################################################
    # Initiate results folder
    ################################################################################
    file_id = file.rsplit("/",1)[1].replace(".tsv","")
    outdir_file = f"{outdir_tables}{file_id}_mean.tsv"

    ################################################################################
    # Create data table
    ################################################################################
    if not os.path.isfile(outdir_file):
        print("--- Generating stats table.")
        df = pd.read_table(file, sep='\t', index_col=0)
        try:
            df = df.drop(columns=['0','1','2'])
        except:
            try:
                df = df.drop(columns=['start', 'end', '3',])
            except:
                try:
                    df = df.drop(columns=['start', 'end'])
                except:
                    df = pd.read_table(file, sep='\t',
                                       )
                    df = df.drop(columns=['0', '1', '2'])

        #############################################################################
        # Annotate metadata
        #############################################################################
        mean_df = df.mean().to_frame()
        mean_df["sample"] = mean_df.index.astype(str)
        mean_df["sample"] = mean_df["sample"
        ].str.split(".", expand=True)[0].str.split("_merge", expand=True)[0]
        mean_df["condition"] = mean_df["sample"].str.split("_", expand=True)[1]
        mean_df["day"] = mean_df["sample"].str.rsplit("DNA_", expand=True)[1]
        mean_df["day"] = mean_df["day"].str.rsplit("D", expand=True)[1]
        mean_df["day"] = mean_df["day"].str.rsplit("replace", expand=True)[0].astype(int)
        mean_df["biol_rep"] = mean_df["sample"
        ].str.rsplit("_",expand=True)[0].map({"FoPHepN2":2, "FoPHepN3":3,
                                              "FoPHepN4":4, "FoPHepN1":1})
        mean_df.rename(columns={0: "mean_meth"}, inplace=True)
        mean_df.to_csv(outdir_file, sep='\t')

    ################################################################################
    # Plotting
    ################################################################################
    if os.path.isfile(outdir_file):
        df = pd.read_table(outdir_file, sep='\t', index_col=0)

        #############################################################################
        # ED 3B - Delta-plot per replicate
        #############################################################################
        g = sns.FacetGrid(data=df, col="biol_rep",palette=palette,
                          hue="condition",sharex=False, sharey=True, aspect=0.8)
        g.map(sns.lineplot, "day", "mean_meth")
        plt.ylim(0,1)
        delta_dict = {}
        # Annotate the avergae delta methylation on the plot
        for i, biol_rep in enumerate(df["biol_rep"].unique()):
             sub_df = df[df["biol_rep"] == biol_rep]
             df2 = sub_df.pivot_table(index="day", columns="condition", values="mean_meth",
                                aggfunc='mean')
             df2["delta"] = (df2["cfDNA"]-df2["gDNA"]).abs()
             mean_delta = round(df2["delta"].mean(), 3)
             delta_dict[str(biol_rep)] = str(mean_delta)
        plt.suptitle(f"Delta meth: {str(delta_dict)}")
        plt.tight_layout()
        plt.legend()
        plt.savefig(f"{outdir_plots}{file_id}_delta.pdf")
        plt.close()

        #############################################################################
        # ED 3D - Statistical comparison (average of all three replicates)
        #############################################################################
        cond_dict={}
        for cond in df["condition"].unique():
            print(f"--- {cond}")
            stats = run_stats(df[df["condition"]==cond], category="day",
                              variable="mean_meth")
            stats.to_csv(outdir_file.replace('.tsv',
                                             f'_stats_{cond}.csv'), sep='\t')
            p_val = stats["p_value"]
            cond_dict[cond] = p_val
        sns.violinplot(data=df, x="day", y="mean_meth", hue="condition", dodge=True,
                       color="white", palette=["white"], density_norm="count",)
        sns.swarmplot(data=df, x="day", y="mean_meth", hue="condition", dodge=True,
                      palette=palette, s=10, edgecolor="black", linewidth=0.2)
        plt.title(f"{cond_dict}")
        plt.xticks(rotation=90)
        plt.ylim(0, 1)
        plt.ylabel("Average CpG Methylation")
        plt.savefig(f"{outdir_plots}{file_id}_stats.pdf")
        plt.close()
        # END OF LOOP
    # END OF FUNCTION


#######################################################################################
# START THE RUN
#######################################################################################
for file in glob.glob(file_dir + "*.tsv"):
    if "mean" in file or "difference" in file or "merge" in file:
        continue
    file_id = file.rsplit('/',1)[1].replace(".tsv","")
    ####################################################################################
    # ED 3B & ED 3D Delta line plot / Violin plots with statistical comp.
    ####################################################################################
    print(f"---- Delta plots and statistics {file}")
    delta_plot_stats(file, outdir_plots=FIGURE_DIR, outdir_tables=TABLE_DIR)

    ####################################################################################
    # ED 3A & ED 3C : Violin plots - and clustermaps for DNA methylation
    ####################################################################################
    R_violins(file)
    clustermap(file, outdir=FIGURE_DIR)

    ####################################################################################
    # Average of all replicates
    ####################################################################################
    average_file = f"{TABLE_DIR}{file_id}_merge.tsv"

    if not os.path.isfile(average_file):
        df = pd.read_table(file, sep="\t", index_col=0)
        d_20g = [e for e in df.columns if "D20" in e and "gDNA" in e]
        d_10g = [e for e in df.columns if "D10" in e and "gDNA" in e]
        d_0g = [e for e in df.columns if "D0" in e and "gDNA" in e]
        d_20cf = [e for e in df.columns if "D20" in e and "cfDNA" in e]
        d_10cf = [e for e in df.columns if "D10" in e and "cfDNA" in e]
        d_0cf = [e for e in df.columns if "D0" in e and "cfDNA" in e]
        df["D20_cfDNA_avg"] = pd.Series(np.nanmean(df[d_20cf], axis=1))
        df['D10_cfDNA_avg'] = pd.Series(np.nanmean(df[d_10cf], axis=1))
        df['D0_cfDNA_avg'] = pd.Series(np.nanmean(df[d_0cf], axis=1))
        df['D20_gDNA_avg'] = pd.Series(np.nanmean(df[d_20g], axis=1))
        df['D10_gDNA_avg'] = pd.Series(np.nanmean(df[d_10g], axis=1))
        df['D0_gDNA_avg'] = pd.Series(np.nanmean(df[d_0g], axis=1))
        cols_of_interest = [e for e in df.columns if "avg" in e or e=="0" or e=="1" or e=="2"]
        df = df[cols_of_interest]
        df.to_csv(average_file, sep="\t")
    R_violins(average_file)
    clustermap(average_file, outdir=FIGURE_DIR)
    # END OF FUNCTION
# END OF SCRIPT
