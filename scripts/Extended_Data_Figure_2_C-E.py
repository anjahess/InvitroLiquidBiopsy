"""

An easy copy-paste-to-folder script that will systematically plot pycoQC reports to figures

Author: Anja Hess
Date: 2025 OCT 13


"""
import os.path
import matplotlib.pyplot as plt
import json
import glob
import seaborn as sns
import pandas as pd
import numpy as np
from utils.constants import SOURCE_DATA_DIR, FIGURE_DIR, TABLE_DIR


# Set variables and load the dataframe
source_subdir = f"{SOURCE_DATA_DIR}FIG_2EXT/E2C-E/"
source_file = f"{source_subdir}cpg_df.csv"
palette = ["#5f9ea0", "#2c4e60", "salmon", "blue"]
order = ["gDNA", "cell-free DNA", "reference"]
df = pd.read_csv(source_file)
df = df.replace(pd.NA, np.nan)
df["input_ng_normalized_bef-adapter-lig"] = (
    df["input_ng_normalized_bef-adapter-lig"].astype(float))
col_val = "biol_rep"


####################################################################################
# ED 2E - CpG FILTER STATS PER BIOLOGICAL REPLICATE
####################################################################################
fig, ax = plt.subplots(2,2)

g = sns.relplot(
    data=df, x="filtering", y="n_cpg_filtered",
    hue="specimen", hue_order=order,
    style="specimen", style_order=order,
    errorbar=('ci', 100),
    col=col_val, kind="line", palette=palette, legend=False,
    linewidth=0.5, height=5, aspect=.75, facet_kws=dict(sharex=False),
)
g.map(sns.scatterplot,"filtering","n_cpg_filtered", "specimen",
      hue_order=order,
      palette=palette,
      legend=True
      ).add_legend()
plt.ylim(0, None)

plt.savefig(f"{FIGURE_DIR}Nanopore_CpG_by-{col_val}_line.pdf",bbox_inches="tight")
plt.close()


####################################################################################
# ED 2C - BAR PLOTS FOR NANOGRAMS INPUT
####################################################################################
df = df[df["specimen"] != "reference"]
x="input_ng_normalized_bef-adapter-lig"
y="specimen"
# Sub-select df @ 1X (input is constant)
df = df[df["filtering"] == "1X"]
g = sns.catplot(
    data=df, palette=[ "#2d435b","#85ada3",],
    x=x, y=y, orient="h", col=col_val, kind="bar", sharex=False,
    margin_titles=True,errwidth=0.1, capsize=0.3, errorbar="ci",
    edgecolor="black",aspect=1)
g.map(sns.scatterplot, x,y,y, palette=["#85ada3", "#2d435b"],
      hue_order=order, s=70, edgecolor="black", alpha=1).add_legend()
plt.savefig(f"{FIGURE_DIR}Nanopore_inputs_by-{col_val}_bar.pdf",bbox_inches="tight")
plt.close()


####################################################################################
# E2 C(lower) & D - Nanopore metrics
####################################################################################
frag_df_path = f"{TABLE_DIR}Nanopore_fragment_df.csv"
cov_df_path = f"{TABLE_DIR}Nanopore_coverage_df.csv"
palette = ["#355269", "#7db1b3"]
bp_term = "bp_size"
xlim = 14000
palette = ["#2d435b","#85ada3"]

# Set empty data points
header = []
mean_cov = []
mean_cov_all = []

####################################################################################
# STEP 1: collect data from json files (generated with pycoQC)
####################################################################################
if not os.path.isfile(frag_df_path) or not os.path.isfile(cov_df_path):
    fragment_df = pd.DataFrame(index=np.arange(50000))
    fragment_df["bp_size"] = fragment_df.index
    for file in glob.glob(source_subdir + "*.json"):
        file_id = file.rsplit("/", 1)[1].rsplit(".json", 1)[0]
        print("--- Gathering ONT metrics for", file_id)
        header.append(file_id)
        ############################################################################
        # Load json
        ############################################################################
        with open(file) as f:
            data =  json.load(f)
        pass_reads = data["Pass Reads"]
        all_reads = data["All Reads"]
        alignment = pass_reads["alignment"]
        alignment_all = all_reads["alignment"]

        ############################################################################
        # Retrieve average coverage and read histograms (aligned data)
        ############################################################################
        for chapter in alignment:
            if chapter == "mean_coverage":
                mean_cov.append(alignment[chapter])
            if "len_hist" in chapter:
                # fragment df
                x = [int(e) for e in alignment[chapter]["x"]]
                y = alignment[chapter]["y"]
                frag_df = pd.DataFrame(y,x)
                frag_df["bp_size"] = frag_df.index
                frag_df.rename(columns={0: file_id}, inplace=True)
                ####################################################################
                # Integrate into the shared fragment dataframe
                ####################################################################
                fragment_df = pd.merge(fragment_df, frag_df, on="bp_size",
                                       how="outer")
        for chapter in alignment_all:
            if chapter == "mean_coverage":
                mean_cov_all.append(alignment_all[chapter])
    fragment_df.to_csv(frag_df_path)


    ################################################################################
    # Save genomic coverage data
    ################################################################################
    cov_df = pd.DataFrame(np.column_stack([mean_cov,mean_cov_all]),
                          index = header, columns = ["mean_coverage_pass_reads",
                                                     "mean_coverage_all"])
    cov_df["sample"] = cov_df.index
    cov_df["specimen"] = cov_df["sample"].str.rsplit(
        'DNA', expand=True)[0].str.split("_", expand=True)[2]
    cov_df["biol_rep"] = cov_df["sample"].str.rsplit('_', expand=True)[1]
    cov_df["biol_rep"] = cov_df["biol_rep"].str.rsplit("N",
                                                       expand=True)[1].astype(int)
    cov_df["day"] = cov_df["sample"
    ].str.rsplit('_', expand=True)[3].str.rsplit('D',
                                                 expand=True)[1].astype(int)
    cov_df.to_csv(cov_df_path, sep="\t", index=False)


####################################################################################
# ED C (lower) - Nanopore read coverage
####################################################################################
cov_df = pd.read_table(cov_df_path)

# PASS READS
g = sns.FacetGrid(data=cov_df, col="biol_rep", sharex=True, sharey=True)
g.map(sns.boxplot, "specimen", "mean_coverage_pass_reads", palette=palette)
g.map(sns.swarmplot,"specimen", "mean_coverage_pass_reads",
      "day",s=9, linewidth=0.5, edgecolor="black", palette="Set2")
legend = plt.legend(bbox_to_anchor=(1.05, 1.0), loc="upper left")
plt.ylim(0, None)
plt.ylabel("Mean Genomic Read Coverage (X)")
plt.savefig(f"{FIGURE_DIR}Nanopore_mean_read_coverage.pdf",
            bbox_inches='tight')
plt.close()


####################################################################################
# ED C (lower) - Nanopore fragment lengths
####################################################################################
df = pd.read_csv(frag_df_path, index_col=0)
cols_of_interest = [e for e in df.columns if ("DNA" in e) or ("bp_size" in e)]
df = df[cols_of_interest]

# Interpolate from PycoQC's ONT histogram
df = df.interpolate(axis=0)
df = df[df[bp_term] < xlim]

# Column normalize
cols_of_interest = [e for e in df.columns if ("DNA" in e)]
df = (df - df.min()) / (df.max() - df.min())
df[bp_term] = df.index

# Long dataframe
cols_of_interest = [e for e in df.columns if ("DNA" in e)]
df_long = df.melt(id_vars=bp_term, value_vars=cols_of_interest)
df_long[bp_term] = df_long[bp_term].astype(int)
df_long["condition"] = df_long["variable"
].str.split('_', expand=True)[3].str.split("D",expand=True)[1].astype(int)
df_long["specimen"] = df_long["variable"].str.split('_', expand=True)[2]
df_long["replicate"] = df_long["variable"].str.split('_', expand=True)[1]
df_long["biological_replicate"] = df_long["replicate"].map({ "FoPHepN2":2,
                                                             "FoPHepN3":3,
                                                             "FoPHepN4":4,
                                                             "FoPHepN1":1})
df_long["cond-day"] = df_long["condition"].astype(str) + df_long["specimen"].astype(str)
df_long.dropna(inplace=True)
df_long.sort_values(by="cond-day",inplace=True)
palette = ["#355269", "#7db1b3", "#4f728f", "#9dcdcf", "#6888a1", "#b6d2d4",]
#start, end, step (resolve small frags w/p having too many bins)
order = list(range(0, 300, 10))
order_2 = list(range(300, xlim, 100))
order = order + order_2
g = sns.FacetGrid(df_long, col="biological_replicate", hue="cond-day",
                  sharey=False, sharex=False, aspect=1.5,
                  palette=palette)
g.map(sns.barplot, "bp_size", "value", errorbar=None,
      order=order).add_legend()
g.set_xticklabels(step=20, fontsize=6, rotation=90)
plt.xlim(0, None)
sns.despine()
plt.savefig(f"{FIGURE_DIR}Nanopore_read_sizes.pdf", dpi=50)

# END OF SCRIPT