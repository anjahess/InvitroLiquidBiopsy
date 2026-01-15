"""

Scripts and source data to reproduce figures of the manuscript:
Hess, Anja et al. 2025: Non-disruptive monitoring of cellular dynamics with cell-free DNA methylation

Max Planck Institute for Molecular Genetics, Berlin, Germany

Author: Anja Hess
Date: 2025 APRIL 07

"""
import scanpy as sc

####################################################################################
# FIG 1F - Principal Component Analysis
####################################################################################

source=f"./sourcedata/FIG_1/1F/"
sc.settings.figdir = source
adata = sc.read_h5ad(f"{source}WGBS_PCA_DNMT1i_mESC.h5ad")

sc.pl.pca(adata, save=f"mean-meth.pdf", color="mean_meth",
          color_map='GnBu', size=400, vmin=0, vmax=1,
          annotate_var_explained=True)

# END OF SCRIPT

