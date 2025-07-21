# R Script to generate genome distribution visualization in
# karyogram-style from bam files.
# !!! YOU NEED TO RUN the Fig1_E_kary.sh FIRST !!!


# AUTHOR: Anja Hess
# DATE: 20241120
# USAGE: Rscript karyo_per_chromo.R BAM1 BAM2 chr1 mm10 plot_name
# Package: Gel B, Serra E (2017). “karyoploteR : an R / Bioconductor
# package to plot customizable genomes displaying arbitrary data.”
# Bioinformatics, 33(19), 3088-3090. doi:10.1093/bioinformatics/btx346.
# Tutorial: https://bernatgel.github.io/karyoploter_tutorial/


#############################################################################
# 0. Load dependencies
#############################################################################
library(karyoploteR)
library("BSgenome")
library(BiocFileCache)
library("regioneR")
bfc <- BiocFileCache(ask=FALSE)

#############################################################################
# 1. Parse the arguments
#############################################################################
#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

bam_cfDNA <- args[1]
bam_gDNA  <- args[2]
chr <- args[3]
genome <- args[4]
pdf_name <- args[5]
colors <- c("cell-free DNA"="#2c4e60", "gDNA"="#5f9ea0")
print(pdf_name)

#############################################################################
# 2.  Plot karyogram for BAMdensity (genome-wide)
#############################################################################
pdf(pdf_name)
kp <- plotKaryotype(genome = genome, chromosomes = chr)
kpAddBaseNumbers(kp, tick.dist = 1e7, add.units = TRUE)
colors <- c("cell-free DNA"="#2c4e60", "gDNA"="#5f9ea0")
kpPlotBAMDensity(kp, data=bam_cfDNA, col=colors["cell-free DNA"], r0=0, r1=0.5, window.size = 100000)
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density, r0=0, r1=0.5, cex=0.8)
kpAbline(kp, h=mean(kp$latest.plot$computed.values$density), lty=2, ymax=kp$latest.plot$computed.values$max.density, r0=0, r1=0.5)

kpPlotBAMDensity(kp, data=bam_gDNA, col=colors["gDNA"], r0=0.55, r1=1.05,  window.size = 100000)
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density, r0=0.55, r1=1.05,  cex=0.8)
kpAbline(kp, h=mean(kp$latest.plot$computed.values$density), lty=2, ymax=kp$latest.plot$computed.values$max.density, r0=0.55, r1=1.05)

legend(x = "right", legend = names(colors), fill = colors, border=NA)
dev.off()
# END OF SCRIPT
