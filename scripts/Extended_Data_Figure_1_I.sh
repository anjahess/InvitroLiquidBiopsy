#!/bin/bash/bin/bash
# This script requires deeptools installed.


##########################################################################
# 1. Compute matrix (GSE293866)
##########################################################################
#multiBigwigSummary bins -b *.bw -o results.npz

##########################################################################
# 2. Plot correlation (Pearson and Spearman)
##########################################################################

cd ../sourcedata/FIG_1EXT/E1I && \
echo $PWD && \
plotCorrelation -in results.npz \
    --corMethod pearson --skipZeros \
    --whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
    -o heatmap_pearson.pdf   \
    --outFileCorMatrix PearsonCorr_readCounts.tab && \
    
plotCorrelation -in results.npz \
    --corMethod spearman --skipZeros \
    --plotTitle "HEATMAP OF GLOBAL METHYLATION" \
    --whatToPlot heatmap --colorMap RdYlBu --plotNumbers \
    -o heatmap_spearman.pdf   \
    --outFileCorMatrix SpearmanCorr_readCounts.tab
