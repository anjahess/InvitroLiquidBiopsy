suppressPackageStartupMessages(library(circlize))
#suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
#library(ComplexHeatmap)
library(circlize)
library(vioplot)

plotvio <- function(df, outdir, title, order, subvars) {
  cols = 5
  level_order <- order
  n_samples <- length(order)
  subvar_len <- length(subvars)
  df$feature <- df$"3"
  palette_blues <- colorRampPalette(colors = c("#bfcfcd","#2d435b"))(n_samples)
  #####################################################################
  # Draw split plot
  #####################################################################
  if (subvar_len == 2) {
    pdf(outdir, width=20, height=12)
    print(subvars)
    var_1 <- subvars[1]
    var_2 <- subvars[2]
    vioplot(value ~ variable, data = subset(df, feature == var_1),
            xlab = "", ylab = "Mean methylation", ylim = c(0,1),
            main = title, las = 2, plotCentre = "line",
            side = "left", col = "#bfcfcd")
    vioplot(value ~ variable,  data = subset(df, feature == var_2),
            xlab = "", ylab = "Mean methylation", ylim = c(0,1),
            add = TRUE, plotCentre = "line",
            side = "right", col = "#2d435b")
    legend("topright", fill = c("#bfcfcd","#2d435b"),
    legend = c(var_1, var_2), title = "Feature")
  }

  else {
    if (subvar_len > 2 && subvar_len < 20) {
        rows = ceiling(subvar_len/cols)
        pdf(outdir, width = cols*10, height = rows*8)
        par(mfrow=c(rows,cols))
        counter = 0
        for (x in subvars) {
            counter = counter + 1
            vioplot(value ~ variable,
            data = subset(df, feature == x),
            main=x,
            col = brewer.pal(10, "Paired")[counter],
            las = 2,  xlab = "", ylab = "Mean CpG methylation",
            level = level_order)
            }
        par(mfrow=c(1,1))
    }
    else {
        pdf(outdir, width=20, height=12)
        vioplot(value ~ variable, main=title, data = df,
            col = palette_blues,
            las = 2,  xlab = "", ylab = "Mean CpG methylation",
            level = level_order)
    }
    }
  ggsave(outdir, limitsize = FALSE) #A4 size in inches 7x9
  dev.off()
  }

plotsmoothscatter <- function(df, outdir, title, a, b) {
    head(df)
    print(a)
    print(b)
    pdf(outdir, width = 10, height = 5)
    layout(matrix(1:2, ncol = 2, byrow = TRUE))
    smoothScatter(df[,a], df[,b], xlab = a, ylab = b, xlim=c(0,1), ylim=c(0,1),
                  cex.axis = 1.2, cex.lab = 1.2, cex.sub=1.2, nrpoints = 0,
                  colramp=colorRampPalette(brewer.pal(9,"Blues")),
    main=title, useRaster=TRUE)
    abline(0.1, 1, lty = 2)
    abline(-0.1, 1, lty = 2)
    smoothScatter(df[,b], df[,a], xlab = b, ylab = a, xlim=c(0,1), ylim=c(0,1),
                  cex.axis = 1.2, cex.lab = 1.2, cex.sub=1.2, nrpoints = 0,
                  colramp=colorRampPalette(brewer.pal(9,"Blues")), useRaster=TRUE)
    abline(0.1, 1, lty = 2)
    abline(-0.1, 1, lty = 2)
    dev.off()
  }
  # END OF FUNCTION
