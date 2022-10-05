library(readr)
GSE183548_series_matrix <- read_delim("GSE183548_series_matrix.txt", delim = "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE, skip = 28)
mouseColdata = data.frame(t(GSE183548_series_matrix[1:39,-1]))
colnames(mouseColdata) = sub("!Sample_", "", GSE183548_series_matrix$X1[1:39])
rownames(mouseColdata) = unlist(GSE183548_series_matrix[18,-1])
mouseColdata = mouseColdata[,-18]
mouseColdata = mouseColdata[-1:-18,]

GSE183548_Normalized_counts_3LL_samples <- read_delim("GSE183548_Normalized_counts_3LL_samples.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
mouseCounts = data.frame(GSE183548_Normalized_counts_3LL_samples)
rownames(mouseCounts) = GSE183548_Normalized_counts_3LL_samples$Row.names
mouseCounts = mouseCounts[,-1:-3]

all(rownames(mouseColdata) %in% colnames(mouseCounts))

library("DESeq2")
mousedds <- DESeqDataSetFromMatrix(countData = round(mouseCounts),
                              colData = mouseColdata,
                              design = ~ title)
mousevsd <- vst(mousedds)
#rld <- rlog(dds)
plotPCA(mousevsd, intgroup=c("title"))

plot.mouse <- function(x, labels,
         main="A UMAP visualization of the Iris dataset",
         colors=c("#ff7f00", "#e377c2", "#17becf"),
         pad=0.1, cex=0.6, pch=19, add=FALSE, legend.suffix="",
         cex.main=1, cex.legend=0.85) {

  layout <- x
  if (is(x, "umap")) {
    layout <- x$layout
  }

  xylim <- range(layout)
  xylim <- xylim + ((xylim[2]-xylim[1])*pad)*c(-0.5, 0.5)
  if (!add) {
    par(mar=c(0.2,0.7,1.2,0.7), ps=10)
    plot(xylim, xylim, type="n", axes=F, frame=F)
    rect(xylim[1], xylim[1], xylim[2], xylim[2], border="#aaaaaa", lwd=0.25)
  }
  points(layout[,1], layout[,2], col=colors[as.integer(labels)],
         cex=cex, pch=pch)
  mtext(side=3, main, cex=cex.main)

  labels.u <- unique(labels)
  legend.pos <- "topleft"
  legend.text <- as.character(labels.u)
  if (add) {
    legend.pos <- "bottomleft"
    legend.text <- paste(as.character(labels.u), legend.suffix)
  }

  legend(legend.pos, legend=legend.text, inset=0.03,
         col=colors[as.integer(labels.u)],
         bty="n", pch=pch, cex=cex.legend)
}

library(umap)
mouse.umap <- umap(mouseCounts)
plot.mouse(mouse.umap, mouseColdata[1:2,1])
