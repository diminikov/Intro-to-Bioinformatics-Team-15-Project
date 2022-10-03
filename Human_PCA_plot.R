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

mouseCounts

dds <- DESeqDataSetFromMatrix(countData = round(mouseCounts),
                              colData = mouseColdata,
                              design = ~ title)
dds

vsd <- vst(dds)
rld <- rlog(dds)
plotPCA(vsd, intgroup=c("title"))