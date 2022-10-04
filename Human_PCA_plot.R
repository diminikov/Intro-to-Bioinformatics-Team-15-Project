library(readr)
GSE183516_series_matrix <- read_delim("GSE183516_series_matrix.txt", delim = "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE, skip = 28)
HumanColdata = data.frame(t(GSE183516_series_matrix[1:39,-1]))
colnames(HumanColdata) = sub("!Sample_", "", GSE183516_series_matrix$X1[1:39])
rownames(HumanColdata) = unlist(GSE183516_series_matrix[19,-1])
#HumanColdata = HumanColdata[,-1]
HumanColdata = HumanColdata[-37:-60,]


GSE183516_Normalized_counts_H23_H358_A549_cell_lines <- read_delim("GSE183516_Normalized_counts_H23_H358_A549_cell_lines.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
HumanCounts<-GSE183516_Normalized_counts_H23_H358_A549_cell_lines[,order(colnames(GSE183516_Normalized_counts_H23_H358_A549_cell_lines))]
GSE183516_Normalized_counts_H23_H358_A549_cell_lines[ , order(colnames(GSE183516_Normalized_counts_H23_H358_A549_cell_lines))]

HumanCounts = data.frame(GSE183516_Normalized_counts_H23_H358_A549_cell_lines)
rownames(HumanCounts) = GSE183516_Normalized_counts_H23_H358_A549_cell_lines$...1
HumanCounts = HumanCounts[,-1]

all(rownames(HumanColdata) %in% colnames(HumanCounts))

HumanCounts

dds <- DESeqDataSetFromMatrix(countData = round(HumanCounts),
                              colData = HumanColdata,
                              design = ~ title)
dds

vsd <- vst(dds)
rld <- rlog(dds)
plotPCA(vsd, intgroup=c("title"))