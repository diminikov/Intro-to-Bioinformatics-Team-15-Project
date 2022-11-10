#Part 3!!!
if (!("DESeq2" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("DESeq2", update = FALSE)
}
if (!("EnhancedVolcano" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("EnhancedVolcano", update = FALSE)
}
if (!("apeglm" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("apeglm", update = FALSE)
}
# Attach the DESeq2 library
library(DESeq2)

# Attach the ggplot2 library for plotting
library(ggplot2)

# We will need this so we can use the pipe: %>%
library(magrittr)
set.seed(12345)



#4.2
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

library("stringr")
HumanCounts = HumanCounts[,str_sort(names(HumanCounts), numeric = TRUE)]

all(rownames(HumanColdata) %in% colnames(HumanCounts))


library("DESeq2")
humandds <- DESeqDataSetFromMatrix(countData = round(HumanCounts),
                                   colData = HumanColdata,
                                   design = ~ characteristics_ch1)
deseq_object <- DESeq(humandds)
deseq_results <- results(deseq_object)
deseq_results <- lfcShrink(
  deseq_object, # The original DESeq2 object after running DESeq()
  coef = 3, # The log fold change coefficient used in DESeq(); the default is 2.
  res = deseq_results # The original DESeq2 results table
)
head(deseq_results)
deseq_df <- deseq_results %>%
  # make into data.frame
  as.data.frame() %>%
  # the gene names are row names -- let's make them a column for easy display
  tibble::rownames_to_column("Gene") %>%
  # add a column for significance threshold results
  dplyr::mutate(threshold = padj < 0.05) %>%
  # sort by statistic -- the highest values will be genes with
  # higher expression in RPL10 mutated samples
  dplyr::arrange(dplyr::desc(log2FoldChange))
head(deseq_df)
plotCounts(humandds, gene = "ENSG00000169059", intgroup = "characteristics_ch1")
human_volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  deseq_df,
  lab = deseq_df$Gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.01 # Loosen the cutoff since we supplied corrected p-values
)
human_volcano_plot


#MOUSE

library(readr)
GSE183548_series_matrix <- read_delim("GSE183548_series_matrix.txt", delim = "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE, skip = 28)
mouseColdata = data.frame(t(GSE183548_series_matrix[1:39,-1]))
colnames(mouseColdata) = sub("!Sample_", "", GSE183548_series_matrix$X1[1:39])
rownames(mouseColdata) = unlist(GSE183548_series_matrix[18,-1])
mouseColdata = mouseColdata[,-18]
mouseColdata = mouseColdata[-19:-27,]

GSE183548_Normalized_counts_3LL_samples <- read_delim("GSE183548_Normalized_counts_tumour_samples.txt", delim = "\t", escape_double = FALSE, trim_ws = TRUE)
mouseCounts = data.frame(GSE183548_Normalized_counts_3LL_samples)
rownames(mouseCounts) = GSE183548_Normalized_counts_3LL_samples$Row.names
mouseCounts = mouseCounts[,-1:-3]

mouseColdata<-mouseColdata[ order(rownames(mouseColdata)), ]

all(rownames(mouseColdata) %in% colnames(mouseCounts))

library("DESeq2")
mousedds <- DESeqDataSetFromMatrix(countData = round(mouseCounts),
                                   colData = mouseColdata,
                                   design = ~characteristics_ch1.2)
deseq_object <- DESeq(mousedds)
deseq_results <- results(deseq_object)
deseq_results <- lfcShrink(
  deseq_object, # The original DESeq2 object after running DESeq()
  coef = 3, # The log fold change coefficient used in DESeq(); the default is 2.
  res = deseq_results # The original DESeq2 results table
)
head(deseq_results)
deseq_df <- deseq_results %>%
  # make into data.frame
  as.data.frame() %>%
  # the gene names are row names -- let's make them a column for easy display
  tibble::rownames_to_column("Gene") %>%
  # add a column for significance threshold results
  dplyr::mutate(threshold = padj < 0.05) %>%
  # sort by statistic -- the highest values will be genes with
  # higher expression in RPL10 mutated samples
  dplyr::arrange(dplyr::desc(log2FoldChange))
head(deseq_df)
mouse_volcano_plot <- EnhancedVolcano::EnhancedVolcano(
  deseq_df,
  lab = deseq_df$Gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.01 # Loosen the cutoff since we supplied corrected p-values
)
mouse_volcano_plot


