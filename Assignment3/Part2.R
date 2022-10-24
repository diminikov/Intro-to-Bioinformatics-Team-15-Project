source("C:/Users/Stephen Smith/Documents/Intro-to-Bioinformatics-Team-15-Project/Human_PCA_plot.R")
#order human counts by variation
var_genes <- apply(HumanCounts, 1, var)
#head(var_genes)
#sort by decreasing and get gene names of the top 5000 rows
select_var <- names(sort(var_genes, decreasing = TRUE))[1:5000]
#head(select_var)
# use top5k for cluster algorithm 
top5k <- HumanCounts[select_var,]
#view(top5k)

#rerunning cluster method with 10, 100, 1000, and 10000 genes
select_var <- names(sort(var_genes, decreasing = TRUE))[1:10]
top10 <- HumanCounts[select_var,]

select_var <- names(sort(var_genes, decreasing = TRUE))[1:100]
top100 <- HumanCounts[select_var,]

select_var <- names(sort(var_genes, decreasing = TRUE))[1:1000]
top1k <- HumanCounts[select_var,]

select_var <- names(sort(var_genes, decreasing = TRUE))[1:10000]
top10k <- HumanCounts[select_var,]

