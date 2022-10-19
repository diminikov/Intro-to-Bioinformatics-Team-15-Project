source("C:/Users/tnels/Documents/Intro-to-Bioinformatics-Team-15-Project/Human_PCA_plot.R")
#order human counts by variation
var_genes <- apply(HumanCounts, 1, var)
head(var_genes)
#sort by decreasing and get gene names of the top 5000 rows
select_var <- names(sort(var_genes, decreasing = TRUE))[1:5000]
head(select_var)
top5k <- HumanCounts[select_var,]
view(top5k)


