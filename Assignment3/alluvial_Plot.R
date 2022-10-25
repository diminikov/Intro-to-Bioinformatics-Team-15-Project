source("C:/Users/Stephen Smith/Documents/Intro-to-Bioinformatics-Team-15-Project/Human_PCA_plot.R")

library(highcharter)
library(htmlwidgets)
library(dplyr)

#create a dataframe that combines the cluster columns from the top 10 - 10000 genes
f10 <- pamResult_top10$clustering
f100 <- pamResult_top100$clustering[1:10]
f1k <- pamResult_top1k$clustering[1:10]
f10k <- pamResult_top10k$clustering[1:10]
cluster_df <- data.frame(cbind(f10, f100, f1k, f10k))

#creating Sankey diagram
hcart(data_to_sankey(cluster_df), "sankey", name = "Cluster Results based on number of Genes")

