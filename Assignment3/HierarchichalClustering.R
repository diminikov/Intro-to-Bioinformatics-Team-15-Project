source("~/Intro-to-Bioinformatics-Team-15-Project/Assignment3/Part2.R")

top5k_scaled <- as.data.frame(scale(top5k))
top1k_scaled <- as.data.frame(scale(top1k))
top100_scaled <- as.data.frame(scale(top100))
top10_scaled <- as.data.frame(scale(top10))

dist_mat5k <- dist(top5k_scaled, method = 'euclidean')
hclust_avg5k <- hclust(dist_mat5k, method = 'average')
dist_mat1k <- dist(top1k_scaled, method = 'euclidean')
hclust_avg1k <- hclust(dist_mat1k, method = 'average')
dist_mat100 <- dist(top100_scaled, method = 'euclidean')
hclust_avg100 <- hclust(dist_mat100, method = 'average')
dist_mat10 <- dist(top10_scaled, method = 'euclidean')
hclust_avg10 <- hclust(dist_mat10, method = 'average')

cut_avg5k <- cutree(hclust_avg5k, k = 36)
cut_avg1k <- cutree(hclust_avg1k, k = 36)
cut_avg100 <- cutree(hclust_avg100, k = 36)
cut_avg10 <- cutree(hclust_avg10, k = 3)

suppressPackageStartupMessages(library(dendextend))
avg_dend_obj5k <- as.dendrogram(hclust_avg5k)
avg_col_dend5k <- color_branches(avg_dend_obj5k, h = 36)
plot(avg_col_dend5k)

suppressPackageStartupMessages(library(dplyr))
top5k_cl <- mutate(top5k, cluster = cut_avg5k)
count(top5k_cl,cluster)


BiocManager::install("ComplexHeatmap")

hclustMatrix <- data.matrix(top5k_cl)

chisq.test(unlist(top5k_cl$cluster), simulate.p.value = TRUE)

heatmap(hclustMatrix)
