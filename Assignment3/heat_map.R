source("C:/Users/tnels/Documents/Intro-to-Bioinformatics-Team-15-Project/Assignment3/Kmeans.R")

# Heatmap
# Human_kmeans.re is kmeansObj
# datamatrix is df
#par(mfrom = c(1,10))
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")
library(ComplexHeatMap)
hm <- heatmap(df, show_row_names = FALSE)
head(df)

ordered <- Human_kmeans.re[order(Human_kmeans.re$cluster),1]

view(ordered)


clustMatrix <- data.matrix(Human_kmeans.re$cluster) 
view(clustMatrix)
view(Human_kmeans.re$cluster)
hm <- heatmap(df, reorderfun = Human_kmeans.re$cluster, show_row_names = FALSE)


