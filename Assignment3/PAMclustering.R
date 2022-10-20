source("C:/Users/Stephen Smith/Documents/Intro-to-Bioinformatics-Team-15-Project/Human_PCA_plot.R")

library(cluster)
library(factoextra)

pamResult <- pam(top5k, k = 2)
pamResult
#pam(x, k, metric="euclidean", stand=FALSE)




