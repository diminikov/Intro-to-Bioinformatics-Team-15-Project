source("C:/Users/Stephen Smith/Documents/Intro-to-Bioinformatics-Team-15-Project/Human_PCA_plot.R")
library(cluster)
library(factoextra)

scaleddata = scale(top5k)
fviz_nbclust(scaleddata, pam, method ="silhouette")+theme_minimal()


pamResult <- pam(top5k, k = 25)
pamResult




