source("C:/Users/Stephen Smith/Documents/Intro-to-Bioinformatics-Team-15-Project/Human_PCA_plot.R")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")

#convert out cluster data frame into a matrix
pamResultMatrix <- data.matrix(PAM_top5k)

heatmap(pamResultMatrix)