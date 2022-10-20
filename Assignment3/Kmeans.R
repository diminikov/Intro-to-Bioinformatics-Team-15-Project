source("C:/Users/tnels/Documents/Intro-to-Bioinformatics-Team-15-Project/Assignment3/Part2.R")

# Installing Packages
install.packages("ClusterR")
install.packages("cluster")

# Loading package
library(ClusterR)
library(cluster)

# Fitting K-Means clustering Model 
# to training dataset
set.seed(240) # Setting seed
Human_kmeans.re <- kmeans(top5k, centers = 9, nstart = 20)
Human_kmeans.re

# Cluster identification for 
# each observation
Human_kmeans.re$cluster

#view(Human_kmeans.re$cluster)
# Confusion Matrix
#colnames(top5k)[0] <- "genes"
#view(row.names(top5k))
Human_cm <- table(row.names(top5k), Human_kmeans.re$cluster)
Human_cm
