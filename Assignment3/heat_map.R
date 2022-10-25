source("C:/Users/tnels/Documents/Intro-to-Bioinformatics-Team-15-Project/Assignment3/Kmeans.R")

# Heatmap
# Human_kmeans.re is kmeansObj
# datamatrix is df
#par(mfrom = c(1,10))
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")
library(ComplexHeatMap)


set.seed(100)
m = matrix(rnorm(10), 100, 5)
km = kmeans(m, 10)
m2 <- cbind(m,km$cluster)
view(m2)

ordered <- Human_kmeans.re[order(Human_kmeans.re$cluster),1]

view(M1)

# M1 = km
# m = df

m2 <- cbind(df, M1$cluster)
view(m2)
o <- order(m2[,37])
view(o)
m3 <- m2[o,]
view(m3)

#heatmap(df)
#heatmap(m3)
heatmap(m3[,1:36])

clusterOutput <- m3[,1:36]



