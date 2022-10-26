
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library(ClusterR)
library(cluster)
BiocManager::install("ConsensusClusterPlus")
d <- data.matrix(top5k,rownames.force = NA)
d = sweep(d,1,apply(d, 1, median,na.rm=T))
library(ConsensusClusterPlus)

title=tempdir()
results = ConsensusClusterPlus(d,maxK=25,reps=1000,pItem=0.8,pFeature=1,title=title,clusterAlg="hc",distance="pearson",seed=1262118388.71279,plot="png")
results[[2]][["consensusMatrix"]]]
#consensusTree - hclust object
results[[2]][["consensusTree"]]
hclust(d = as.dist(1 - fm), method = finalLinkage)

icl = calcICL(results,title=title,plot="png")
results <- data.matrix(icl[["clusterConsensus"]])

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")

#convert out cluster data frame into a matrix
ResultMatrixCCP <- data.matrix(results)
view(ResultMatrixCCP)

heatmap(ResultMatrixCCP)




