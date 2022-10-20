if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ConsensusClusterPlus")
d <- data.matrix(top5k,rownames.force = NA)
d = sweep(d,1,apply(d, 1, median,na.rm=T))
library(ConsensusClusterPlus)
title=tempdir()
results = ConsensusClusterPlus(d,maxK=25,reps=1000,pItem=0.8,pFeature=1,title=title,clusterAlg="hc",distance="pearson",seed=1262118388.71279,plot="png")
icl = calcICL(results,title=title,plot="png")
icl[["clusterConsensus"]]