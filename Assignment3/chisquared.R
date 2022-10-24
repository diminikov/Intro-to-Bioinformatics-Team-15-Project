source("C:/Users/tnels/Documents/Intro-to-Bioinformatics-Team-15-Project/Assignment3/heat_map.R")

length(clusterOutput)
length(df)


dfM <- matrix(unlist(df))
coM <- matrix(unlist(clusterOutput))

dfM <- abs(dfM)
coM <- abs(coM)

chisq.test( dfM)  # human counts 
chisq.test( coM)  # cluster info
#chisq.test(dfM, coM)


