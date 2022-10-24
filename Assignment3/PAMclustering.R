source("C:/Users/Stephen Smith/Documents/Intro-to-Bioinformatics-Team-15-Project/Human_PCA_plot.R")
library(cluster)
library(factoextra)

#for now, create a temp top5k so it doesn't interfere
#with the other cluster methods
PAM_top5k <- top5k

#scale the data
scaleddata = scale(PAM_top5k)

#scale data for top 10, 100, and 1000 genes
scaleddata_top10 <- scale(top10)
scaleddata_top100 <- scale(top100)
scaleddata_top1k <- scale(top1k)
scaleddata_top10k <- scale(top10k)

 
fviz_nbclust(scaleddata, pam, method ="silhouette")+theme_minimal()
    
#create a PAM clustering using k = 6
pamResult <- pam(scaleddata, k = 6)
pamResult

#creating a PAM cluster with top10, top100, and top1000 genes
pamResult_top10 <- pam(scaleddata_top10, k = 6)
pamResult_top100 <- pam(scaleddata_top100, k = 6)
pamResult_top1k <- pam(scaleddata_top1k, k = 6)
pamResult_top10k <- pam(scaleddata_top10k, k = 6)

#bind the cluster data to the top5k
PAM_top5k$cluster = pamResult$clustering
head(PAM_top5k)

fviz_cluster(pamResult, 
             palette =c("#DAF7A6","#FFC300", "#FF5733", "#C70039", "#900C3F", "#581845"),
             ellipse.type ="euclid",
             repel =TRUE,
             ggtheme =theme_minimal())







