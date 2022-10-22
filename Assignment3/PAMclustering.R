source("C:/Users/Stephen Smith/Documents/Intro-to-Bioinformatics-Team-15-Project/Human_PCA_plot.R")
library(cluster)
library(factoextra)

#for now, create a temp top5k so it doesn't interfere
#with the other cluster methods
PAM_top5k <- top5k

#scale the data
scaleddata = scale(PAM_top5k)
fviz_nbclust(scaleddata, pam, method ="silhouette")+theme_minimal()
    
#create a PAM clustering using k = 25
pamResult <- pam(scaleddata, k = 2)
pamResult

#bind the cluster data to the top5k
PAM_top5k$cluster = pamResult$clustering
head(PAM_top5k)

fviz_cluster(pamResult, 
             palette =c("#007892","#D9455F"),
             ellipse.type ="euclid",
             repel =TRUE,
             ggtheme =theme_minimal())







