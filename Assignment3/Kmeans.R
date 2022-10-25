source("C:/Users/tnels/Documents/Intro-to-Bioinformatics-Team-15-Project/Assignment3/Part2.R")

# Installing Packages
install.packages("ClusterR")
install.packages("cluster")

# Loading package
library(ClusterR)
library(cluster)

head(top5k)
df <- scale(top5k)
head(df, n = 3)

install.packages("factoextra")
library(factoextra)


fviz_nbclust(df, kmeans, method = "wss") +
  geom_vline(xintercept = 4, linetype = 2)+
  labs(subtitle = "Elbow method")

# Fitting K-Means clustering Model 
# to training dataset
set.seed(123) # Setting seed
Human_kmeans.re <- kmeans(df, 5, nstart = 25)
print(Human_kmeans.re)




