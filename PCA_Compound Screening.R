# PCA analysis of chemical structures to check the variances based on descriptors

# Loading libraries for reading data frame, PCA analysis, and visualisation

library(factoextra)
library(FactoMineR)
library(tidyr)
library(dplyr)
library(cluster)

# Setting up your working directory
setwd("C:/Users/arnab/xx/xx/xx/")

# Loading/reading the data frame
df <- read.csv("C:/Users/arnab/xx/xx/data.csv", header = T, sep = ",", row.names = 1)

#summary of the object (df)
str(df)

#PCA analysis
pca <- prcomp(df, scale. = T)

#Summary of the PCA. Determining the variances in terms of dimensions/Principle components
summary(pca)

# Measuring the covariance of the variables applied to eigenvectors to determine the magnitude in contribution
eig.val <- get_eigenvalue(pca)
eig.val

# Scree plot for determing the variance percentage in each dimensions
fviz_eig(pca, addlabels = TRUE, ylim = c( 0, 50))

# Contribution plot of the individual compounds
fviz_pca_ind(pca,
             col.ind = "cos2", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     
)

# Biplot
fviz_pca_biplot(pca, repel = TRUE,
                col.var = "#2E9FDF",
                col.ind = "#696969"
)

# Factor map
fviz_pca_var(pca,
             col.var = "contrib", # representation by contributions to PCs
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE    
)

#Laveraging dimensions of the PCA data
pcacluster <- as.data.frame(-pca$x[,1:4])
head(pcacluster)

# Optimum Number of Clusters using WSS method
fviz_nbclust(pcacluster, kmeans, method = 'wss')

#Silhouette method (Validation of optimum number of clusters)
fviz_nbclust(df, kmeans, "silhouette", k.max = 15) + labs(subtitle = "Silhouette method")

# k-Means clustering
set.seed(123)
cluster_pca <- kmeans(pcacluster, centers = 9, nstart = 25, iter.max = 1000)

#summary
print(cluster_pca)
data1 <- cluster_pca$cluster

#save your clusters in .csv
write.csv(data1, file="CLUSTER2.csv")

# k-Means cluster plot
fviz_cluster(cluster_pca, data = pcacluster, repel = TRUE)

#cluster details
pc_cluster_2 <-kmeans(pcacluster, 9)
pc_cluster_2$cluster
pc_cluster_2$centers
pc_cluster_2$size
pc_cluster_2$totss
pc_cluster_2$withinss
pc_cluster_2$betweenss
pc_cluster_2$iter

#Arnab- arnabbiotech.gen@gmail.com