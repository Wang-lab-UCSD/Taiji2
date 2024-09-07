#### findOptimal function is to determine the best distance method and cluster numbers ####
#### arguments: ####
##>>>>>>> 
library(factoextra)
#### main body ####
findOptimal <- function(data_reduced, max_k=60){
  
  # elbow method
  fviz_nbclust(data_reduced, kmeans, method = "wss", k.max = max_k)
  
  # silhouette method
  g <- fviz_nbclust(data_reduced, kmeans, method = "silhouette",k.max = max_k)
  df2 <- g$data
  
  # use pearson distance
  res.dist <- get_dist(data_reduced, method = "pearson")
  g <- fviz_nbclust(data_reduced, kmeans, method = "silhouette",diss=res.dist,k.max = max_k)
  df2 <- rbind(df2,g$data)
  
  # use spearman distance
  res.dist <- get_dist(data_reduced, method = "spearman")
  g <- fviz_nbclust(data_reduced, kmeans, method = "silhouette",diss=res.dist,k.max = max_k)
  df2 <- rbind(df2,g$data)
  
  # use kendall distance
  res.dist <- get_dist(data_reduced, method = "kendall")
  g <- fviz_nbclust(data_reduced, kmeans, method = "silhouette",diss=res.dist,k.max = max_k)
  df2 <- rbind(df2,g$data)
  
  # use manhattan distance
  res.dist <- get_dist(data_reduced, method = "manhattan")
  g <- fviz_nbclust(data_reduced, kmeans, method = "silhouette",diss=res.dist,k.max = max_k)
  df2 <- rbind(df2,g$data)
  
  # visualization
  df2$method <- c(rep("euclidean",max_k),rep('pearson',max_k),rep('spearman',max_k),rep('kendall',max_k),rep('manhattan',max_k))
  df2[,'clusters'] <- as.numeric(as.character(df2[,'clusters']))
  df2 <- df2[df2$clusters>=1,]
  return(df2)

}

