plotPCA <- function(x){
  data.pca <- prcomp(x)
  
  # plot the cumulative variance explained vs # of components
  s = (data.pca$sdev)^2/sum((data.pca$sdev)^2)
  df2 <- data.frame(x=c(1:length(data.pca$sdev)),y=cumsum(s))
  p2 <- ggplot(df2)+aes(x=x,y=y)+geom_point()+geom_line()+
    labs(x="Principal Component",y="Cumulative Proportion of Variance Explained")+
    theme(axis.title = element_text(size = 15),
          axis.text  = element_text(size = 15))
  png("pca.png",units="in", width=4, height=5, res=300)
  print(p2)
  dev.off()
}
