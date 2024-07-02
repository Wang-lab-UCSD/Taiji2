#  bubble is a function to draw the bubble plot of specific TFs -----------------

# -------------- arguments: -----------------------------------------------------
#  Tcell: specific TF names
#  log_exp: log expression levels
#  rank: normalized pagerank scores
#  always generate the vertical plot

# -------------- main -----------------------------------------------------------
library(ggplot2)
library(RColorBrewer)
library(grDevices)
library(grid)
library(ggpubr)

bubble <- function(Tcell,log_exp,rank,
                   angle=90, show.legend=T,
                   bubble.color = c("blue","white","red"),
                   plot.margin=c(4,4.5,0.5,1),med=0.5,
                   text.size=4,
                   color.title="Normalized rank score"){
  samples <- colnames(rank)
  
  # map pagerank scores to colors-----------------------------------
  # map2color<-function(x,pal,limits=NULL){
  #   if(is.null(limits)) limits=range(x)
  #   pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
  # }
  # mypal <- colorRampPalette(bubble.color)(20)
  # mypal <- map2color(as.numeric(unlist(rank)),mypal)
  a <- as.numeric(unlist(log_exp))
  
  # the default is vertical plot, we make sure height > width, i.e. y > x
  # if (nrow(log_exp) < ncol(log_exp)){ 
  #   df <- expand.grid(x=seq(nrow(log_exp)),y=seq(ncol(log_exp)))
  #   plotMargin <- c(2,8,1,1) # top, right, bottom, left. Tcell on top, samples on right
  # }else{
  #   df <- expand.grid(y=seq(nrow(log_exp)),x=seq(ncol(log_exp)))
  #   plotMargin <- c(3.5,4.5,1,1)
  # }
  
  # set plot grid and plot margin---------------------------------
  df <- expand.grid(y=seq(nrow(log_exp)),x=seq(ncol(log_exp)))
  df <- data.frame(df,"log expression level" = a)
  
  # ggplot -------------------------------------------------------
  color <- as.numeric(unlist(rank))
  p <- ggplot(data=df) + aes(x=x, y=y, size=a) + 
    geom_point(aes(color=color))+
    # scale_color_gradient2(low = bubble.color[1], 
    #                       mid = bubble.color[2], 
    #                       high = bubble.color[3], 
    #                       midpoint = quantile(color,med))+
    scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100))+
    theme_classic()+theme(plot.margin = unit(plot.margin, "lines"),
                          axis.line = element_blank(),axis.text = element_blank(),
                          axis.ticks = element_blank(),axis.title = element_blank(),
                          legend.position = "left",legend.title = element_text(size = text.size),
                          legend.text = element_text(size = text.size))+
    labs(size = "Log expression level",color = color.title)

  for (i in length(Tcell):1)  {
    p <- p + annotation_custom(
      grob = text_grob(label = Tcell[i], hjust = 0, size = text.size, face = "italic"),
      ymin = i,      # Vertical position of the textGrob
      ymax = i,
      xmin = length(samples)+0.5,         # Note: The grobs are positioned outside the plot area
      xmax = length(samples)+0.5)
  }
  for (i in length(samples):1)  {
    p <- p + annotation_custom(
      grob = text_grob(label = samples[i], hjust = 0, size = text.size, rot = angle),
      ymin = length(Tcell)+0.5,
      ymax = length(Tcell)+0.5,
      xmin = i,
      xmax = i)
  }
  
  if (show.legend==F){
    p <- p + theme(legend.position = "none")
  }
  # Code to override clipping-----------------------------
  gt <- ggplot_gtable(ggplot_build(p))
  gt$layout$clip[gt$layout$name == "panel"] <- "off"
  return(gt)
}
