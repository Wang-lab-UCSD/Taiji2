wave <- function(df){
  stroke = 2; 
  p <- ggplot(data=df) + aes(x=x, y=y,fill = fill,color=color) + 
    geom_link2(data=df[c(1,2),],group=1, size = stroke, alpha = 0.6,
               lineend = 'round', n = 500) + 
    geom_link2(data=df[c(1,3,5,7),],group=1, size = stroke, alpha = 0.6,
               lineend = 'round', n = 500) + 
    geom_link2(data=df[c(1,4,6),],group=1, size = stroke, alpha = 0.6,
               lineend = 'round', n = 500) +
    geom_link2(data=df[c(4,8),],group=1, size = stroke, alpha = 0.6,
               lineend = 'round', n = 500) + 
    geom_link2(data=df[c(4,9),],group=1, size = stroke, alpha = 0.6,
               lineend = 'round', n = 500) + 
    geom_link2(data=df[c(3,7),],group=1, size = stroke, alpha = 0.6,
               lineend = 'round', n = 500) +
    geom_link2(data=df[c(6,8),],group=1, size = stroke, alpha = 0.6,
               lineend = 'round', n = 500) +
    geom_link2(data=df[c(6,9),],group=1, size = stroke, alpha = 0.6,
               lineend = 'round', n = 500) +
    geom_link2(data=df[c(8,9),],group=1, size = stroke, alpha = 0.6,
               lineend = 'round', n = 500) +
    scale_color_gradient(low = "white", high = "red", guide = "none")+
    geom_point(shape = 21, alpha=0.6, size=5, stroke = 1, color = "gray")+
    scale_fill_gradient(low = "white", high = "red")+ 
    geom_text(aes(label=samplename,x=labelposx, y = labelposy),color = 'black',size = 4, hjust=0, check_overlap = TRUE)+
    # geom_text_repel(aes(label=samplename),color='black',size=2)+
    theme_classic()+theme(plot.margin = unit(c(1,1,1,1), "lines"),
                          axis.line = element_blank(),axis.text = element_blank(),
                          axis.ticks = element_blank(),axis.title = element_blank(),
                          legend.position = "left",legend.title = element_text(size = 6),
                          legend.text = element_text(size=6))+
    labs(fill = "Normalized rank score") + coord_equal() + ylim(c(-1,6)) + xlim(c(0,7))
  return(p)
}
