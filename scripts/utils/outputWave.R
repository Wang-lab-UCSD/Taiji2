outputWave <- function(x, df, df2, Da, Group){
  # x: cluster name
  # df: wave plot layout
  # df2: clusters of TF result
  # Da: normalized PageRank score
  # Group: group of samples 

  geneL <- rownames(subset(df2, cluster == x))
  DT <- as.data.table(t(as.matrix(Da[geneL,])))
  DT <- data.table(m = rowMeans(DT))
  DT$group <- Group$Group
  DT <- DT[,.(mm=mean(m)),by=group] 
  DT <- DT[match(df$samplename, DT$group),]
  
  rank = DT$mm
  df$color = range01(as.numeric(unlist(rank)))
  df$fill = range01(as.numeric(unlist(rank)))
  
  # png(filename = paste0("c",x,"_",length(geneL),".png"),units="in", width=5, height=5, res=300)
  pdf(paste0("c",x,".pdf"), w=5, h=5)
  p <- wave(df)
  print(p)
  dev.off()
  writeLines(geneL,paste0("c",x,".txt"))
  return(p)
}
