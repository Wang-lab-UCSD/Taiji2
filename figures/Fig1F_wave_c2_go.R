library(ggplot2)
library(xlsx)
library(dplyr)

setwd("~/Documents/allcombn/differentiationWave/wave7/")
df <- read.xlsx("c2_GO_KEGG.xlsx", sheetIndex = 1) %>% filter(ID %in% c("GO:0045165","mmu05235","mmu05321","mmu05161","mmu05169")) 
p1 <- ggplot(df, aes(x=Count,y=Description, fill=-log10(p.adjust))) + geom_bar(stat = "identity") +
      # scale_fill_gradientn(limits = c(0,0.05), colours = c("red","blue"), name = "-log10(p.adjust)") +
      scale_fill_gradient(high = "red", low = "blue")+
      theme_classic() + xlab("count") + ylab("") 
pdf(file = paste0("c2_selected_GO_KEGG.pdf"), height = 2)
print(p1)
dev.off()

p1 <- ggplot(df, aes(x=-log10(p.adjust),y=Description, fill=-log10(p.adjust))) + geom_bar(stat = "identity") +
  # scale_fill_gradientn(limits = c(0,0.05), colours = c("red","blue"), name = "-log10(p.adjust)") +
  scale_fill_gradient(high = "red", low = "blue")+
  theme_classic() + xlab("-log10(p.adjust)") + ylab("") 
pdf(file = paste0("c2_selected_GO_KEGG_v2.pdf"), height = 2)
print(p1)
dev.off()

#### plot for other transcriptional waves ----------------
meta <- read.xlsx("../../fromKay/SuppFig4 wave GSEA pathways.xlsx", sheetIndex = 1)
lapply(unique(meta$wave.cluster), function(x){
  df <- meta %>% filter(wave.cluster==x) 
  p1 <- ggplot(df, aes(x=-log10(p.adjust),y=Description, fill=-log10(p.adjust))) + geom_bar(stat = "identity") +
    # scale_fill_gradientn(limits = c(0,0.05), colours = c("red","blue"), name = "-log10(p.adjust)") +
    scale_fill_gradient(high = "red", low = "blue")+
    theme_classic() + xlab("-log10(p.adjust)") + ylab("") 
  pdf(file = paste0(x,"_selected_GO_KEGG.pdf"), height = 3)
  print(p1)
  dev.off()
  
})

#### plot the pathways all together ----------------
library(pheatmap)
df <- read.xlsx("../../fromKay/SuppFig4 wave GSEA pathways.xlsx", sheetIndex = 1) %>% dplyr::select(wave.cluster, Description, p.adjust) %>%
      tidyr::pivot_wider(names_from = wave.cluster, values_from = p.adjust) %>% as.data.frame() %>% replace(is.na(.), 1) %>%
      tibble::column_to_rownames("Description") 

get_heatmap <- function(dt, pathway_name="all", clustering_method="hclust", cutoff=0.05, k=5){
  
  ## set cut-off as p-value max
  dt2 <- apply(dt,1:2, function(x) min(x, cutoff)) %>% as.data.frame()
  
  ## remove the rows with no variation
  dt3 <- dt2[rowSums(dt2)!=cutoff*ncol(dt2),]
  print(dim(dt3))
  ## heatmap
  if (clustering_method=="hclust"){
    p1 <- pheatmap(dt3, fontsize = 7, show_rownames = T,
                   angle_col = 45, show_colnames = T,
                   cellwidth = 8, cellheight = 8,
                   cluster_rows = T, cluster_cols = T,
                   clustering_distance_rows = "correlation",
                   clustering_distance_cols = "correlation",
                   clustering_method = "average",
                   # annotation_col = annotation_col,
                   # annotation_colors = mycolors, 
                   border_color = NA,
                   color = colorRampPalette(c("red","blue"))(20),
                   filename = paste0("hp_transcriptional_waves_",pathway_name,"_pathway_summary_",cutoff,"_max_",clustering_method,"_cluster.pdf"))
    
  }else if(clustering_method=="kmeans"){
    prefix <- paste0("hp_transcriptional_waves_",pathway_name,"_pathway_summary_",cutoff,"_max_",clustering_method,"_",k,"_cluster")
    p <- pheatmap(dt3, kmeans_k = k,
                  fontsize = 7, show_rownames = T,
                  angle_col = 45, show_colnames = T,
                  cellwidth = 8, cellheight = 8,
                  cluster_rows = T, cluster_cols = T,
                  clustering_distance_rows = "correlation",
                  clustering_distance_cols = "correlation",
                  clustering_method = "average",
                  # annotation_col = annotation_col,
                  # annotation_colors = mycolors, 
                  border_color = NA,
                  color = colorRampPalette(c("red","blue"))(20),
                  filename = paste0(prefix,".pdf"))
    gene.clust <- data.frame(cluster=p$kmeans$cluster) %>% tibble::rownames_to_column(var = "pathway") %>%
      arrange(cluster, pathway)
    write.csv(gene.clust,paste0(prefix,"_pathwayClusters.txt"), quote=F, row.names=F)
  }else if(clustering_method=="none"){
    p1 <- pheatmap(dt3, fontsize = 7, show_rownames = T,
                   angle_col = 45, show_colnames = T,
                   cellwidth = 8, cellheight = 8,
                   cluster_rows = F, cluster_cols = T,
                   clustering_distance_rows = "correlation",
                   clustering_distance_cols = "correlation",
                   clustering_method = "average",
                   # annotation_row = annotation_row,
                   # annotation_colors = mycolors,
                   border_color = NA,
                   color = colorRampPalette(c("red","blue"))(20),
                   filename = paste0("hp_transcriptional_waves_",pathway_name,"_pathway_summary_",cutoff,"_max_",clustering_method,"_cluster.pdf"))
    
  }
}
get_heatmap(df, clustering_method = "hclust")
get_heatmap(df, clustering_method = "none")
