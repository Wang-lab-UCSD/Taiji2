library(ggplot2)
library(xlsx)
library(dplyr)
library(readr)
library(RColorBrewer)
library(data.table)
library(pheatmap)
maindir <- '~/allcombn/'
setwd(maindir)
source(paste0(maindir, 'scripts/functions/utils.R'))
source(paste0(maindir, 'scripts/functions/firstUp.R'))
source(paste0(maindir,"scripts/functions/enrich_all.R"))

# genesets <- lapply(paste0('c',1:5),function(x) read.xlsx('Fig4. network GSEA analysis.xlsx', sheetName = x, header = F) %>% na.omit())
# tex_file <- read.xlsx(paste0(maindir,'network/community/TexTerm/clust_5_v2/community_regulatees_edge_weight_top2000_p.adj_1e-07_enrich.xlsx'), sheetIndex = 1) %>% mutate(cluster=paste0(cluster,'_TEX'))
# trm_file <- read.xlsx(paste0(maindir,'network/community/TRM/clust_5/community_regulatees_edge_weight_top2000_p.adj_1e-07_enrich.xlsx'), sheetIndex = 1) %>% mutate(cluster=paste0(cluster,'_TRM'))
tex_file <- read.xlsx(paste0(maindir,'network/community/TexTerm/clust_5_v2/community_regulatees_edge_weight_top2000_p.adj_1e-10_enrich.xlsx'), sheetIndex = 1) %>% mutate(cluster=paste0(cluster,'_TEX'))
trm_file <- read.xlsx(paste0(maindir,'network/community/TRM/clust_5/community_regulatees_edge_weight_top2000_p.adj_1e-10_enrich.xlsx'), sheetIndex = 1) %>% mutate(cluster=paste0(cluster,'_TRM'))
df <- rbind(tex_file, trm_file) %>% mutate(Description = gsub(' $','',Description))
df2 = df %>% filter(Description %in% ps)
df2 = df %>% filter(Description %in% ps) %>% dplyr::select(Description, p.adjust, cluster)
  
get_heatmap <- function(x){
  ps = genesets[[x]]$X1
  df2 = df %>% filter(Description %in% ps) %>% dplyr::select(Description, p.adjust, cluster) %>% 
    tidyr::pivot_wider(names_from = cluster, values_from = p.adjust, values_fill=1) %>% tibble::column_to_rownames('Description')
  df2 = -log10(df2[intersect(ps, rownames(df2)),names(df2)[order(names(df2))]])
  p1 <- pheatmap(df2, show_rownames = T,
                 angle_col = 90, show_colnames = T,
                 cellwidth = 8, cellheight = 8,
                 fontsize = 8, 
                 cluster_rows = F, cluster_cols = F,
                 clustering_distance_rows = 'correlation',
                 clustering_method = 'average',
                 treeheight_row = 2, treeheight_col = 2,
                 border_color = NA,
                 color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(20),
                 filename = paste0("fig4_hp_community_",x,"_regulatees_ew_2000_p_1e-10_pathways.pdf")
  )
}
lapply(1:5, get_heatmap)


# combined heatmap
genesets <- read.xlsx('Fig4. network GSEA analysis_20240229.xlsx', sheetName = 'combined', header = F) %>% na.omit()
names(genesets) <- c('genesets','categories')
category_colors <- read.xlsx('Fig4. network GSEA analysis_20240229.xlsx', sheetName = 'categories', header = F)
tex_file <- read.xlsx(paste0(maindir,'network/community/TexTerm/clust_5_v2/community_regulatees_edge_weight_top2000_p.adj_1e-07_enrich.xlsx'), sheetIndex = 1) %>% mutate(cluster=paste0(cluster,'_TEX'))
trm_file <- read.xlsx(paste0(maindir,'network/community/TRM/clust_5/community_regulatees_edge_weight_top2000_p.adj_1e-07_enrich.xlsx'), sheetIndex = 1) %>% mutate(cluster=paste0(cluster,'_TRM'))
df <- rbind(tex_file, trm_file) %>% mutate(Description = gsub(' $','',Description))
ps <- genesets$genesets
df2 = df %>% filter(Description %in% ps)
df2 = df %>% filter(Description %in% ps) %>% dplyr::select(Description, p.adjust, cluster) %>% 
  tidyr::pivot_wider(names_from = cluster, values_from = p.adjust, values_fill=1) %>% tibble::column_to_rownames('Description')
df2 = -log10(df2[intersect(ps, rownames(df2)),names(df2)[order(names(df2))]])

# set annotation
annotation_row <- genesets %>% tibble::column_to_rownames('genesets')
mycolors <- list(categories = category_colors$X2)
names(mycolors$categories) = category_colors$X1

p1 <- pheatmap(df2, show_rownames = T,
               angle_col = 90, show_colnames = T,
               cellwidth = 8, cellheight = 8,
               fontsize = 8, 
               cluster_rows = T, cluster_cols = F,
               clustering_distance_rows = 'correlation',
               clustering_method = 'average',
               treeheight_row = 10, treeheight_col = 2,
               border_color = NA,
               annotation_row = annotation_row, annotation_colors = mycolors,
               color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdBu")))(20),
               filename = paste0("fig4_combined_hp_community_regulatees_ew_2000_p_1e-07_pathways.pdf")
)



