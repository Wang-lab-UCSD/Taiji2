
library(Seurat)
library(pheatmap)
library(tibble)
library(dplyr)
library(purrr)
library(gridExtra)
library(gtable)
library(grid)
library(ggplot2)
library(fgsea)
library(viridis)
library(RColorBrewer)

setwd('/Users/alexanderjambor/workspace/wang/manuscripts/chung_mcd8t/Fig4D')

# z_score = function(x) {
#   (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
# }

CreateNESMatrix = function(fname, gene_col_order) {
  # k=16 clusters
  tmp = read.csv('kmeans_cluster_k_16_log2FC_btw_TEX_and_TRM_mean_edge_weight_top500_subset_TFs_v2_v2_ordered_by_group.csv')
  cluster_row_order = c('8', '1', '9', '12', '7', '11', '16', '2', '10', '4', '13', '6', '14', '3', '5', '15')
  
  # k=20 clusters
  # tmp = read.csv('kmeans_cluster_k_20_log2FC_btw_TEX_and_TRM_mean_edge_weight_subset_TFs_v2_v2_ordered_by_group.csv')
  # cluster_row_order = c('19','20','5','12','6','7','15','10','18','9','17','8','1','13','2','11','3','16','4','14')

  cluster_ids = unique(tmp$kmeans_cluster)
  
  signature_groups = c('Trm', 
                       'Multitasker','Multitasker','Multitasker','Multitasker','Multitasker','Multitasker','Multitasker',
                       'Tex','Tex','Tex','Tex','Tex','Tex','Tex','Tex','Tex')
  
  if (grepl('Arm', fname)) {
    title = "Arm NES heatmap"
    col_gaps = c(1,8)
  } else {
    title = "Cl13 NES heatmap"
    col_gaps = c(7)
  }

  sobj = readRDS(fname) # 7538 cells x 32285 genes (2000 var.); 22 perturbed genes
  sobj = sobj[, sobj$perturbed_gene %in% c('gScramble', gene_col_order)] # 16 perturbed genes (missing: gFosb)

  custom_gene_sets = list()
  for (id in cluster_ids) {
    custom_gene_sets[[id]] = sapply(tmp$gene[tmp$kmeans_cluster == id], function(x) tools::toTitleCase(tolower(x)))
  }
  names(custom_gene_sets) = cluster_ids

  all_genes = unlist(custom_gene_sets, use.names=FALSE) # 21432 genes
  unique_genes = unique(all_genes) # 21432 genes, so can just randomly select gene sets

  # Ok, just use Seurat FindMarkers, get log2FC, ranked gene list, etc.
  de_results = list()
  background_cells = WhichCells(sobj, expression = perturbed_gene == 'gScramble')

  for (gene in gene_col_order) {
    query_cells = WhichCells(sobj, expression = perturbed_gene == gene)
    res = FindMarkers(sobj, ident.1 = query_cells, ident.2 = background_cells, min.pct = 0.1)
    de_results[[gene]]= res[order(res$avg_log2FC, decreasing = TRUE), ]
  }

  gsea_results = list()

  # For now, just rank by log2FC but may want
  for (gene in names(de_results)) {
    gene_ranks = de_results[[gene]]$avg_log2FC
    names(gene_ranks) = rownames(de_results[[gene]])
    gsea_results[[gene]] = fgsea(pathways=custom_gene_sets, stats=gene_ranks, minSize=0, maxSize=5000)
  }

  for (x in gsea_results) {
    print(length(x$NES))
  }

  ##### Violin plot of cluster 37 NES across all gRNA KOs

  # nes_data <- data.frame(nes = numeric(), perturbed_gene = factor())
  #
  # for (gene in names(gsea_results)) {
  #   cluster_nes = gsea_results[[gene]][gsea_results[[gene]]$pathway == 37, "NES"]
  #   temp_df = data.frame(nes = cluster_nes, perturbed_gene = gene)
  #   nes_data = rbind(nes_data, temp_df)
  # }
  #
  # nes_data$perturbed_gene = factor(nes_data$perturbed_gene)
  #
  # Pre-compute the jitter values
  # nes_data$jitter_x <- 1  # All points initially at 1 on the x-axis for alignment with the violin
  # nes_data$jitter_x <- jitter(nes_data$jitter_x, amount = 0.3)  # Add jitter

  # Create the ggplot
  # ggplot(nes_data, aes(x = jitter_x, y = NES)) +
  #   geom_violin(data = data.frame(x = 1, y = nes_data$NES), aes(x = x, y = y), trim = FALSE, width = 0.75) +
  #   geom_point(aes(color = perturbed_gene), size = 2, alpha = 0.6) +
  #   geom_text(aes(label = perturbed_gene), vjust = -0.5, hjust = 0, check_overlap = FALSE) +  # Adjust text alignment
  #   theme_bw() +  # Use a theme with a white background
  #   theme(axis.text.x = element_blank(),  # Hide x-axis text
  #         axis.ticks.x = element_blank(),  # Hide x-axis ticks
  #         axis.title.x = element_blank(),  # Hide x-axis title
  #         legend.position = "none") +  # Hide the legend
  #   labs(y = "Normalized Enrichment Score", title = "NES Across gRNAs for Cluster 37 (Arm)")

  ##### Create heatmap

  nes_matrix_list = list()

  for (gRNA in names(gsea_results)) {
    gene_sets = gsea_results[[gRNA]]$pathway
    nes_values = gsea_results[[gRNA]]$NES
    nes_vector = setNames(nes_values, gene_sets)
    nes_matrix_list[[gRNA]] = nes_vector
  }

  nes_matrix = do.call(rbind, nes_matrix_list)
  nes_matrix = t(nes_matrix)
  nes_matrix = nes_matrix[cluster_row_order, gene_col_order]

  ann_colors = list(
    signature_group = c(Trm = "#139C3B", Multitasker = "#7F139C", Tex = "#548135"))

  gene_annots = data.frame(signature_group = as.factor(signature_groups))
  gene_annots$signature_group = factor(gene_annots$signature_group, levels=unique(gene_annots$signature_group))
  
  annot_gRNAs = c('gFosb', 'gNr4a2','gStat3', 'gNfil3', 'gHic1', 'gPrdm1',  'gGfi1', 'gIkzf3','gTfdp1',
                  'gIrf8', 'gArid3a',  'gEtv5', 'gJdp2', 'gZfp410', 'gZscan20', 'gZfp324', 'gNfatc1')
  
  rownames(gene_annots) = annot_gRNAs

  norm_nes_matrix = t(scale(t(nes_matrix)))
  #norm_nes_matrix = nes_matrix
  
  val_min = min(norm_nes_matrix)
  val_max = max(norm_nes_matrix)

  heatmap = pheatmap(norm_nes_matrix, margins = c(5,5,0,5), annotation_col = gene_annots, annotation_colors = ann_colors, 
                     gaps_row = c(4,7),
                     cluster_cols = FALSE, annotation_names_col = FALSE, show_colnames = TRUE, cellheight = 10, cellwidth = 14,
                     border_color = NA, treeheight_row = 0, treeheight_col = 0, fontsize_col = 8, angle_col = '315', 
                     cluster_rows = FALSE,
                     #cluster_rows = TRUE,
                     silent = TRUE, breaks = seq(val_min, val_max, length.out = 101), main = title, 
                     gaps_col = col_gaps,
                     #gaps_col = c(1, col_gaps+1), # w/gScramble
                     #color = viridis(100, option = "D"))
                     #color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100))
                     )

  return(list(heatmap = heatmap, data = nes_matrix))
}
    





CreateModuleScoreMatrix = function(fname, offset, gene_col_order) {
  #browser()
  
  # k=16 clusters
  tmp = read.csv('kmeans_cluster_k_16_log2FC_btw_TEX_and_TRM_mean_edge_weight_top500_subset_TFs_v2_v2_ordered_by_group.csv')
  cluster_row_order = c('8', '1', '9', '12', '7', '11', '16', '2', '10', '4', '13', '6', '14', '3', '5', '15')
  
  # k=20 clusters
  # tmp = read.csv('kmeans_cluster_k_20_log2FC_btw_TEX_and_TRM_mean_edge_weight_subset_TFs_v2_v2_ordered_by_group.csv')
  # cluster_row_order = c('19','20','5','12','6','7','15','10','18','9','17','8','1','13','2','11','3','16','4','14')
  
  cluster_ids = unique(tmp$kmeans_cluster)
  
  signature_groups = c('Trm', 
                       'Multitasker','Multitasker','Multitasker','Multitasker','Multitasker','Multitasker','Multitasker',
                       'Tex','Tex','Tex','Tex','Tex','Tex','Tex','Tex','Tex')
  
  if (grepl('Arm', fname)) {
    title = "Arm Module Scores"
    col_gaps = c(1,8)
  } else {
    title = "Cl13 Module Scores"
    col_gaps = c(7)
  }
  
  sobj = readRDS(fname) # 7538 cells x 32285 genes (2000 var.); 22 perturbed genes
  #sobj = subset(sobj, subset = sample == 'Arm_SPL')  # OPTIONAL
  
  #perturbed_genes = read.csv('perturbed-genes.csv', header=FALSE) # 17 perturbed genes
  sobj = sobj[, sobj$perturbed_gene %in% gene_col_order]

  
  custom_gene_sets = list()
  for (id in cluster_ids) {
    custom_gene_sets[[id]] = sapply(tmp$gene[tmp$kmeans_cluster == id], function(x) tools::toTitleCase(tolower(x)))
  }
  names(custom_gene_sets) = cluster_ids
  
  all_genes = unlist(custom_gene_sets, use.names=FALSE) # 21432 genes
  unique_genes = unique(all_genes) # 21432 genes, so can just randomly select gene sets
  
  sobj = AddModuleScore(object = sobj, 
                        features = custom_gene_sets,
                        name = "ModuleScore")
  
  library(dplyr)
  library(tidyr)
  
  # Convert Seurat object to a data frame for easier manipulation
  df_scores = sobj@meta.data %>% 
    dplyr::select(starts_with("ModuleScore"), perturbed_gene)
  
  
  # Gather scores into long format for grouping
  df_long = pivot_longer(df_scores, cols = starts_with("ModuleScore"), names_to = "GeneSet", values_to = "Score")
  
  # Calculate average score for each gene set and perturbed gene
  average_scores = df_long %>%
    group_by(GeneSet, perturbed_gene) %>%
    summarize(AvgScore = mean(Score, na.rm = TRUE)) %>%
    ungroup()
  
  #browser()
  
  # Spread the average scores back into a wide format (gene set id x perturbed gene) matrix
  heatmap_matrix = pivot_wider(average_scores, names_from = perturbed_gene, values_from = AvgScore)
  
  # Convert to matrix if needed
  heatmap_matrix = as.matrix(heatmap_matrix[,-1]) # Assuming the first column is GeneSet names
  rownames(heatmap_matrix) = average_scores$GeneSet[!duplicated(average_scores$GeneSet)]
  rownames(heatmap_matrix) = c(unlist(lapply(rownames(heatmap_matrix), function (x) gsub("ModuleScore", "", x))))
  heatmap_matrix = heatmap_matrix[cluster_row_order, gene_col_order]
  
  if (offset == 'none') {
    norm_matrix = heatmap_matrix
  } else if (offset == 'gScramble') {
    norm_matrix = heatmap_matrix - heatmap_matrix[, 'gScramble', drop = TRUE]
    norm_matrix = norm_matrix[, !(colnames(norm_matrix) %in% 'gScramble')]
  }
  
  scaled_norm_matrix = t(scale(t(norm_matrix))) # scale across gRNAs
  
  val_min = min(scaled_norm_matrix)
  val_max = max(scaled_norm_matrix)
  
  ann_colors = list(
    signature_group = c(Trm = "#139C3B", Multitasker = "#7F139C", Tex = "#548135"))
  
  gene_annots = data.frame(signature_group = as.factor(signature_groups))
  gene_annots$signature_group = factor(gene_annots$signature_group, levels=unique(gene_annots$signature_group))
  
  annot_gRNAs = c('gFosb', 'gNr4a2','gStat3', 'gNfil3', 'gHic1', 'gPrdm1',  'gGfi1', 'gIkzf3','gTfdp1',
                  'gIrf8', 'gArid3a',  'gEtv5', 'gJdp2', 'gZfp410', 'gZscan20', 'gZfp324', 'gNfatc1')
  
  rownames(gene_annots) = annot_gRNAs
  
  heatmap = pheatmap(scaled_norm_matrix, margins = c(5,5,0,5), annotation_col = gene_annots, annotation_colors = ann_colors, 
                     #gaps_row = c(7,14),
                     gaps_row = c(4, 7),
                     cluster_cols = FALSE, annotation_names_col = FALSE, show_colnames = TRUE, cellheight = 10, cellwidth = 14,
                     border_color = NA, treeheight_row = 0, treeheight_col = 0, fontsize_col = 8, angle_col = '315',
                     cluster_rows = FALSE, 
                     #cluster_rows = TRUE, 
                     silent = TRUE, breaks = seq(val_min, val_max, length.out = 101), main = title, 
                     gaps_col = col_gaps,
                     #gaps_col = c(1, col_gaps+1), # w/gScramble
                    #color = viridis(100, option = "D"))
                    color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100))
                    #)

                    
  return(list(heatmap = heatmap, data = norm_matrix))
}

### NES calculations
# 
# # Arm data
# gene_col_order = c('gFosb', 'gNr4a2','gStat3', 'gNfil3', 'gHic1', 'gPrdm1',  'gGfi1', 'gIkzf3','gTfdp1',
#                    'gIrf8', 'gArid3a',  'gEtv5', 'gJdp2', 'gZfp410', 'gZscan20', 'gZfp324', 'gNfatc1')
# 
# res = CreateNESMatrix('input/TP16_Arm_clean_simpler4clustersn_20231020.rds', gene_col_order)
# Arm_nes_heatmap = res$heatmap
# Arm_nes_matrix = res$data
#  
# # Cl13 data
# gene_col_order = c('gNr4a2','gStat3', 'gNfil3', 'gHic1', 'gPrdm1',  'gGfi1', 'gIkzf3','gTfdp1',
#                    'gIrf8', 'gArid3a',  'gEtv5', 'gJdp2', 'gZfp410', 'gZscan20', 'gZfp324', 'gNfatc1')
# 
# res = CreateNESMatrix('input/Cl13_23TP04_20231010.rds', gene_col_order)
# Cl13_nes_heatmap = res$heatmap
# Cl13_nes_matrix = res$data
# 
# ###
# 
# gene_col_order = c('gNr4a2','gStat3', 'gNfil3', 'gHic1', 'gPrdm1',  'gGfi1', 'gIkzf3','gTfdp1',
#                    'gIrf8', 'gArid3a',  'gEtv5', 'gJdp2', 'gZfp410', 'gZscan20', 'gZfp324', 'gNfatc1')
# 
# signature_groups = c( 'Multitasker','Multitasker','Multitasker','Multitasker','Multitasker','Multitasker','Multitasker',
#                       'Tex','Tex','Tex','Tex','Tex','Tex','Tex','Tex','Tex')
# 
# Arm_nes_matrix = Arm_nes_matrix[, gene_col_order]
# Cl13_nes_matrix = Cl13_nes_matrix[, gene_col_order]
# 
# gene_annots = data.frame(signature_group = as.factor(signature_groups))
# gene_annots$signature_group = factor(gene_annots$signature_group, levels=unique(gene_annots$signature_group))
# rownames(gene_annots) = gene_col_order
# 
# ann_colors = list(
#   signature_group = c(Multitasker = "#7F139C", Tex = "#548135"))
# 
# 
# diff_mtx = abs(Cl13_nes_matrix - Arm_nes_matrix)
# diff_mtx = t(scale(t(diff_mtx)))
# 
# 
# heatmap = pheatmap(diff_mtx, margins = c(5,5,0,5), annotation_col = gene_annots, annotation_colors = ann_colors, gaps_row = c(4,7),
#                    cluster_cols = FALSE, annotation_names_col = FALSE, show_colnames = TRUE, cellheight = 10, cellwidth = 14,
#                    border_color = NA, treeheight_row = 0, treeheight_col = 0, fontsize_col = 8, angle_col = '315', cluster_rows = FALSE, 
#                    silent = TRUE, breaks = seq(min(diff_mtx), max(diff_mtx), length.out = 101), main = "Abs(NES Difference) Cl13/Arm heatmap", gaps_col = c(7),
#                    #color = viridis(100, option = "D"))
#                    #color = colorRampPalette(brewer.pal(11, "BrBG"))(100))
#                    color = colorRampPalette(brewer.pal(11, "YlOrBr"))(100))
# #)
# heatmap


# Arm data
gene_col_order = c('gFosb', 'gNr4a2','gStat3', 'gNfil3', 'gHic1', 'gPrdm1',  'gGfi1', 'gIkzf3','gTfdp1',
                   'gIrf8', 'gArid3a',  'gEtv5', 'gJdp2', 'gZfp410', 'gZscan20', 'gZfp324', 'gNfatc1')

res = CreateModuleScoreMatrix('input/TP16_Arm_clean_simpler4clustersn_20231020.rds', 'none', gene_col_order)
Arm_mod_heatmap = res$heatmap
Arm_mod_matrix = res$data

# Cl13 data
gene_col_order = c('gNr4a2','gStat3', 'gNfil3', 'gHic1', 'gPrdm1',  'gGfi1', 'gIkzf3','gTfdp1',
                   'gIrf8', 'gArid3a',  'gEtv5', 'gJdp2', 'gZfp410', 'gZscan20', 'gZfp324', 'gNfatc1')

res = CreateModuleScoreMatrix('input/Cl13_23TP04_20231010.rds', 'none', gene_col_order)
Cl13_mod_heatmap = res$heatmap
Cl13_mod_matrix = res$data


# module scores applied to ea. cell which is why log2FC isn't used w.r.t. control

gene_col_order = c('gNr4a2','gStat3', 'gNfil3', 'gHic1', 'gPrdm1',  'gGfi1', 'gIkzf3','gTfdp1',
                   'gIrf8', 'gArid3a',  'gEtv5', 'gJdp2', 'gZfp410', 'gZscan20', 'gZfp324', 'gNfatc1')

signature_groups = c( 'Multitasker','Multitasker','Multitasker','Multitasker','Multitasker','Multitasker','Multitasker',
                      'Tex','Tex','Tex','Tex','Tex','Tex','Tex','Tex','Tex')

Arm_mod_matrix = Arm_mod_matrix[, gene_col_order]
Cl13_mod_matrix = Cl13_mod_matrix[, gene_col_order]


gene_annots = data.frame(signature_group = as.factor(signature_groups))
gene_annots$signature_group = factor(gene_annots$signature_group, levels=unique(gene_annots$signature_group))
rownames(gene_annots) = gene_col_order

ann_colors = list(
  signature_group = c(Multitasker = "#7F139C", Tex = "#548135"))


diff_mtx = abs(Cl13_mod_matrix - Arm_mod_matrix)
diff_mtx = t(scale(t(diff_mtx)))


heatmap = pheatmap(diff_mtx, margins = c(5,5,0,5), annotation_col = gene_annots, annotation_colors = ann_colors, gaps_row = c(4,7),
                   cluster_cols = FALSE, annotation_names_col = FALSE, show_colnames = TRUE, cellheight = 10, cellwidth = 14,
                   border_color = NA, treeheight_row = 0, treeheight_col = 0, fontsize_col = 8, angle_col = '315', cluster_rows = FALSE, 
                   silent = TRUE, breaks = seq(min(diff_mtx), max(diff_mtx), length.out = 101), main = "Abs(Module Score Difference) Cl13/Arm heatmap", gaps_col = c(7),
                   #color = viridis(100, option = "D"))
                   #color = colorRampPalette(brewer.pal(11, "BrBG"))(100))
                   color = colorRampPalette(brewer.pal(11, "YlOrBr"))(100))
                   #)
heatmap



