
library(Seurat)
library(pheatmap)
library(tibble)
library(dplyr)
library(purrr)
library(gridExtra)
library(gtable)
library(grid)
#library(ComplexHeatmap)

library(RColorBrewer)
library(viridis)

setwd('/Users/alexanderjambor/workspace/wang/manuscripts/chung_mcd8t/Fig3L')
sobj = readRDS('input/TP16_Arm_clean_simpler4clustersn_20231020.rds') # 7538 cells x 32285 genes (2000 var.)
group_guides = read.csv('query-guides.csv', header=TRUE)

guide_groups = group_guides['group']
guides = group_guides['gRNA.ID']
names(guides) = 'V1'

signatures = read.csv('3h-signatures.heatmap2.csv', header=FALSE)
rownames(signatures) = signatures$V1
signatures = signatures[-1]
signatures = signatures[, colSums(is.na(signatures)) == 0]

query_genes = unique(unlist(signatures))
query_genes = query_genes[query_genes != ""] # 73 genes

sobj = sobj[,sobj@meta.data$guide_ID %in% guides$V1] # 16 of 19 guide IDs found
#sobj = sobj[,sobj$sample == 'Arm_SPL']
#sobj = sobj[,sobj$sample == 'Arm_IEL']
sobj = sobj[query_genes, ] # 3407 cells, 107 genes (60 var.)

rna_data = LayerData(sobj, assay='RNA', layer='data') # dgCMatrix
rna_df = as.data.frame(t(as.data.frame(rna_data)))
rna_df$cell = rownames(rna_df)
rownames(rna_df) = NULL

guides_df = as.data.frame(sobj@meta.data[, c("guide_ID")])
colnames(guides_df) = 'guide_ID'
guides_df$cell = Cells(sobj)

rna_df <- left_join(rna_df, guides_df, by = "cell")

avg_exp = rna_df %>%
  group_by(guide_ID) %>%
  summarise(across(-c(cell), \(x) mean(x, na.rm=TRUE)))


avg_exp = as.data.frame(avg_exp)
rownames(avg_exp) = avg_exp$guide_ID
avg_exp = subset(avg_exp, select=-c(guide_ID))

norm_exp = scale(avg_exp)


gene_col_order = c()
signature_groups = c()
for (r in rownames(signatures)) {
  tmp = as.character(signatures[r,])
  tmp = tmp[tmp != ""]
  signature_groups = c(signature_groups, rep(r, times=length(tmp)))
  gene_col_order = c(gene_col_order, tmp)
}


norm_exp = norm_exp[guides$V1, gene_col_order]

# now visualize

gene_annots = data.frame(signature_group = as.factor(signature_groups))
gene_annots$signature_group = factor(gene_annots$signature_group, levels=unique(gene_annots$signature_group))

rownames(guide_groups) = guides$V1

ann_colors = list(
  signature_group = c(Trm_Up = "#548135", Trm_Down = "#bf1d28"),
  group = c(Control = "#7f7f7f", Trm_selective = "#0432ff", Multi = "#d59ea6", non_selective = "#000000", Tex_selective = "#8eaadb")
)

rownames(gene_annots) = colnames(norm_exp)
change_indices <- which(gene_annots$signature_group[-1] != gene_annots$signature_group[-length(gene_annots$signature_group)])

# threshold to [-2.5, +2.5]
#tmp = apply(norm_exp, 1:2, function(x) min(x, 2.5))
#norm_exp2 = as.data.frame(apply(tmp,1:2, function(x) max(x, -2.5)))
norm_exp2 = norm_exp

ctrl_guides = c('gScramble-4')

#non_cluster_subset = norm_exp2[ctrl_guides, ]
#rows_to_exclude = rownames(norm_exp2) %in% ctrl_guides
#cluster_subset = norm_exp2[!rows_to_exclude, ]

#clustered_rows = hclust(dist(cluster_subset))
#ordered_cluster_subset = cluster_subset[clustered_rows$order, ]
#final_norm_exp = rbind(non_cluster_subset, ordered_cluster_subset)

final_norm_exp = norm_exp2



# original attempt
pheatmap_out = pheatmap(final_norm_exp, annotation_col = gene_annots, annotation_colors = ann_colors, gaps_row = c(1),
                        cluster_cols = FALSE, annotation_names_col = FALSE, show_colnames = TRUE, cellheight = 10, cellwidth = 10,
                        border_color = NA, treeheight_row = 0, gaps_col = change_indices, fontsize_col = 8, angle_col = '315', cluster_rows = FALSE, 
                        annotation_row = guide_groups,
                        #color = viridis(100, option = "D"))
                        #color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100))
                        )

###### return to pheatmap

# library(gridExtra)
# 
# non_cluster_subset = norm_exp2[ctrl_guides, ]
# rows_to_exclude = rownames(norm_exp2) %in% ctrl_guides
# cluster_subset = norm_exp2[!rows_to_exclude, ]
# 
# overall_min <- min(min(cluster_subset), min(non_cluster_subset))
# overall_max <- max(max(cluster_subset), max(non_cluster_subset))
# 
# color_palette <- colorRampPalette(c("blue", "white", "red"))(100)
# 
# heatmap_1 = pheatmap(non_cluster_subset, annotation_col = gene_annots, annotation_colors = ann_colors, gaps_row = c(2),
#                      cluster_cols = FALSE, annotation_names_col = FALSE, show_colnames = FALSE, cellheight = 10, cellwidth = 10,
#                      border_color = NA, gaps_col = change_indices, fontsize_col = 0, cluster_rows = FALSE, silent = TRUE, 
#                      annotation_legend = TRUE, legend = FALSE, margins = c(0,5,5,5),
#                      breaks = seq(overall_min, overall_max, length.out = 101),
#                      #color = viridis(100, option = "D"))
#                      color = colorRampPalette(brewer.pal(11, "RdBu"))(100))
# 
# 
# heatmap_2 = pheatmap(cluster_subset, gaps_col = change_indices, margins = c(5,5,0,5),
#                      cluster_cols = FALSE, annotation_names_col = FALSE, show_colnames = TRUE, cellheight = 10, cellwidth = 10,
#                      border_color = NA, treeheight_row = 25, fontsize_col = 8, angle_col = '315', cluster_rows = TRUE, silent = TRUE,
#                      breaks = seq(overall_min, overall_max, length.out = 101),
#                      #color = viridis(100, option = "D"))
#                      color = colorRampPalette(brewer.pal(11, "RdBu"))(100))
# 
# 
# 
# library(grid)
# 
# grob1 <- heatmap_1$gtable
# grob2 <- heatmap_2$gtable
# 
# max_width <- grid::unit.pmax(grob1$widths, grob2$widths)
# grob1$widths <- max_width
# grob2$widths <- max_width
# 
# combined_grob <- gtable_rbind(grob1, grob2, size = "first")
# grid.newpage()
# grid::grid.draw(combined_grob)

