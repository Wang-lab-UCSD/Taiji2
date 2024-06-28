
library(Seurat)
library(pheatmap)
library(tibble)
library(dplyr)
library(purrr)
library(gridExtra)
library(gtable)
library(grid)

setwd('/Users/alexanderjambor/workspace/wang/manuscripts/chung_mcd8t/Fig2I')

sobj = readRDS('input/Cl13_23TP04_20231010.rds') # 7538 cells x 32285 genes (2000 var.)

guides = read.csv('query-guides.csv', header=FALSE)
signatures = read.csv('2f-signatures_aggr.csv', header=FALSE)
rownames(signatures) = signatures$V1
signatures = signatures[-1]

query_genes = unique(unlist(signatures, use.names = FALSE))
query_genes = query_genes[query_genes != ""] # 107 genes

sobj = sobj[,sobj@meta.data$guide_ID %in% guides$V1] 
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

z_score = function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

norm_tibble = avg_exp %>%
  mutate(across(where(is.numeric), z_score))

norm_exp = as.data.frame(norm_tibble)

rownames(norm_exp) = norm_exp$guide_ID
norm_exp = subset(norm_exp, select=-c(guide_ID))

# why is Ifnb1 all zeros? There's some sort of weird discrepancy, return to later
# norm_exp = subset(norm_exp, select=-c(Ifnb1))

gene_col_order = c()
signature_groups = c()
for (r in rownames(signatures)) {
  tmp = as.character(signatures[r,])
  tmp = tmp[tmp != ""]
  signature_groups = c(signature_groups, rep(r, times=length(tmp)))
  gene_col_order = c(gene_col_order, tmp)
}

norm_exp = norm_exp[, gene_col_order]

# now visualize





gene_annots = data.frame(signature_group = as.factor(signature_groups))
gene_annots$signature_group = factor(gene_annots$signature_group, levels=unique(gene_annots$signature_group))

ann_colors = list(
  signature_group = c(TexTerm = "#548135",
                      TexEffLike = "#8eaadb",
                      Tem = "#00b0f0", 
                      Teff = "#0432ff",
                      Tnaive = "#7f7f7f",
                      Tmem = "#bf1d28",
                      TexProg = "#d59ea6",
                      Ca = "#000000", 
                      Tcr = "#929292",
                      Ctk = "#d0cece"))

rownames(gene_annots) = colnames(norm_exp)

tmp = apply(norm_exp, 1:2, function(x) min(x, 2.5))
norm_exp2 = as.data.frame(apply(tmp,1:2, function(x) max(x, -2.5)))

change_indices <- which(gene_annots$signature_group[-1] != gene_annots$signature_group[-length(gene_annots$signature_group)])

non_cluster_subset = norm_exp2[c('gScramble-2', 'gScramble-4'), ]

rows_to_exclude = rownames(norm_exp2) %in% c('gScramble-2', 'gScramble-4')
cluster_subset = norm_exp2[!rows_to_exclude, ]

clustered_rows = hclust(dist(cluster_subset))
ordered_cluster_subset = cluster_subset[clustered_rows$order, ]
final_norm_exp = rbind(non_cluster_subset, ordered_cluster_subset)

library(gridExtra)

#non_cluster_subset = norm_exp2[c('gScramble-2', 'gScramble-4'), ]
#rows_to_exclude = rownames(norm_exp2) %in% c('gScramble-2', 'gScramble-4')
#cluster_subset = norm_exp2[!rows_to_exclude, ]

#overall_min <- min(min(cluster_subset), min(non_cluster_subset))
#overall_max <- max(max(cluster_subset), max(non_cluster_subset))

library(viridis)
library(RColorBrewer)

guide_order = as.vector(guides)$V1
norm_exp2 = norm_exp2[guide_order,]


dev.off()
heatmap1 = pheatmap(norm_exp2, annotation_col = gene_annots, annotation_colors = ann_colors, gaps_row = c(2, 17),
                     cluster_cols = FALSE, annotation_names_col = FALSE, show_colnames = TRUE, cellheight = 10, cellwidth = 10,
                     border_color = NA, gaps_col = change_indices, fontsize_col = 8, cluster_rows = FALSE, silent = TRUE,
                     annotation_legend = TRUE, legend = FALSE, margins = c(0,5,5,5), angle_col = '315',
                     breaks = seq(min(norm_exp2), max(norm_exp2), length.out = 101),
                     color = viridis(100, option = "D"))
                     #color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100))
                     #)

heatmap1

# heatmap_1 = pheatmap(non_cluster_subset, annotation_col = gene_annots, annotation_colors = ann_colors, gaps_row = c(2),
#                      cluster_cols = FALSE, annotation_names_col = FALSE, show_colnames = FALSE, cellheight = 10, cellwidth = 10,
#                      border_color = NA, gaps_col = change_indices, fontsize_col = 0, cluster_rows = FALSE, silent = TRUE, 
#                      annotation_legend = TRUE, legend = FALSE, margins = c(0,5,5,5),
#                      breaks = seq(overall_min, overall_max, length.out = 101),
#                      #color = viridis(100, option = "D"))
#                      color = colorRampPalette(brewer.pal(11, "RdBu"))(100))


# heatmap_2 = pheatmap(cluster_subset, gaps_col = change_indices, margins = c(5,5,0,5),
#                      cluster_cols = FALSE, annotation_names_col = FALSE, show_colnames = TRUE, cellheight = 10, cellwidth = 10,
#                      border_color = NA, treeheight_row = 25, fontsize_col = 8, angle_col = '315', cluster_rows = TRUE, silent = TRUE,
#                      breaks = seq(overall_min, overall_max, length.out = 101),
#                      #color = viridis(100, option = "D"))
#                      color = colorRampPalette(brewer.pal(11, "RdBu"))(100))

################
library(grid)

grob1 <- heatmap_1$gtable
grob2 <- heatmap_2$gtable

max_width <- grid::unit.pmax(grob1$widths, grob2$widths)
grob1$widths <- max_width
grob2$widths <- max_width

combined_grob <- gtable_rbind(grob1, grob2, size = "first")


# grid.arrange(grob1, grob2, ncol = 1, heights = c(0.5, 0.5)) # doesn't work, too much vertical space
grid.newpage()
grid::grid.draw(combined_grob)

