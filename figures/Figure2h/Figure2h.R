library(Seurat)
library(pheatmap)
library(tibble)
library(dplyr)
library(purrr)
library(gridExtra)
library(gtable)
library(grid)
library(viridis)
library(RColorBrewer)
library(fgsea)
library(reshape2)
library(ggplot2)
library(tidyr)
library(grid)
library(this.path)

##### MAIN #####

setwd(this.path::here())

GO_term_to_description = list(GO_term = c("GO:0031331",	"GO:0031647",	"GO:1903052",	"GO:0043161",	"GO:0061919",	"GO:0097193",	"GO:0032623",	"GO:0032743",	"GO:0032663",	"GO:0071559",	"GO:0007015",	"GO:2000112",	"GO:0043087",	"GO:0007159"),
                              Description = c('Catabolism', 'Protein stability', 'Proteolysis', 'UPS', 'Autophagy', 'Intrinsic apoptosis', 'IL-2 production', 'Pos. reg. of IL-2 production', 'Reg. of IL-2 production', 'TGF-Beta', 'Actin', 'Reg. of cellular macromolecule biosynthetic process', 'Reg. of GTPase activity', 'Leukocyte cell-cell adhesion'))

signatures = read.table('SubGenesets_TRM_TEX_TFcommunity_20240419.csv', header=FALSE, sep=',')

rownames(signatures) = signatures$V1
signatures = signatures[-1]
custom_gene_sets = list()

for (r in rownames(signatures)) {
  tmp = as.character(signatures[r,])
  tmp = tmp[tmp != ""]
  custom_gene_sets[[r]] = toupper(tmp)
}

custom_gene_sets = lapply(custom_gene_sets, function (x) x[x != "NA"])
gene_set_order = names(custom_gene_sets)

### Part 1: PanCancer scRNA-seq GSEA
sobj = readRDS('Pancancer_CD8_simple_7cluster_220518.rds')

de_results = list()
background_cells = WhichCells(sobj, idents = 'TRM')
query_cells = WhichCells(sobj, idents = 'Tex')

de_res = FindMarkers(sobj, ident.1 = query_cells, ident.2 = background_cells, min.pct = 0.1, assay = 'RNA', slot = 'data')
de_res = de_res[order(de_res$avg_log2FC, decreasing = TRUE), ]

gene_ranks = setNames(de_res$avg_log2FC, rownames(de_res))
gsea_res = fgsea(pathways=custom_gene_sets, stats=gene_ranks, minSize=0, maxSize=5000)
# note: GSEA typically does not benefit from filtering

nes_vec = setNames(gsea_res$NES, gsea_res$pathway)
padj_vec = setNames(gsea_res$padj, gsea_res$pathway)
tmp_df = as.data.frame(cbind(padj_vec, nes_vec))
colnames(tmp_df) = c('padj', 'NES')
tmp_df$Condition = "PanCancer"

tmp_df$padj_size = (-log10(tmp_df$padj))^0.5
tmp_df = tmp_df[gene_set_order, ]
tmp_df$GeneSet = rownames(tmp_df)
tmp_df$GeneSet = sapply(tmp_df$GeneSet, function(x) paste0(x, ' (', GO_term_to_description$Description[GO_term_to_description$GO_term == x], ')'))
tmp_df$GeneSet = factor(tmp_df$GeneSet, levels = unique(tmp_df$GeneSet))
rownames(tmp_df) = NULL

tmp_df1 = tmp_df

### Part 2: LCMV scRNA-seq GSEA
library(edgeR)

cts_mtx = read.csv('Taiji_counts_matrix.tsv', sep='\t', row.names=1)
exclude_idx = grep('*lung*', colnames(cts_mtx))
cts_mtx = cts_mtx[, -exclude_idx]

cell_type = factor(unname(sapply(colnames(cts_mtx), function(x) strsplit(x, '_')[[1]][[1]])))
experiment = factor(unname(sapply(colnames(cts_mtx), function(x) strsplit(x, '_')[[1]][[2]])))

dge = DGEList(counts = cts_mtx)
dge$samples$group = relevel(cell_type, 'TRM')
dge$samples$experiment = experiment

des = model.matrix(~1 + group, data=dge$samples)
# des = model.matrix(~0 + group, data=dge$samples)
# contrast = makeContrasts(groupTexTerm - groupTRM, levels = des)

dge = calcNormFactors(dge)
print(dim(dge))

keep = filterByExpr(dge, des)
dge = dge[keep, , keep.lib.sizes = FALSE]
print(dim(dge))

dge = estimateDisp(dge, des)

fit = glmQLFit(dge, des)
glm_res = glmQLFTest(fit)
# glm_res = glmQLFTest(fit, contrast=contrast)

adj_pvals = p.adjust(glm_res$table$PValue, method='BY')

res_table = as.data.frame(topTags(glm_res, n=Inf)) # n=Inf to include all genes
res_table = res_table[order(res_table$logFC, decreasing=TRUE), ]

gene_ranks = setNames(res_table$logFC, toupper(rownames(res_table)))

gsea_res = fgsea(pathways=custom_gene_sets, stats=gene_ranks, minSize=0, maxSize=5000)

nes_vec = setNames(gsea_res$NES, gsea_res$pathway)
padj_vec = setNames(gsea_res$padj, gsea_res$pathway)
tmp_df = as.data.frame(cbind(padj_vec, nes_vec))
colnames(tmp_df) = c('padj', 'NES')
tmp_df$Condition = "LCMV"

tmp_df$padj_size = (-log10(tmp_df$padj))^0.5
tmp_df = tmp_df[gene_set_order, ]
tmp_df$GeneSet = rownames(tmp_df)
tmp_df$GeneSet = sapply(tmp_df$GeneSet, function(x) paste0(x, ' (', GO_term_to_description$Description[GO_term_to_description$GO_term == x], ')'))
tmp_df$GeneSet = factor(tmp_df$GeneSet, levels = unique(tmp_df$GeneSet))
rownames(tmp_df) = NULL

tmp_df2 = tmp_df

####### Composite dot plot

combined_df = rbind(tmp_df1, tmp_df2)

legend_padj_values = c(5e-2, 1e-1, 5e-1)
legend_sizes = (-log10(legend_padj_values))^0.5
legend_labels = c('5e-2', '1e-1', '5e-1')

extra_points = data.frame(GeneSet = rep(tail(combined_df$GeneSet, 1), length(legend_sizes)),
                          Condition = rep(tail(combined_df$Condition, 1), length(legend_sizes)),
                          NES = rep(min(combined_df$NES), length(legend_sizes)),  # Use a consistent value for NES
                          padj_size = legend_sizes)

fig = ggplot(combined_df, aes(y = GeneSet, x = Condition, color = NES, size = padj_size)) +
  geom_point() +
  geom_point(data = extra_points, shape=21, alpha=0) + 
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  scale_size_continuous(breaks = legend_sizes, labels = legend_labels) + 
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(color = "NES",
       size = "Adjusted p-value",
       title = "GO terms")

fig
