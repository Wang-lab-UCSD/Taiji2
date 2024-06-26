library(ggplot2)
library(xlsx)
library(dplyr)
library(RColorBrewer)

maindir <- '~/allcombn/'
setwd(maindir)

### read data
df <- read.xlsx('gsea_intersection_armSPLIEL_heatmap_highlight genesets 20240212.xlsx', sheetIndex = 1)

### transform data
df2 <- df %>% dplyr::select(!starts_with('regulated')) %>% tidyr::pivot_longer(!path, names_to = c('category','gRNA'), names_sep = '_', values_to = 'value') %>% 
  mutate(TF = gsub('\\..*','',gsub('g','',gRNA))) %>% tidyr::pivot_wider(names_from = 'category', values_from = 'value') %>% 
  mutate(size = as.factor(ifelse(adjP > 0.05, 1, ifelse(adjP>0.005, 2, 3))))

## plot
p <- ggplot(df2, aes(TF, path)) + geom_point(aes(color = NES, size = size)) + 
  scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100)) + 
  theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
  scale_size_discrete(name = "p.adjust", labels = c("ns", "<0.05", "<0.005"))
pdf('fig3I.pdf', width = 14, height = 93)
print(p)
dev.off()

### selected pathways and TFs
TFs <- gsub('g','',readLines('fig3I_tfs.txt'))
pathways <- readLines('fig3I_pathways.txt')

df3 <- df2 %>% filter(TF %in% TFs & path %in% pathways)

# custom order
p <- ggplot(df3, aes(x=factor(TF, levels = TFs), y=factor(path, levels=pathways))) + geom_point(aes(color = NES, size = size)) + 
  scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100)) + 
  theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.title = element_blank())+
  scale_size_discrete(name = "p.adjust", labels = c("ns", "<0.05", "<0.005"))
pdf('fig3I_v2.pdf', width = 9.5, height = 6.5)
print(p)
dev.off()


# hclust order by NES
df4 <- df3 %>% select(path, TF, NES) %>% tidyr::pivot_wider(names_from = 'TF', values_from = 'NES') %>% tibble::column_to_rownames('path')
p <- pheatmap::pheatmap(df4, cellwidth = 10, cellheight = 10, filename = 'fig3I_hp_hclust_by_NES.pdf')

pathway_order <- rownames(df4)[p$tree_row$order]
tf_order <- names(df4)[p$tree_col$order]
p <- ggplot(df3, aes(x=factor(TF, levels = tf_order), y=factor(path, levels = pathway_order))) + geom_point(aes(color = NES, size = size)) + 
  scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100)) + 
  theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.title = element_blank())+
  scale_size_discrete(name = "p.adjust", labels = c("ns", "<0.05", "<0.005"))
pdf('fig3I_order_by_NES.pdf', width = 9.5, height = 6.5)
print(p)
dev.off()

# hclust order by p-value
df4 <- df3 %>% select(path, TF, adjP) %>% tidyr::pivot_wider(names_from = 'TF', values_from = 'adjP') %>% tibble::column_to_rownames('path')
p <- pheatmap::pheatmap(df4, cellwidth = 10, cellheight = 10, filename = 'fig3I_hp_hclust_by_adjP.pdf')

pathway_order <- rownames(df4)[p$tree_row$order]
tf_order <- names(df4)[p$tree_col$order]
p <- ggplot(df3, aes(x=factor(TF, levels = tf_order), y=factor(path, levels = pathway_order))) + geom_point(aes(color = NES, size = size)) + 
  scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100)) + 
  theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.title = element_blank())+
  scale_size_discrete(name = "p.adjust", labels = c("ns", "<0.05", "<0.005"))
pdf('fig3I_order_by_pvalue.pdf', width = 9.5, height = 6.5)
print(p)
dev.off()
