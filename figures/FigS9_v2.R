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

### selected pathways and TFs
TFs <- gsub('g','',readLines('fig3N_tfs.txt'))
pathways <- readLines('figS9_pathways.txt')
df3 <- df2 %>% filter(TF %in% TFs & path %in% pathways)

# custom order
p <- ggplot(df3, aes(x=factor(TF, levels = TFs), y=factor(path, levels=pathways))) + geom_point(aes(color = NES, size = size)) + 
  scale_color_stepsn(colours = c("#1984c5", "#22a7f0", "#63bff0", "#a7d5ed", "#e2e2e2", "#e1a692", "#de6e56", "#e14b31", "#c23728"), breaks = c(-2.5,-2,-1.5,-1,-0.5,0.5,1,1.5,2,2.5))+
  theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.title = element_blank())+
  scale_size_discrete(name = "p.adjust", labels = c("ns", "<0.05", "<0.005"))
pdf('figS9.pdf', width = 8.5, height = 4)
print(p)
dev.off()


