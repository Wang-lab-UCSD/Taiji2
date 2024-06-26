library(ggplot2)
library(xlsx)
library(dplyr)
library(RColorBrewer)

maindir <- '~/allcombn/'
setwd(maindir)

### Taiji prediction 
### read data
L <- list.files(path = '21_targets_csv/',pattern = '.*csv')
df <- do.call('rbind',lapply(L, function(x) read.csv(paste0('21_targets_csv/',x)) %>% mutate(TF=gsub('\\.csv','',x))))

### transform data
df2 <- df %>% mutate(size = as.factor(ifelse(pvalue > 0.05, 1, ifelse(pvalue>0.005, 2, 3)))) %>% filter(TF %in% TFs & Description %in% pathways)

## plot
p <- ggplot(df2, aes(x=factor(TF, levels = TFs), y=factor(Description, levels=pathways))) + geom_point(aes(color = NES, size = size)) + 
  scale_color_stepsn(colours = c("#e1a692", "#de6e56", "#e14b31"), breaks = c(0.5,1,1.5,2))+
  # scale_color_distiller(palette = "RdYlBu")+
  theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.title = element_blank())+
  scale_size_discrete(name = "p.adjust", labels = c("ns", "<0.05", "<0.005"))
# pdf('fig5A_rdylbu.pdf', width = 9.5, height = 3)
pdf('fig5A_v2.pdf', width = 9.5, height = 3)
print(p)
dev.off()





### perturb-seq data
### read data
df <- read.csv('gseaimmune_intersection_cl13_noNA_23TP04.csv')


### transform data
df2 <- df %>% dplyr::select(!starts_with('regulated')) %>% tidyr::pivot_longer(!path, names_to = c('category','gRNA'), names_sep = '_', values_to = 'value') %>% 
  mutate(TF = gsub('\\..*','',gsub('g','',gRNA))) %>% tidyr::pivot_wider(names_from = 'category', values_from = 'value') %>% 
  mutate(size = as.factor(ifelse(adjP > 0.05, 1, ifelse(adjP>0.005, 2, 3))))

### selected pathways and TFs
gTFs <- gsub('-','.',readLines('fig5b_tfs.txt'))
TFs <- readLines('fig5_tfs.txt')
pathways <- rev(readLines('fig5_pathways.txt'))
df3 <- df2 %>% filter(gRNA %in% gTFs & path %in% pathways)

# custom order
p <- ggplot(df3, aes(x=factor(TF, levels = TFs), y=factor(path, levels=pathways))) + geom_point(aes(color = NES, size = size)) + 
  # scale_color_stepsn(colours = c("#1984c5", "#22a7f0", "#63bff0", "#a7d5ed", "#e2e2e2","#e2e2e2", "#e1a692", "#de6e56", "#e14b31", "#c23728"), breaks = c(-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5))+
  scale_color_distiller(palette = "RdYlBu")+
  theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.title = element_blank())+
  scale_size_discrete(name = "p.adjust", labels = c("ns", "<0.05", "<0.005"))
pdf('fig5B_rdylbu.pdf', width = 9.5, height = 2.5)
print(p)
dev.off()
