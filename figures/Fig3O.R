library(ggplot2)
library(xlsx)
library(dplyr)
library(RColorBrewer)

maindir <- '~/allcombn/'
setwd(maindir)

### read data
df <- read.xlsx('network/figS9_trm_regulatees/figS9_cut2000_ew_regulatees_gsea.xlsx', sheetIndex = 1)

### transform data
df2 <- df %>% mutate(size = as.factor(ifelse(pvalue > 0.05, 1, ifelse(pvalue>0.005, 2, 3))))

### selected pathways and TFs
TFs <- gsub('g','',readLines('fig3O_tfs.txt'))
pathways <- readLines('fig3N_pathways.txt')
df3 <- df2 %>% filter(TF %in% TFs & Description %in% pathways) %>% select(Description, TF, size, NES)

tmp1 <- data.frame(Description = rep(pathways, 5), TF = rep(TFs, each=10), size = 1, NES = 0)
df4 <- tmp1 %>% left_join(df3, by = c('Description','TF')) %>% mutate(size = factor(ifelse(is.na(size.y), 1, size.y)), NES = ifelse(is.na(NES.y), 0, NES.y))

# custom order
p <- ggplot(df4, aes(x=factor(TF, levels = TFs), y=factor(Description, levels=pathways))) + geom_point(aes(color = NES, size = size)) + 
  scale_color_stepsn(colours = c("#63bff0", "#a7d5ed", "#e2e2e2","#e2e2e2", "#e1a692", "#de6e56", "#e14b31", "#c23728"), breaks = c(-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5))+
  theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.title = element_blank())+
  scale_size_discrete(name = "p-value", labels = c("ns", "<0.05", "<0.005"))
pdf('fig3O.pdf', width = 6, height = 3.5)
print(p)
dev.off()


