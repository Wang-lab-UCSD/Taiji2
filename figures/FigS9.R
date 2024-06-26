library(ggplot2)
library(xlsx)
library(dplyr)
library(readr)
library(RColorBrewer)
library(data.table)
maindir <- '~/allcombn/'
setwd(maindir)
source(paste0(maindir, 'scripts/functions/getRegulatees.R'))
source(paste0(maindir, 'scripts/functions/utils.R'))
source(paste0(maindir, 'scripts/functions/firstUp.R'))
source(paste0(maindir,"scripts/functions/enrich_all.R"))

# import PageRank and group data
Data <- read.table("GeneRanks.tsv", check.names = F)
rownames(Data) <- firstup(tolower(rownames(Data)))
Data <- Data[!rownames(Data) %in% c("Dnmt1"),] # remove DNMT1

# group info
group_sorted <- read.table('group_file.txt',header = T, row.names = 1)
group_sorted$oldName <- rownames(group_sorted)
rownames(group_sorted) <- group_sorted$newName
names(Data) <- group_sorted$newName[match(names(Data),group_sorted$oldName)]

# remove two TRM.lung samples from Hayward dataset -------
group_sorted <- group_sorted[!grepl("TRM.*Hayward",rownames(group_sorted)),] # remove two TRM.lung samples from Hayward dataset 
# group_sorted <- group_sorted[!grepl("TRM.*Hayward|liver|TE_Kaech|MP_Kaech",rownames(group_sorted)),] # remove Kaech and Hayward


#### 1. import rna_seq ---------------
rna_seq <- read.table("expression_profile.tsv", check.names = F)
# row.names(rna_seq) <- toupper(rownames(rna_seq))
rna_seq <- rna_seq[complete.cases(rna_seq),]
names(rna_seq) <- group_sorted$newName[match(names(rna_seq),group_sorted$oldName)]
rna_seq <- rna_seq[,rownames(group_sorted)]
rna_seq <- rna_seq[rowSums(rna_seq)>0,]

## read TF list
TF <- toupper(sub('g','',readLines('fig3I_tfs.txt')))


## get regulatees
setwd(paste0(maindir,'network/'))
trm_samples <- group_sorted[group_sorted$Group=='TRM',] %>% pull(oldName)
network <- read_delim('allTFs_TRM_regulatees.tsv', delim = ' ')
lapply(TF, function(x) network %>% 
         dplyr::filter(`:START_ID`==toupper(x) & filename %in% trm_samples) %>% 
         write.table(file = paste0(x,"_TRM_regulatees.tsv"), quote = FALSE, row.names = FALSE))

## manually move the contents to a new folder
setwd(paste0(maindir, 'network/figS9_trm_regulatees/'))
genesets <- lapply(TF, function(x) getRegulatees('TRM', x, cut = 2000, input_dir = "./"))
genesets <- lapply(TF, function(x) getRegulatees('TRM', x, cut = 1000, input_dir = "./"))


## gsea analysis
library(clusterProfiler)
library(ReactomePA)
library(org.Mm.eg.db)
library(msigdbr)
source(paste0(maindir,"scripts/functions/gsea_c6_c7_hallmark_GO.R"))

top <- 2000
# label <- 'rna'
label <- 'edge_weight'

# enrichment analysis
# L <- list.files(path="./", pattern = paste0(".*cut",top,".*",label,".txt"))
# lapply(L, function(x) enrich_all(x=x, is.plot = F)) 
# source(paste0(maindir,"scripts/functions/gomf_gsea_c6_c7_hallmark.R"))
# lapply(L, function(x) enrich_all_v2(x=x, is.plot = F)) 

# GSEA analysis, which needs ranked gene list
L <- list.files(path="./", pattern = paste0(".*cut",top,".*",label,".tsv"))
lapply(L, function(x) {
  df <- read.delim(x, sep = ' ')
  out <- df %>% group_by(regulatees) %>% summarise(mean = mean(!!rlang::sym(names(df)[3]))) %>% arrange(-mean) 
  x1 <- sub('.tsv','_ranked_for_gsea.csv',x)
  write.csv(out, x1, row.names = F)
  gsea_c6_c7_hallmark_GO(x1, pvalue = 1)
})


### load selected pathways from scRNA-seq analysis
selected_pathways <- readLines(paste0(maindir,'fig3I_pathways.txt'))
# L <- list.files(path="./", pattern = paste0(".*cut",top,".*",label,"_enrich_all_v2.xlsx"))
L <- list.files(path="./", pattern = paste0(".*cut",top,".*",label,".*gsea_result.xlsx"))
df <- do.call('rbind', lapply(L, function(x) read.xlsx(x, sheetIndex = 1) %>% filter(Description %in% selected_pathways) %>% mutate(TF = firstup(tolower(gsub('_.*','',x))))))
write.xlsx(df, 'figS9_cut2000_ew_regulatees_gsea.xlsx', row.names = F)

### dot plot
df <- read.xlsx('figS9_cut2000_ew_regulatees_gsea.xlsx', sheetIndex = 1)

### transform data
df2 <- df %>% mutate(size = as.factor(ifelse(pvalue > 0.05, 1, ifelse(pvalue>0.005, 2, 3))))

## plot
l2 <- unique(df$Description)
selected_pathways_v2 <- c(intersect(selected_pathways, l2), setdiff(l2, selected_pathways))
TFs <- sub('g','',readLines(paste0(maindir,'fig3I_tfs.txt')))

p <- ggplot(df2, aes(x=factor(TF, levels = TFs), y=factor(Description, levels = selected_pathways_v2))) + geom_point(aes(color = NES, size = size)) + 
  scale_color_gradientn(colours = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100)) + 
  theme_classic() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.title = element_blank())+
  scale_size_discrete(name = "p.value", labels = c("ns", "<0.05", "<0.005"))
pdf('figS9.pdf', width = 9.5, height = 6.5)
print(p)
dev.off()


### plot the significant pathways
L <- list.files(path="./", pattern = paste0(".*cut",top,".*",label,".*gsea_result.xlsx"))
df <- do.call('rbind', lapply(L, function(x) read.xlsx(x, sheetIndex = 1) %>% filter(p.adjust <= 0.1) %>% mutate(TF = firstup(tolower(gsub('_.*','',x))))))
write.xlsx(df, 'figS9_cut2000_ew_regulatees_gsea.xlsx', row.names = F)

